#!/panvol1/simon/bin/python2.7

from __future__ import division

import argparse
import subprocess
import genobox_modules
import os
import subprocess
import sys
import re
import logging
import gzip

def get_genome(chr_file):
   '''Read chromosome file into list of list'''
   fh = open(chr_file, 'r')
   d = {}
   for line in fh:
      line = line.rstrip()
      s = line.split('\t')
      d[s[0]] = s
   return d


def filter_vcf(f, genome, rms_mq=25.0, Q=30.0, allelic_balance=0.2):
   '''Filter vcf from stdin for several expressions'''
   
   # open input file
   if f == 'stdin':
      fh = sys.stdin
   elif f.endswith('.gz'):
         fh = gzip.open(f, 'rb')
   else:
      fh = open(f, 'r')
      
   # open output processes
   
   f_base = os.path.split(f)[1]
   f_indels = 'tmp/tmp.%s.indels.vcf' % f_base
   f_snps = 'tmp/tmp.%s.snps.vcf' % f_base
   f_ref = 'tmp/tmp.%s.ref.vcf' % f_base
   
   fh_indels = open(f_indels, 'w')
   fh_snps = open(f_snps, 'w')
   fh_ref = open(f_ref, 'w')
   
   
   # read genome file and create regular exp
   genome_d = get_genome(genome)
   reg_dp = re.compile('DP=(\d+)')
   reg_mq = re.compile('MQ=(\d+)')
   reg_sb = re.compile('SB=(.?\d+\.\d+e.\d+)')
   reg_sb2 = re.compile('SB=(.?\d+\.\d+)')
   reg_rpos = re.compile('ReadPosRankSum=(.?\d+\.\d+)')
   
   for line in fh:
      # headers
      if line.startswith('#'):
         fh_indels.write(line)
         fh_snps.write(line)
         fh_ref.write(line)
         continue
      
      # check that chrom is in genome-file
      fields = line.split('\t')
      if genome_d.has_key(fields[0]):
         pass
      else:
         raise ValueError('%s is not in genome-file\n' % fields[0])
      
      # check for depth, mapping qual, qual filters
      info = fields[7]
      chr = fields[0]
      qual = float(fields[5])
      dp = int(reg_dp.search(info).group(1))
      mq = float(reg_mq.search(info).group(1))
      
      if dp < int(genome_d[chr][4]):
         continue
      elif dp > int(genome_d[chr][5]):
         continue
      elif qual < Q:
         continue
      elif mq < rms_mq:
         continue
      
      # now check for ref/variant
      if fields[4] == '.': fh_ref.write(line)
      else:
         # check for additional if variant
         
         # strand bias, ReadPosRankSum
         sb_m = reg_sb.search(info)
         rpos_m = reg_rpos.search(info)
         if not sb_m:
            sb_m2 = reg_sb2.search(info)
            if sb_m2:
               if float(sb_m2.group(1)) > 0.1:
                  continue
         else:
            if float(sb_m.group(1)) > 0.1:
               continue
         
         if rpos_m:
            if float(rpos_m.group(1)) <= -8.0:
               continue
         
         # heterozygote/homozygote
         tags = fields[9]
         tag_fields = tags.split(':')
         
         # allelic balance
         ab_split = tag_fields[1].split(',')
         ab_sum = sum(map(int, ab_split))          # may happen with deletions
         if ab_sum == 0:
            pass
         else:
            ab = float(ab_split[0])/(sum(map(int, ab_split)))
         
         if ab > 0.5: ab = 1-ab
         
         state = tag_fields[0]
         if state == '0/1':
            if genome_d[chr][3] == 'haploid':
               continue
            if ab < allelic_balance:
               continue
         #if state == '1/1':
         #   if 1-ab < allelic_balance:
         #      continue
         
         # write out
         if len(fields[3]) > 1:
            fh_indels.write(line)
         elif len(fields[4]) > 1:
            if fields[4].find(',') > -1:
               fh_snps.write(line)
            else:
               fh_indels.write(line)
         else:
            fh_snps.write(line)
   
   # end Popens
   fh_snps.close()
   fh_indels.close()
   fh_ref.close()
   
   return [f_indels, f_snps, f_ref]

def gatk_filtration(filtered_vcfs, fa, rmsk, prune, logger):
   '''Filter indels using VariantFiltration'''
   
   import genobox_modules
   import os
   
   paths = genobox_modules.setSystem()
   gatk_cmd = paths['GATK_home'] + 'GenomeAnalysisTK.jar'
   java_cmd = 'java -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms3000m -Xmx3000m -jar '
   cmd = java_cmd + gatk_cmd
   bed_cmd = '%sintersectBed' % paths['bedtools_home']
   ec = []
   
   # header file (because bedtools does not print header)
   header_vcf = filtered_vcfs[0].replace('.raw.vcf.gz.indels.vcf', '.header.vcf')
   grep_cmd = '''head -n 2000 %s | grep "#" > %s''' % (filtered_vcfs[0], header_vcf)
   ec.append(subprocess.call(grep_cmd, shell=True))
   
   # indels
   f_indels_pass = filtered_vcfs[0].replace('.raw.vcf.gz.indels.vcf', '.indels.pass.vcf')
   if rmsk:
      arg = """ -T VariantFiltration -R %s --variant %s --clusterSize 2 --clusterWindowSize %i -o /dev/stdout | perl -ne 'if ($_ =~ m/^INFO/ or $_ =~ m/^WARN/) {} else {print $_}' | cat - | %s -v -a stdin -b %s | grep "PASS" | cat %s - > %s""" % (fa, filtered_vcfs[0], prune, bed_cmd, rmsk, header_vcf, f_indels_pass)
   else:
      arg = ' -T VariantFiltration -R %s --variant %s --clusterSize 2 --clusterWindowSize %i -o %s' % (fa, filtered_vcfs[0], prune, f_indels_pass)
   #logger.info(cmd+arg)
   ec.append(subprocess.call(cmd+arg, shell=True))
   
   # snps
   # header file (because bedtools does not print header)
   header_vcf = filtered_vcfs[1].replace('.raw.vcf.gz.snps.vcf', '.header.vcf')
   grep_cmd = '''head -n 2000 %s | grep "#" > %s''' % (filtered_vcfs[0], header_vcf)
   ec.append(subprocess.call(grep_cmd, shell=True))
   
   f_snps_pass = filtered_vcfs[1].replace('.raw.vcf.gz.snps.vcf', '.snps.pass.vcf')
   if rmsk:
      arg = ''' -T VariantFiltration -R %s --variant %s --clusterSize 2 --clusterWindowSize %i --mask %s -maskName indel -o /dev/stdout | perl -ne 'if ($_ =~ m/^INFO/ or $_ =~ m/^WARN/) {} else {print $_}'  | cat - | %s -v -a stdin -b %s | grep "PASS" | cat %s - > %s''' % (fa, filtered_vcfs[1], prune, f_indels_pass, bed_cmd, rmsk, header_vcf, f_snps_pass)
   else:
      arg = ' -T VariantFiltration -R %s --variant %s --clusterSize 2 --clusterWindowSize %i --mask %s -maskName indel -o %s' % (fa, filtered_vcfs[1], prune, f_indels_pass, f_snps_pass)
   #logger.info(cmd+arg)
   ec.append(subprocess.call(cmd+arg, shell=True))
   
   # ref
   # header file (because bedtools does not print header)
   header_vcf = filtered_vcfs[2].replace('.raw.vcf.gz.ref.vcf', '.header.vcf')
   grep_cmd = '''head -n 2000 %s | grep "#" > %s''' % (filtered_vcfs[0], header_vcf)
   ec.append(subprocess.call(grep_cmd, shell=True))
   
   # dont grep PASS because it is not written for non-variants
   f_ref_pass = filtered_vcfs[2].replace('.raw.vcf.gz.ref.vcf', '.ref.pass.vcf')
   if rmsk:
      if f_ref_pass.find('.MT.') > -1:
         arg = ''' -T VariantFiltration -R %s --variant %s --mask %s -maskName indel -o /dev/stdout | perl -ne 'if ($_ =~ m/^INFO/ or $_ =~ m/^WARN/) {} else {print $_}' | cat - | %sgenobox_vcffilter_mt.py --rmsk %s | cat %s - > %s''' % (fa, filtered_vcfs[2], f_indels_pass, paths['genobox_home'], rmsk, header_vcf, f_ref_pass)
      else:
         arg = ''' -T VariantFiltration -R %s --variant %s --mask %s -maskName indel -o /dev/stdout | perl -ne 'if ($_ =~ m/^INFO/ or $_ =~ m/^WARN/) {} else {print $_}' | cat - | %s -v -a stdin -b %s | cat %s - > %s''' % (fa, filtered_vcfs[2], f_indels_pass, bed_cmd, rmsk, header_vcf, f_ref_pass)
   else:
      arg = ' -T VariantFiltration -R %s --variant %s --mask %s -maskName indel -o %s' % (fa, filtered_vcfs[2], f_indels_pass, f_ref_pass)
   #logger.info(cmd+arg)
   ec.append(subprocess.call(cmd+arg, shell=True))
   
   
   # remove non-pass lines and join to one file #
   if rmsk:
      final_vcf = filtered_vcfs[0].replace('tmp/tmp.', 'genotyping/').replace('.raw.vcf.gz.indels.vcf', '.filtered.rmsk.vcf.gz')
   else:
      final_vcf = filtered_vcfs[0].replace('tmp/tmp.', 'genotyping/').replace('.raw.vcf.gz.indels.vcf', '.filtered.vcf.gz')
   arg = ''' -T CombineVariants -U LENIENT_VCF_PROCESSING --assumeIdenticalSamples -R %s --variant %s --variant %s --variant %s -o /dev/stdout | perl -ne 'if ($_ =~ m/^INFO/ or $_ =~ m/^WARN/) {} else {print $_}' | gzip -c - > %s''' % (fa, f_indels_pass, f_snps_pass, f_ref_pass, final_vcf)
   #logger.info(cmd+arg)
   ec.append(subprocess.call(cmd+arg, shell=True))
   
   # copy snp and indel calls
   if rmsk:
      final_indel = filtered_vcfs[0].replace('tmp/tmp.', 'genotyping/').replace('.raw.vcf.gz.indels.vcf', '.filtered.rmsk.indels.vcf.gz')
   else:
      final_indel = filtered_vcfs[0].replace('tmp/tmp.', 'genotyping/').replace('.raw.vcf.gz.indels.vcf', '.filtered.indels.vcf.gz')
   arg = 'gzip -c %s > %s' % (f_indels_pass, final_indel)
   #logger.info(arg)
   ec.append(subprocess.call(arg, shell=True))
   
   if rmsk:
      final_snp = filtered_vcfs[1].replace('tmp/tmp.', 'genotyping/').replace('.raw.vcf.gz.snps.vcf', '.filtered.rmsk.snps.vcf.gz')
   else:
      final_snp = filtered_vcfs[1].replace('tmp/tmp.', 'genotyping/').replace('.raw.vcf.gz.snps.vcf', '.filtered.snps.vcf.gz')
   arg = 'gzip -c %s > %s' % (f_snps_pass, final_snp)
   #logger.info(arg)
   ec.append(subprocess.call(arg, shell=True))
   

def clean(f):
   '''Clean up tmp and raw files'''
   
   import genobox_modules
   import os
      
   # finding files to delete
   f_base = os.path.split(f)[1]
   f_base = f_base.replace('.raw.vcf.gz', '')
   f_base = 'tmp/tmp.' + f_base
   files_to_delete = []
   files_to_delete.append(f_base+'.header.vcf')
   files_to_delete.append(f_base+'.indels.pass.vcf')
   files_to_delete.append(f_base+'.indels.pass.vcf.idx')
   files_to_delete.append(f_base+'.raw.vcf.gz.indels.vcf')
   files_to_delete.append(f_base+'.raw.vcf.gz.indels.vcf.idx')
   files_to_delete.append(f_base+'.raw.vcf.gz.ref.vcf')
   files_to_delete.append(f_base+'.raw.vcf.gz.ref.vcf.idx')
   files_to_delete.append(f_base+'.raw.vcf.gz.snps.vcf')
   files_to_delete.append(f_base+'.raw.vcf.gz.snps.vcf.idx')
   files_to_delete.append(f_base+'.ref.pass.vcf')
   files_to_delete.append(f_base+'.ref.pass.vcf.idx')
   files_to_delete.append(f_base+'.snps.pass.vcf')
   files_to_delete.append(f_base+'.snps.pass.vcf.idx')
   
   # deleting files
   genobox_modules.rm_files(files_to_delete)


def main(args, logger):
   '''Filter raw vcfs'''
   
   filtered_vcfs = filter_vcf(args.vcf, args.genome, 25.0, args.Q, args.ab)
   gatk_filtration(filtered_vcfs, args.fa, args.rmsk, args.prune, logger)
   clean(args.vcf)


if __name__ == '__main__':
   
   # create the parser
   parser = argparse.ArgumentParser(description=
   '''Filter vcf created from gatk (hard-filter)''')
   
   # add the arguments
   parser.add_argument('--vcf', help='input vcf files')
   parser.add_argument('--fa', help='reference fasta')
   parser.add_argument('--genome', help='file containing chromosomes to analyse, format: chrom\tchrom_len\tchrom_short_name\ploidy\tlow_d\thigh_d', default=None)
   parser.add_argument('--Q', help='minimum quality score', type=float, default=20.0)
   parser.add_argument('--rmsk', help='rmsk to use', default=None)
   parser.add_argument('--ab', help='allelic balance threshold [0] (no filter)', type=float, default=0)
   parser.add_argument('--prune', help='distance (nt) to prune within [0] (no filter)', type=int, default=0)
   parser.add_argument('--log', help='log level [info]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args('--vcf genotyping/BI16.flt.sort.rmdup.realign.22.raw.vcf.gz --fa /panvol1/simon/projects/arctic/references/build37/hs.build37.1.fa --genome /panvol1/simon/projects/arctic/references/build37/hs.build37.1.fa.male.hc.genome --Q 40.000000 --rmsk /panvol1/simon/projects/arctic/references/build37/rmsk.build37 --ab 0.200000 --prune 5'.split())
   #args = parser.parse_args('--vcf /home/projects5/pr_99009/people/simon/projects/stanford_genomes/analysis/BI16/genotyping/BI16.flt.sort.rmdup.realign.MT.raw.vcf.gz --fa /panvol1/simon/projects/arctic/references/build37/hs.build37.1.fa --genome /panvol1/simon/projects/arctic/references/build37/hs.build37.1.fa.male.hc.genome --Q 40.000000 --rmsk /panvol1/simon/projects/arctic/references/build37/rmsk.build37 --ab 0.200000 --prune 5'.split())
   
   # set logging
   logger = logging.getLogger('genobox.py')
   hdlr = logging.FileHandler('genobox.log')
   formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
   hdlr.setFormatter(formatter)
   logger.addHandler(hdlr)
   if args.log == 'info':
      logger.setLevel(logging.INFO)
   
   # set paths
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   
   main(args, logger)
   
