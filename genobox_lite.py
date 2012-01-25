#!/panvol1/simon/bin/python

# index fasta
# run fastqc, parse seq to trim
# run trimming
# map trimmed reads using bwa
# get mapped q30 and sort
# rmdup
# merge to final
# genotype using samtools mpileup
# filter based on pruning
# filter based on ploidy

import argparse
import os
import subprocess

bwa = '/panvol1/simon/bin/bwa'
shrimp = '/panvol1/simon/bin/gmapper-cs'
fastqc = '/panvol1/simon/bin/fastqc'
genobox_trim_se = '/panvol1/simon/bin/genobox/genobox_trim_h.py'
genobox_trim_pe = '/panvol1/simon/bin/genobox/genobox_trim_pe.py'
fastx_detect_fq = '/panvol1/simon/bin/pipeline/fastx_detect_fq.py'
fq2allstd = '/panvol1/simon/bin/fq2allstd.pl'

## PIPELINE FUNCTIONS ##

def bwa_index(ref):
   '''Index reference fasta with bwa'''
   
   cmd = bwa + ' index %s' % ref
   exitcode = subprocess.call(cmd, shell=True)
   if exitcode: raise SystemExit('Failed at bwa_index: %s\n' % cmd)

def shrimp_index(ref, n):
   '''Index reference for colorspace'''
   
   cmd = shrimp + ' -N %i -S %s.cs %s' % (n, ref, ref)
   exitcode = subprocess.call(cmd, shell=True)
   if exitcode: raise SystemExit('Failed at shrimp_index: %s\n' % cmd)

def fastqc_run(se, pe1, pe2):
   '''Run fastqc on input reads'''
   
   if not os.path.exists('fastqc'):
      os.makedirs('fastqc')
   
   cmd = fastqc + ' -o fastqc %s %s %s'  (' '.join(se), ' '.join(pe1), ' '.join(pe2))
   exitcode = subprocess.call(cmd, shell=True)
   if exitcode: raise SystemExit('Failed at fastqc_run: %s\n' % cmd)
   
def fastqc_parse(se, pe1, pe2):
   '''Parse adaptors to clip in reads'''
   
   # parse out 
   return dict_of_files->list_of_adaptors

def illumina_trim(se, pe1, pe2, files_adaptors_dict, gz):
   '''Trim reads'''
      
   def single_trim(se, files_adaptors_dict, gz):
      '''Create single end trim calls'''
      
      calls = []
      outfiles_se = []
      for i,f in enumerate(se):
         if gz: outfile_se = 'trimmed/' + os.path.split(f)[1] + '.trim.fq.gz'
         else: outfile_se = 'trimmed/' + os.path.split(f)[1] + '.trim.fq'
         outfiles_se.append(outfile_se)
         cmd = genobox_trim_se + ' --i %s --min_baseq 15 --o %s --min_adaptor_match 15 --adaptors %s' % (f, outfile_se, ' '.join(files_adaptors_dict[f]))
         if gz: arg = arg + ' --gz' 
         calls.append(cmd+arg)
      return (calls, outfiles_se)
   
   def paired_trim(pe1, pe2, files_adaptors_dict, gz):
      '''Create paired end trim calls'''
      
      if len(args.pe1) != len(args.pe2): raise ValueError('same number of files must be given to --pe1 and --pe2')
      
      calls = []
      outfiles_pe1 = []
      outfiles_pe2 = []
      for i,f in enumerate(pe1):
         if gz:
            outfile_pe1 = 'trimmed/' + os.path.split(args.pe1[i])[1] + '.trim.fq.gz'
            outfile_pe2 = 'trimmed/' + os.path.split(args.pe2[i])[1] + '.trim.fq.gz'
         else:
            outfile_pe1 = 'trimmed/' + os.path.split(args.pe1[i])[1] + '.trim.fq'
            outfile_pe2 = 'trimmed/' + os.path.split(args.pe2[i])[1] + '.trim.fq'
         outfiles_pe1.append(outfile_pe1)
         outfiles_pe2.append(outfile_pe2)
         arg = ' --i %s %s --min_baseq 15 --min_adaptor_match 15 --adaptors %s' % (pe1[i], pe2[i], ' '.join(files_adaptors_dict[f]))
         if gz: arg = arg + ' --gz'
         calls.append(cmd+arg)
      return (calls, outfiles_pe1, outfiles_pe2)
   
   # main
   if not os.path.exists('trimmed'):
      os.makedirs('trimmed')
   
   se_calls, se_new = single_trim(se, files_adaptors_dict, gz)
   pe_calls, pe1_new, pe2_new = paired_trim(pe1, pe2, files_adaptors_dict, gz)
   
   calls = se_calls + pe_calls
   for c in calls:
      exitcode = subprocess.call(c, shell=True)
      if exitcode: raise SystemExit('Failed at illumina_trim: %s\n' % cmd)
   
   return (se_new, pe1_new, pe2_new)


def bwa_pe_aln(pe1, pe2, ref, n):
   '''BWA alignment of paired reads'''
   
   
   pe1_fq = map(setfqtype, pe1)
   pe2_fq = map(setfqtype, pe2)
   
   sys.stderr.write('Running alignments\n')
   cmds = []
   sais1 = []
   sais2 = []
   for i in range(len(pe1)):
      sai1 = os.path.split(pe1[i])[1] + '.sai'
      sai2 = os.path.split(pe2[i])[1] + '.sai'
      if pe1_fq[i][1] == 'Illumina':
         cmd1 = bwa + ' aln -t %i -I %s %s -f %s' % (n, ref, pe1_fq[i][0], sai1)
         cmd2 = bwa + ' aln -t %i -I %s %s -f %s' % (n, ref, pe2_fq[i][0], sai2)
      else:
         cmd1 = bwa + ' aln -t %i %s %s -f %s' % (n, ref, pe1_fq[i][0], sai1)
         cmd2 = bwa + ' aln -t %i %s %s -f %s' % (n, ref, pe2_fq[i][0], sai2)
      
      cmds.extend([cmd1, cmd2])
      sais1.extend([sai1])
      sais2.extend([sai2])   
      
   for c in cmds:
      sys.stderr.write(c+'\n')
      p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      stdout,stderr = p.communicate()
      #print stderr
      
   return (sais1, sais2)


def pe_sai2bam(sais1, sais2, pe1, pe2, ref):
   '''Convert sai to bam by sampe. Capture stderr to get isizes'''
   
   sys.stderr.write('Converting alignments\n')
   def average(values):
      '''Computes the arithmetic mean of a list of numbers.'''
      return sum(values, 0.0) / len(values)
   
   cmds = []
   bams = []
   for i in range(len(sais1)):
      bam = os.path.split(pe1[i])[1] + 'pe.q30.sort.bam'
      bams.append(bam)
      cmd = bwa + ' sampe %s %s %s %s %s | %s view -Sb -F 4 -q 30 - | samtools sort - > %s' % (ref, sais1[i], sais2[i], pe1[i], pe2[i], samtools, bam)
      cmds.append(cmd)
   
   stderr_L = []
   for c in cmds:
      sys.stderr.write(c+'\n')
      p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      stdout, stderr = p.communicate()
      stderr_L.append(stderr)
      
   return (bams)


def single_bwa(se, ref, n):
   '''Alignment of single end reads'''
   
   sys.stderr.write('Running single alignments\n')
   se_fq = map(setfqtype, se)

   cmds = []
   sais = []
   bams = []
   for i in range(len(se)):
      sai = os.path.split(se[i])[1] + '.sai'
      bam = os.path.split(se[i])[1] + '.se.q30.sort.bam'
      if se_fq[i][1] == 'Illumina': cmd1 = bwa + ' aln -I -t %i %s %s -f %s' % (n, contigs, se[i], sai)
      else: cmd1 = bwa + ' aln -t %i %s %s -f %s' % (n, contigs, se[i], sai)
      cmd2 = bwa + ' samse %s %s %s | %s view -Sb -F 4 -q 30 - | samtools sort - > %s' % (contigs, sai, se[i], samtools, bam)
      
      cmds.extend([cmd1, cmd2])
      sais.append(sai)
      bams.append(bam)
   
   for c in cmds:
      sys.stderr.write(c+'\n')
      p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      stdout,stderr = p.communicate()
      #print stderr
   
   return bams
      
   
   

def bwa_bwasw():
   '''Map long reads using bwasw'''

def shrimp_aln():
   '''Map cs reads using shrimp'''

def bamprocess():
   '''Filter alignment for Q30, sort, rmdup and merge'''
   
def samtools_genotype():
   '''Genotype using samtools'''

def vcf_filter():
   '''Filter vcf calls (quality, ploidy)'''

def vcf_prune():
   '''Prune vcf calls'''



## HELPER FUNCTIONS ##

def setfqtype(f):
   '''Detect fqtype and convert to Sanger if Solexa (bwa cant handle Solexa quals)'''
   
   if f.endswith('.gz'): cmd = 'gzip -dc %s | %s' % (f, fastx_detect_fq)
   else: cmd = '%s --i %s' % (fastx_detect_fq, f)
   fqtype = subprocess.check_output(cmd, shell=True).rstrip()
   
   if fqtype == 'Solexa':
      sys.stderr.write('converting %s from solexa to sanger\n' % f)
      if f.endswith('.gz'):
         f_sanger = os.path.split(f)[1]+'.sanger.fq.gz'
         cmd = 'gzip -dc %s | %s sol2std | gzip -c > %s' % (f, fq2allstd, f_sanger)
         ec = subprocess.call(cmd, shell=True)
      else:
         f_sanger = os.path.split(f)[1]+'.sanger.fq'
         cmd = 'cat %s | %s sol2std > %s' % (f, fq2allstd, f_sanger)
         ec = subprocess.call(cmd, shell=True)
      
      return (f_sanger, 'Sanger')
   else:
      return (f, fqtype)


## MAIN PROGRAM ##

def main():
   '''Main program'''
   
   # illumina mapping #
   # output q30 sorted bam #
   
   if args.se:
   fastqc
   # bwa se
   se_bams = single_bwa(se, args.ref, args.n)
   
   if args.pe1:
   # bwa pe
   sais1, sais2 = bwa_pe_aln(pe1, pe2, args.ref, args.n)
   pe_bams = pe_sai2bam(sais1, sais2, pe1, pe2, args.ref)
   
   
   # bamprocessing #

if __name__ == '__main__':
   
   # create the parser
   parser = argparse.ArgumentParser(prog='genobox_lite.py', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=130), usage='%(prog)s assembler [options]', description='''Trimming, Mapping, Genotyping and variant filtering''')
   
   parser.add_argument('--se', help='input single end reads(s)', nargs='+')
   parser.add_argument('--pe1', help='input pair1 read(s)', nargs='+', action=set_abspath())
   parser.add_argument('--pe2', help='input pair2 read(s)', nargs='+', action=set_abspath())
   parser.add_argument('--long', help='input long read(s) (454, Ion Torrent)', nargs='+', action=set_abspath())
   parser.add_argument('--cs', help='input color space read(s)', nargs='+', action=set_abspath())
   parser.add_argument('--ref', help='input reference', required=True, action=set_abspath())
   parser.add_argument('--sample', help='outdir', required=True)
   parser.add_argument('--n', help='number of cpus [4]', default=4, type=int)
   parser.add_argument('--gz', help='input files are gzipped [False]', default=False, action='store_true')
   
   args = parser.parse_args()
   #args = parser.parse_args(''.split())
   
   # create dir
   if not os.path.exists(args.sample):
      os.makedirs(args.sample)
   
   os.chdir(args.sample)
   
   main(args)
