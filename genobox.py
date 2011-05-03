#!/panvol1/simon/bin/python2.7

# Simon Rasmussen, CBS - DTU

# Genobox is a toolbox for running Illumina mapping and genotyping:

#  alignment         Single and paired end alignment using bwa on a reference sequence
#  bam processing    Sort, filter, merge to libraries, rmdup and final merge of bam-files
#  genotyping        Genotyping using samtools on bam-file
#  vcf-filtering     Filter variant genotyping (vcf) for read depth, rmsk (repeat masking), allelic balance, ploidy, proximity pruning
#  vcf-annotation    Annotate vcf file with dbSNP information
#  bcf2ref           Extract high confidence same-as-reference positions filtering read depth, rmsk (repeat masking) and annotate using dbSNP
#  abgv              Run alignment, bamprocess, genotyping, vcf-filter, vcf-annotation and bcf2ref on input fastqs

import argparse
import genobox_moab
import subprocess
import logging
import os


def set_abspath():
   '''Returns absolute path of file to argparse'''
   class SetAbspath(argparse.Action):
      def __call__(self, parser, args, filenames, option_string=None):
         import os
         if type(filenames) == str:
            f_abs = os.path.abspath(filenames)
            setattr(args, self.dest, f_abs)
         elif type(filenames) == list:
            new_list = []
            for f in filenames:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, filenames)
   return SetAbspath

def required_interval(nmin,nmax):
   '''Enforces input numbers (int/float) to be within min and max'''
   class RequiredInterval(argparse.Action):
      def __call__(self, parser, args, value, option_string=None):
         if not nmin<=value<=nmax:
            msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(
               f=self.dest,nmin=nmin,nmax=nmax)
            raise argparse.ArgumentTypeError(msg)
         setattr(args, self.dest, value)
   return RequiredInterval


parser = argparse.ArgumentParser(prog='genobox.py', description='''Genobox is a toolbox for mapping and genotyping using Illumina read data''', usage='%(prog)s module [options]')

# general (parent parser)
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--sample', help='name of run and output directory', default=None)
parent_parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
parent_parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cbs]', default='cbs')
parent_parser.add_argument('--log', help='log level [INFO]', default='info')

# create subparsers
subparsers = parser.add_subparsers(dest='module')

# alignment
parser_alignment = subparsers.add_parser('alignment', help='Single and paired end alignment using bwa', parents=[parent_parser])
parser_alignment.add_argument('--se', help='input single end fqfiles', nargs='+', action=set_abspath(), default=[])
parser_alignment.add_argument('--pe1', help='input paired end fqfiles PE1', nargs='+', action=set_abspath(), default=[])
parser_alignment.add_argument('--pe2', help='input paired end fqfiles PE2', nargs='+', action=set_abspath(), default=[])
parser_alignment.add_argument('--fa', help='input fasta to map against', action=set_abspath())
parser_alignment.add_argument('--libfile', help='input parameter file', action=set_abspath())
parser_alignment.add_argument('--libs', help='if --libfile is not given, library names for each input fq', nargs='+', default=['lib'])
parser_alignment.add_argument('--a', help='maximum insert size for bwa sampe (-a) [500]', default=500, type=int)
parser_alignment.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int, action=required_interval(0,1000))

# bam process
parser_bamprocess = subparsers.add_parser('bamprocess', help='Sort, filter, merge, rmdup and final merge of bams', parents=[parent_parser])
parser_bamprocess.add_argument('--libfile', help='input parameter file', action=set_abspath())
parser_bamprocess.add_argument('--bam', help='input bam-files', nargs='+', action=set_abspath())
parser_bamprocess.add_argument('--mapq', help='mapping quality threshold', type=int, nargs='+', default=[30])
parser_bamprocess.add_argument('--libs', help='library name for each bam file', nargs='+', default=['lib'])
parser_bamprocess.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/', action=set_abspath())
parser_bamprocess.add_argument('--outbam', help='output file [alignment/final.flt.sort.rmdup.bam]', default='alignment/final.flt.sort.rmdup.bam')

# bam stats
parser_bamstats = subparsers.add_parser('bamstats', help='Statistics of alignment', parents=[parent_parser])
parser_bamstats.add_argument('--bam', help='input bam-file', action=set_abspath())
parser_bamstats.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath())

# genotyping
parser_genotyping = subparsers.add_parser('genotyping', help='Genotyping using samtools on bam-file', parents=[parent_parser])
parser_genotyping.add_argument('--bam', help='input file', action=set_abspath())
parser_genotyping.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath())
parser_genotyping.add_argument('--fa', help='reference fasta', action=set_abspath())
parser_genotyping.add_argument('--prior', help='type of prior to use for bayesian model (full, flat, cond2) [flat]', default='flat')
parser_genotyping.add_argument('--pp', help='posterior probability cutoff [0.001]', default=0.001, type=float, action=required_interval(0,1))
parser_genotyping.add_argument('--o', help='output bcffile [genotyping/snpcalls.bcf]', default='genotyping/snpcalls.bcf')

# vcffilter
parser_vcffilter = subparsers.add_parser('vcffilter', help='Filter variant genotyping (vcf)', parents=[parent_parser])
parser_vcffilter.add_argument('--bcf', help='input bcf var file', action=set_abspath())
parser_vcffilter.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath())
parser_vcffilter.add_argument('--caller', help='samtools or gatk [samtools]', default='samtools')
parser_vcffilter.add_argument('--Q', help='minimum quality score', type=float, default=20.0, action=required_interval(0,10000))
parser_vcffilter.add_argument('--ex', help='exhange chromosome names using file [None]', default=None)
parser_vcffilter.add_argument('--rmsk', help='rmsk to use', default=None, action=set_abspath())
parser_vcffilter.add_argument('--ab', help='allelic balance threshold [0.5] (no filter)', type=float, default=0.50, action=required_interval(0,1))
parser_vcffilter.add_argument('--prune', help='distance (nt) to prune within [0] (no filter)', type=int, default=0, action=required_interval(0,100000000))
parser_vcffilter.add_argument('--o', help='output vcffile [snpcalls.vcf.gz]', default='snpcalls.vcf.gz')

# dbsnp
parser_dbsnp = subparsers.add_parser('dbsnp', help='Annotate vcf file with dbSNP information', parents=[parent_parser])
parser_dbsnp.add_argument('--vcf', help='input vcf file', action=set_abspath())
parser_dbsnp.add_argument('--ex', help='exhange chromosome names using file [None]', default=None, action=set_abspath())
parser_dbsnp.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None, action=set_abspath())
parser_dbsnp.add_argument('--o', help='output vcf.gz [snpcalls.dbsnp.vcf.gz]', default='snpcalls.vcf.gz')

# bcf2ref
parser_bcf2ref = subparsers.add_parser('bcf2ref', help='Extract high confidence same-as-reference positions', parents=[parent_parser])
parser_bcf2ref.add_argument('--bcf', help='input bcf var file', action=set_abspath())
parser_bcf2ref.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--Q', help='minimum quality score', type=float, default=20.0, action=required_interval(0,10000))
parser_bcf2ref.add_argument('--ex', help='exhange chromosome names using file [None]', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--rmsk', help='rmsk to use', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--indels', help='indels vcf to remove', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--o', help='output prefix for ref.ann.vcf.gz files [refcalls]', default='refcalls')


# abgv
parser_abgv = subparsers.add_parser('abgv', help='Perform alignment, bam-processing, genotyping and filtering', parents=[parent_parser])
parser_abgv.add_argument('--se', help='input single end fqfiles', nargs='+', action=set_abspath(), default=[])
parser_abgv.add_argument('--pe1', help='input paired end fqfiles PE1', nargs='+', action=set_abspath(), default=[])
parser_abgv.add_argument('--pe2', help='input paired end fqfiles PE2', nargs='+', action=set_abspath(), default=[])
parser_abgv.add_argument('--fa', help='bwa index to map against', action=set_abspath())
parser_abgv.add_argument('--libfile', help='input parameter file', action=set_abspath())
parser_abgv.add_argument('--a', help='maximum insert size for bwa sampe (-a) [500]', default=500, type=int)
parser_abgv.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int, action=required_interval(0,1000))
parser_abgv.add_argument('--mapq', help='mapping quality threshold', type=int, nargs='+', default=[30])
parser_abgv.add_argument('--libs', help='if --libfile is not given, library names for each input fq', nargs='+', default=['lib'])
parser_abgv.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/', action=set_abspath())
parser_abgv.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath())
parser_abgv.add_argument('--prior', help='type of prior to use for bayesian model (full, flat, cond2) [flat]', default='flat')
parser_abgv.add_argument('--pp', help='posterior probability cutoff [0.001]', default=0.001, type=float, action=required_interval(0,1))
parser_abgv.add_argument('--caller', help='samtools or gatk [samtools]', default='samtools')
parser_abgv.add_argument('--Q', help='minimum quality score', type=float, default=20.0, action=required_interval(0,10000))
parser_abgv.add_argument('--rmsk', help='rmsk to use', default=None, action=set_abspath())
parser_abgv.add_argument('--ab', help='allelic balance threshold [0.5] (no filter)', type=float, default=0.50, action=required_interval(0,1))
parser_abgv.add_argument('--prune', help='distance (nt) to prune within [0] (no filter)', type=int, default=0, action=required_interval(0,100000000))
parser_abgv.add_argument('--ex', help='exhange chromosome names using file [None]', default=None, action=set_abspath())
parser_abgv.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None, action=set_abspath())
parser_abgv.add_argument('--indels', help='indels vcf to remove', default=None, action=set_abspath())
parser_abgv.add_argument('--ovar', help='output variant vcf-file [genotyping/\"name\".var.flt.vcf.gz]', default='genotyping/var.flt.vcf.gz')
parser_abgv.add_argument('--oref', help='output reference vcf-file [genotyping/\"name\".ref.flt.vcf.gz]', default='genotyping/ref.flt.vcf.gz')



# parse args
args = parser.parse_args()
#args = parser.parse_args('alignment --se SRR002081se.recal.fastq SRR002082se.recal.fastq --pe1 SRR002137pe_1.recal.fastq SRR002138pe_1.recal.fastq --pe2 SRR002137pe_2.recal.fastq SRR002138pe_2.recal.fastq --n 16 --sample NA12891 --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_rCRS.fa'.split(' '))
#args = parser.parse_args('alignment --se Kleb-10-213361_2_1_sequence.trim.fq --fa kleb_pneu.fa args.sample kleb_10_213361se --n 16'.split(' '))
#args = parser.parse_args('bamprocess --bam alignment/SRR002081se.recal.fastq.bam alignment/SRR002082se.recal.fastq.bam alignment/SRR002137pe_1.recal.fastq.bam alignment/SRR002138pe_1.recal.fastq.bam --mapq 30 --libs libA'.split(' '))
#args = parser.parse_args('bamprocess --bam alignment/Kleb-10-213361_2_1_sequence.trim.fq.bam --mapq 30 --outbam alignment/kleb_10_213361.flt.sort.rmdup.bam --n 16'.split(' '))
#args = parser.parse_args('bamstats --bam alignment/kleb_10_213361.flt.sort.rmdup.bam --genome ../kleb_pneu.genome'.split(' '))
#args = parser.parse_args('genotyping --bam alignment/kleb_10_213361.flt.sort.rmdup.bam --genome kleb_pneu.genome --fa ../kleb_pneu.fa --prior full --pp 0.01 --var --o genotyping/kleb_10_213361.bcf'.split(' '))
#args = parser.parse_args('vcffilter --bcf genotyping/kleb_10_213361.bcf --Q 40 --prune 5 --ab 0.20 --genome kleb_pneu.genome --o genotyping/kleb_10_213361.vcf'.split(' '))
#args = parser.parse_args('abgv --pe1 Kleb-10-213361_2_1_sequence.trim.fq --pe2 Kleb-10-213361_2_2_sequence.trim.fq --fa kleb_pneu.fa --genome kleb_pneu.genome --ab 0.2 --prune 5 args.sample kleb_10_213361'.split(' '))
#args = parser.parse_args('abgv --se data/*.trim --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa --r @RG\\tID:AA\\tSM:AusAboriginal\\tLB:ACD\\tPL:ILLUMINA\\tCN:BGI --lib_file libs --genome build37_rCRS.genome --prior flat --pp 0.0001 --Q 40 --rmsk rmsk_build37_rCRS.number.sort.genome --ex gi2number.build37_rCRS --ab 0.2 --prune 5 --n 16 --ovar genotyping/AusAboriginal.var.flt.vcf.gz --oref genotyping/AusAboriginal.ref.flt.vcf.gz'.split())

# If working dir is given, create and move to working directory else run where program is invoked
if args.sample:
   if not os.path.exists(args.sample):
      os.makedirs(args.sample)
   os.chdir(args.sample)
else:
   pass

# create log dir
if not os.path.exists('log'):
   os.makedirs('log')

# start logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# start modules
if args.module == 'alignment':
   from genobox_alignment import *
   bamfiles = start_alignment(args, logger)

elif args.module == 'bamprocess':
   from genobox_bamprocess import *
   final_bam = start_bamprocess(args.lib_file, args.bam, args.mapq, args.libs, args.tmpdir, args.queue, args.outbam, args.sample, logger)

elif args.module == 'bamstats':
   from genobox_bamstats import *
   start_bamstats(args, logger)

elif args.module == 'genotyping':
   from genobox_genotyping import *
   final_bcf = start_genotyping(args.bam, args.genome, args.fa, args.prior, args.pp, args.queue, args.o, logger)

elif args.module == 'vcffilter':
   from genobox_vcffilter import *
   final_vcf = start_vcffilter(args.bcf, args.genome, args.caller, args.Q, args.ex, args.rmsk, args.ab, args.prune, args.o, args.queue, args.sample, logger)

elif args.module == 'dbsnp':
   from genobox_dbsnp import *
   start_dbsnp(args.vcf, args.ex, args.dbsnp, args.o, args.queue, logger)

elif args.module == 'bcf2ref':
   from genobox_bcf2ref import *
   start_bcf2ref(args.bcf, args.genome, args.Q, args.ex, args.dbsnp, args.rmsk, args.indels, args.o, args.queue, args.sample, logger)

elif args.module == 'abgv':
   from genobox_abgv import *
   start_abgv(args, logger)
   
else:
   pass


## testing using kleb data ##

# echo "gi|206564770|gb|CP000964.1|\t5641239\tchr\thaploid\t10\t100" > kleb_pneu.genome
# genobox.py abgv --pe1 Kleb-10-213361_2_1_sequence.trim.fq --pe2 Kleb-10-213361_2_2_sequence.trim.fq --n 16 --sample Kleb-10-213361 --mapq 30 30 --libs A --genome kleb_pneu.genome --fa kleb_pneu.fa --ab 0.2 --prune 5


# test on 1000G data
# genobox.py alignment --se SRR002081se.recal.fastq SRR002082se.recal.fastq --pe1 SRR002137pe_1.recal.fastq SRR002138pe_1.recal.fastq --pe2 SRR002137pe_2.recal.fastq SRR002138pe_2.recal.fastq --n 16 --sample NA12891 --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa

# test on Aus Aboriginal
# create libs file

#echo 'ID\tData\tMAPQ\tLB\tCN\tPL\tSM' > libs.AusAboriginal.txt
#ll data/*.trim | perl -ane 'BEGIN{$count=0} ; if ($F[8] =~ m/.*\/(\w)/) {$count++ ; $id ="AusAboriginal_$count" ; print "$id\t/panvol1/simon/projects/AusAboriginal/$F[8]\t30\tlib$1\tBGI\tILLUMINA\tAusAboriginal\n" } ' >> libs.AusAboriginal.txt

# genobox.py abgv --se data/*.trim --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa  --libfile libs.AusAboriginal.txt --genome build37_rCRS.genome --prior flat --pp 0.0001 --Q 40 --rmsk rmsk_build37_rCRS.number.sort.genome --ex gi2number.build37_rCRS --ab 0.2 --prune 5 --n 4 --ovar genotyping/AusAboriginal.var.flt.vcf.gz --oref genotyping/AusAboriginal.ref.flt.vcf.gz



