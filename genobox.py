#!/panvol1/simon/bin/python2.7

# Simon Rasmussen, CBS - DTU

# Genobox is a toolbox for running Illumina mapping and genotyping:

#  alignment         Single and paired end alignment using bwa on a reference sequence
#  bam processing    Sort, filter, merge to libraries, rmdup and final merge of bam-files
#  genotyping        Genotyping using samtools on bam-file
#  vcf-filtering     Filter variant genotyping (vcf) for read depth, rmsk (repeat masking), allelic balance, ploidy, proximity pruning
#  vcf-annotation    Annotate vcf file with dbSNP information
#  bcf2ref           Extract high confidence same-as-reference positions filtering read depth, rmsk (repeat masking) and annotate using dbSNP


import argparse
import pipelinemod
import subprocess
import logging

parser = argparse.ArgumentParser(prog='genobox', description='''Genobox is a toolbox for running Illumina mapping and genotyping''')

# general
parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cbs]', default='cbs')
parser.add_argument('--log', help='log level [INFO]', default='info')

# create subparsers
subparsers = parser.add_subparsers(help='sub-command help', dest='module')

# alignment
parser_alignment = subparsers.add_parser('alignment', help='Single and paired end alignment using bwa')
parser_alignment.add_argument('--se', help='input single end fqfiles', nargs='+')
parser_alignment.add_argument('--pe1', help='input paired end fqfiles PE1', nargs='+')
parser_alignment.add_argument('--pe2', help='input paired end fqfiles PE2', nargs='+')
parser_alignment.add_argument('--bwaindex', help='bwa indexes to map against')
parser_alignment.add_argument('--a', help='maximum insert size for bwa sampe (-a) [500]', default=500, type=int)
parser_alignment.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int)

# bam process
parser_bamprocess = subparsers.add_parser('bamprocess', help='Sort, filter, merge to libraries, rmdup and final merge of bam-files')
parser_bamprocess.add_argument('--lib_file', help='input file with libraries')
parser_bamprocess.add_argument('--bam', help='input bam-files', nargs='+')
parser_bamprocess.add_argument('--mapq', help='mapping quality threshold', type=int, nargs='+')
parser_bamprocess.add_argument('--libs', help='library name for each bam file', nargs='+')
parser_bamprocess.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/')
parser_bamprocess.add_argument('--outbam', help='output file [alignment/final.flt.sort.rmdup.bam]', default='alignment/final.flt.sort.rmdup.bam')

# genotyping
parser_genotyping = subparsers.add_parser('genotyping', help='Genotyping using samtools on bam-file')
parser_genotyping.add_argument('--bam', help='input file')
parser_genotyping.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None)
parser_genotyping.add_argument('--fa', help='reference fasta')
parser_genotyping.add_argument('--prior', help='type of prior to use for bayesian model (full, flat, cond2) [flat]', default='flat')
parser_genotyping.add_argument('--pp', help='posterior probability cutoff [0.001]', default=0.001, type=float)
parser_genotyping.add_argument('--var', help='variants only, default all [False]', default=False, action='store_true')
parser_genotyping.add_argument('--o', help='output bcffile [genotyping/snpcalls.bcf]', default='genotyping/snpcalls.bcf')

# vcffilter
parser_vcffilter = subparsers.add_parser('vcffilter', help='Filter variant genotyping (vcf)')
parser_vcffilter.add_argument('--bcf', help='input bcf var file')
parser_vcffilter.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None)
parser_vcffilter.add_argument('--caller', help='samtools or gatk [samtools]', default='samtools')
parser_vcffilter.add_argument('--Q', help='minimum quality score', type=float, default=20.0)
parser_vcffilter.add_argument('--rmsk', help='rmsk to use', default=None)
parser_vcffilter.add_argument('--ab', help='allelic balance threshold [0.5] (no filter)', type=float, default=0.50)
parser_vcffilter.add_argument('--prune', help='distance (nt) to prune within [0] (no filter)', type=int, default=0)
parser_vcffilter.add_argument('--o', help='output vcffile [snpcalls.vcf]', default='snpcalls.vcf')

# dbsnp
parser_dbsnp = subparsers.add_parser('dbsnp', help='Annotate vcf file with dbSNP information')
parser_dbsnp.add_argument('--vcf', help='input vcf file')
parser_dbsnp.add_argument('--ex', help='exhange chromosome names using file [None]', default=None)
parser_dbsnp.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None)
parser_dbsnp.add_argument('--o', help='output vcf.gz [snpcalls.dbsnp.vcf.gz]', default='snpcalls.vcf.gz')

# bcf2ref
parser_bcf2ref = subparsers.add_parser('bcf2ref', help='Extract high confidence same-as-reference positions')
parser_bcf2ref.add_argument('--bcf', help='input bcf var file')
parser_bcf2ref.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None)
parser_bcf2ref.add_argument('--Q', help='minimum quality score', type=float, default=20.0)
parser_bcf2ref.add_argument('--ex', help='exhange chromosome names using file [None]', default=None)
parser_bcf2ref.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None)
parser_bcf2ref.add_argument('--rmsk', help='rmsk to use', default=None)
parser_bcf2ref.add_argument('--indels', help='indels vcf to remove', default=None)
parser_bcf2ref.add_argument('--o', help='output prefix for ref.ann.vcf.gz files [refcalls]', default='refcalls')



# parse args
args = parser.parse_args()
#args = parser.parse_args('--alignment --se SRR002081se.recal.fastq SRR002082se.recal.fastq --pe1 SRR002137pe_1.recal.fastq SRR002138pe_1.recal.fastq --pe2 SRR002137pe_2.recal.fastq SRR002138pe_2.recal.fastq --bwaindex /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa'.split(' '))
#args = parser.parse_args('--alignment --pe1 ../Kleb-10-213361_2_1_sequence.trim.fq --pe2 ../Kleb-10-213361_2_2_sequence.trim.fq --bwaindex kleb_pneu.fa'.split(' '))
#args = parser.parse_args('bamprocess --bam alignment/SRR002081se.recal.fastq.bam alignment/SRR002082se.recal.fastq.bam alignment/SRR002137pe_1.recal.fastq.bam alignment/SRR002138pe_1.recal.fastq.bam --mapq 30 --libs libA'.split(' '))

# start logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)


# set modules

if args.module == 'alignment':
   from genobox_alignment import *
   start_alignment(args.se, args.pe1, args.pe2, args.bwaindex, 'alignment/', args.qtrim, args.a, args.n, args.queue, logger)

elif args.module == 'bamprocess':
   from genobox_bamprocess import *
   start_bamprocess(args.lib_file, args.bam, args.mapq, args.libs, args.tmpdir, args.queue, args.outbam, logger)

elif args.module == 'genotyping':
   from genobox_genotyping import *
   start_genotyping(args.bam, args.genome, args.fa, args.prior, args.pp, args.var, args.queue, args.o, logger)

elif args.module == 'vcffilter':
   from genobox_vcffilter import *
   start_vcffilter(args.bcf, args.genome, args.caller, args.Q, args.rmsk, args.ab, args.prune, args.o, args.queue, logger)

elif args.module == 'dbsnp':
   from genobox_dbsnp import *
   start_dbsnp(args.vcf, args.ex, args.dbsnp, args.o, args.queue, logger)

elif args.module == 'bcf2ref':
   from genobox_bcf2ref import *
   start_bcf2ref(args.bcf, args.genome, args.Q, args.ex, args.dbsnp, args.rmsk, args.indels, args.o, args.queue, logger)

elif args.module == 'all':
   # start all
   pass


# genobox.py bcf2ref --bcf genotyping/abcalls.all.bcf --genome build37_rCRS.genome --Q 40 --ex /panvol1/simon/databases/hs_ref37_rCRS/gi2number.build37_rCRS --dbsnp /panvol1/simon/databases/dbsnp/dbsnp132_hg19.vcf.gz --rmsk /panvol1/simon/databases/hs_ref37_rCRS/rmsk/rmsk_build37_rCRS.number.sort.genome --header genotyping/header.vcf --indels genotyping/indels_for_filtering.number.vcf --o abcalls

# genobox.py bcf2ref --bcf genotyping/abcalls.all.bcf --genome tmp.genome --Q 40 --ex /panvol1/simon/databases/hs_ref37_rCRS/gi2number.build37_rCRS --dbsnp /panvol1/simon/databases/dbsnp/dbsnp132_hg19.vcf.gz --rmsk /panvol1/simon/databases/hs_ref37_rCRS/rmsk/rmsk_build37_rCRS.number.sort.genome --indels genotyping/indels_for_filtering.number.vcf --o abcalls
