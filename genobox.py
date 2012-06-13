#!/panvol1/simon/bin/python2.7

# Simon Rasmussen, CBS - DTU

# Genobox is a toolbox for running Illumina mapping and genotyping:

#  trim              Trim and filter single or paired end reads 
#  alignment         Single and paired end alignment using bwa on a reference sequence
#  bam processing    Sort, filter, merge to libraries, rmdup and final merge of bam-files
#  genotyping        Genotyping using samtools on bam-file
#  vcf-filtering     Filter variant genotyping (vcf) for read depth, rmsk (repeat masking), allelic balance, ploidy, proximity pruning
#  vcf-annotation    Annotate vcf file with dbSNP information
#  bcf2ref           Extract high confidence same-as-reference positions filtering read depth, rmsk (repeat masking) and annotate using dbSNP
#  abgv              Run alignment, bamprocess, genotyping, vcf-filter, vcf-annotation and bcf2ref on input fastqs

import argparse
import genobox_modules
import subprocess
import logging
import os

os.environ['PYTHONPATH'] = '/lib/python:/home/people/simon/lib/misc:/panvol1/simon/lib/python:/panvol1/simon/bin/genobox'
os.environ['LD_LIBRARY_PATH'] = '/panvol1/simon/lib'

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

def required_choices(choices):
   '''Enforces input numbers (int/float) to be within min and max'''
   class RequiredChoices(argparse.Action):
      def __call__(self, parser, args, input, option_string=None):
         for _in in input:
            if _in in choices:
               pass
            else:
               msg='argument "{f}" requires values of following: "{choices}"'.format(
               f=self.dest, choices=', '.join(choices))
               raise argparse.ArgumentTypeError(msg)
         setattr(args, self.dest, input)
   return RequiredChoices

def required_nargs_abspath(min,max):
   '''Enforces input to nargs to be between min and max long and sets abspath for input 1:'''
   class RequiredInterval(argparse.Action):
      def __call__(self, parser, args, value, option_string=None):
         if not min<=len(value)<=max:
            msg='argument "{f}" requires between {min} and {max} arguments'.format(
               f=self.dest,min=min,max=max)
            raise argparse.ArgumentTypeError(msg)
         setattr(args, self.dest, value)
         import os
         if type(value) == str:
            f_abs = os.path.abspath(value)
            setattr(args, self.dest, f_abs)
         elif type(value) == list:
            new_list = [value[0]]
            for f in value[1:]:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, value)
   return RequiredInterval


parser = argparse.ArgumentParser(prog='genobox.py',
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=25, width=140),
                                 description='''Genobox is a toolbox for mapping and genotyping using Illumina read data''', 
                                 usage='%(prog)s module [options]')

# general (parent parser)
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('--sample', help='name of run and output directory', default=None)
parent_parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
parent_parser.add_argument('--m', help='memory needed for high-memory jobs [7gb]', default='7gb')
parent_parser.add_argument('--gz', help='input files are gzipped', default=False, action='store_true')
parent_parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cbs]', default='cbs')
parent_parser.add_argument('--partition', help='partition to run on (ALL, cge-cluster, uv) [ALL]', default='ALL')
parent_parser.add_argument('--log', help='log level [INFO]', default='info')

# create subparsers
subparsers = parser.add_subparsers(dest='module')

# trim
parser_trim = subparsers.add_parser('trim', help='Trim and filter single or paired end reads', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py trim [options]')
parser_trim.add_argument('--se', help='input single end fqfiles', nargs='+', action=set_abspath(), default=[])
parser_trim.add_argument('--pe1', help='input paired end fqfiles PE1', nargs='+', action=set_abspath(), default=[])
parser_trim.add_argument('--pe2', help='input paired end fqfiles PE2', nargs='+', action=set_abspath(), default=[])
parser_trim.add_argument('--min_length', help='minimum length of a read to keep pairs [25]', type=int, default=25)
parser_trim.add_argument('--min_baseq', help='chomp bases with quality less than [20]', default=20, type=int)
parser_trim.add_argument('--min_avgq', help='minimum average quality of read [20]', default=20, type=int)
parser_trim.add_argument('--adaptors', help='adaptor sequence to clip [Illumina adaptors]', nargs='+', default=['GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT', 'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT', 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'])
parser_trim.add_argument('--keep_n', help='do not remove sequences containing N', default=False, action='store_true')
parser_trim.add_argument('--min_adaptor_match', help='minimum length of match to adaptor (0=all of adaptor) [20]', default=20, type=int)

# alignment
parser_alignment = subparsers.add_parser('alignment', help='Single and paired end alignment using bwa', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py alignment [options]')
parser_alignment.add_argument('--se', help='input single end fqfiles', nargs='+', action=set_abspath(), default=[])
parser_alignment.add_argument('--pe1', help='input paired end fqfiles PE1', nargs='+', action=set_abspath(), default=[])
parser_alignment.add_argument('--pe2', help='input paired end fqfiles PE2', nargs='+', action=set_abspath(), default=[])
parser_alignment.add_argument('--fa', help='input fasta to map against', required=True, action=set_abspath())
parser_alignment.add_argument('--libfile', help='input parameter file', action=set_abspath())
parser_alignment.add_argument('--libs', help='if --libfile is not given, library names for each input fq (order: se, pe1, pe2)', nargs='+', default=['lib'])
parser_alignment.add_argument('--pl', help='if --libfile is not given, platform for each input fq (order: se, pe1, pe2)', nargs='+', default=['ILLUMINA'], action=required_choices(['CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'PACBIO']))
parser_alignment.add_argument('--a', help='maximum insert size for bwa sampe (-a) [500]', default=500, type=int)
parser_alignment.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int, action=required_interval(0,1000))
parser_alignment.add_argument('--N', help='maximum number of alignments to output in the XA tag [3]', default=3, type=int)
parser_alignment.add_argument('--add_aln', help='additional parameters to bwa aln', default=None)
parser_alignment.add_argument('--bwa6', help='use bwa-6.1 (False)', default=False, action='store_true')


# bam process
parser_bamprocess = subparsers.add_parser('bamprocess', help='Sort, filter, merge, rmdup and final merge of bams', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py bamprocess [options]')
parser_bamprocess.add_argument('--libfile', help='input parameter file', action=set_abspath())
parser_bamprocess.add_argument('--bam', help='input bam-files', nargs='+', action=set_abspath())
parser_bamprocess.add_argument('--mapq', help='mapping quality threshold', type=int, nargs='+', default=[30])
parser_bamprocess.add_argument('--libs', help='library name for each bam file', nargs='+', default=['lib'])
parser_bamprocess.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/', action=set_abspath())
parser_bamprocess.add_argument('--outbam', help='output file [alignment/final.flt.sort.rmdup.bam]', default='alignment/final.flt.sort.rmdup.bam')

# bam stats
parser_bamstats = subparsers.add_parser('bamstats', help='Statistics of alignment', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py bamstats [options]')
parser_bamstats.add_argument('--bam', help='input bam-file', action=set_abspath())
parser_bamstats.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath())

# genotyping
parser_genotyping = subparsers.add_parser('genotyping', help='Genotyping using samtools on bam-file', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=100), usage='genobox.py genotyping [options]')
parser_genotyping.add_argument('--bam', help='input file', action=set_abspath())
parser_genotyping.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath(), required=True)
parser_genotyping.add_argument('--fa', help='reference fasta', action=set_abspath())
parser_genotyping.add_argument('--prior', help='type of prior to use for bayesian model (full, flat, cond2) [flat]', default='flat')
parser_genotyping.add_argument('--pp', help='posterior probability cutoff [0.001]', default=0.001, type=float, action=required_interval(0,1))
parser_genotyping.add_argument('--o', help='output bcffile [genotyping/snpcalls.bcf]', default='genotyping/snpcalls.bcf')

# vcffilter
parser_vcffilter = subparsers.add_parser('vcffilter', help='Filter variant genotyping (vcf)', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py vcffilter [options]')
parser_vcffilter.add_argument('--bcf', help='input bcf var file', action=set_abspath())
parser_vcffilter.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath())
parser_vcffilter.add_argument('--caller', help='samtools or gatk [samtools]', default='samtools')
parser_vcffilter.add_argument('--Q', help='minimum quality score', type=float, default=20.0, action=required_interval(0,10000))
parser_vcffilter.add_argument('--ex', help='exhange chromosome names using file [None]', default=None)
parser_vcffilter.add_argument('--rmsk', help='rmsk to use', default=None, action=set_abspath())
parser_vcffilter.add_argument('--ab', help='allelic balance threshold [0.5] (no filter)', type=float, default=0.50, action=required_interval(0,1))
parser_vcffilter.add_argument('--prune', help='distance (nt) to prune within [0] (no filter)', type=int, default=0, action=required_interval(0,100000000))
parser_vcffilter.add_argument('--o', help='output vcf.gz [genotyping/snpcalls.vcf.gz]', default='genotyping/snpcalls.vcf.gz')

# dbsnp
parser_dbsnp = subparsers.add_parser('dbsnp', help='Annotate vcf file with dbSNP information', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py dbsnp [options]')
parser_dbsnp.add_argument('--vcf', help='input vcf or vcf.gz file', action=set_abspath())
parser_dbsnp.add_argument('--ex', help='exhange chromosome names using file [None]', default=None, action=set_abspath())
parser_dbsnp.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None, action=set_abspath())
parser_dbsnp.add_argument('--o', help='output vcf.gz [genotyping/snpcalls.dbsnp.vcf.gz]', default='genotyping/snpcalls.dbsnp.vcf.gz')

# bcf2ref
parser_bcf2ref = subparsers.add_parser('bcf2ref', help='Extract high confidence same-as-reference positions', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=100), usage='genobox.py bcf2ref [options]')
parser_bcf2ref.add_argument('--bcf', help='input bcf file', action=set_abspath())
parser_bcf2ref.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath(), required=True)
parser_bcf2ref.add_argument('--Q', help='minimum quality score', type=float, default=20.0, action=required_interval(0,10000))
parser_bcf2ref.add_argument('--ex', help='exhange chromosome names using file [None]', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--rmsk', help='rmsk to use', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--indels', help='indels vcf to remove', default=None, action=set_abspath())
parser_bcf2ref.add_argument('--o', help='output prefix for ref.ann.vcf.gz files [genotyping/refcalls.flt.vcf.gz]', default='genotyping/refcalls.flt.vcf.gz')

# velvet
parser_velvet = subparsers.add_parser('velvet', help='Run velvet denovo assembly', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=140), usage='genobox.py velvet [options]', description=
   '''Run Velvet denovo assembly, add_velvetg/add_velveth has to be added with quotes, eg: add_velvetg "-very_clean yes".\nFormat can be: fasta, fastq, raw, fasta.gz, fastq.gz, raw.gz, sam, bam.\nDefault ksizes (auto) will run a range (read_length/3 ... read_length/3*2) with step 4 between \n''')
parser_velvet.add_argument('--short', help='format and short reads', nargs='+', action=required_nargs_abspath(0,3))
parser_velvet.add_argument('--shortPaired', help='format and short paired reads', nargs='+', action=required_nargs_abspath(0,3))
parser_velvet.add_argument('--short2', help='format and short2 reads', nargs='+', action=required_nargs_abspath(0,3))
parser_velvet.add_argument('--shortPaired2', help='format and short paired2 reads', nargs='+', action=required_nargs_abspath(0,3))
parser_velvet.add_argument('--long', help='format and long reads', nargs='+', action=required_nargs_abspath(0,3))
parser_velvet.add_argument('--longPaired', help='format and long paired reads', nargs='+', action=required_nargs_abspath(0,3))
parser_velvet.add_argument('--ksizes', help='kmers to run assemblies for (single (m) or m M s (min, max, step)) [auto]', nargs='+', default=['auto'])
parser_velvet.add_argument('--outpath', help='name of run, also output dir [assembly]', default='assembly')
parser_velvet.add_argument('--trim', help='should input files be trimmed (illumina only) [False]', default=False, action='store_true')
parser_velvet.add_argument('--min_contig_lgth', help='mininum length to report contig [100]', default=100, type=int)
parser_velvet.add_argument('--cov_cut', help='coverage cutoff for removal of low coverage (float) [None]', default=None, type=float)
parser_velvet.add_argument('--exp_cov', help='Expected mean coverage (None, float, auto) [auto]', default='auto')
parser_velvet.add_argument('--ins_length', help='insert size (reads included) [None]', default=None, type=int)
parser_velvet.add_argument('--add_velveth', help='additional parameters to velveth', default=None)
parser_velvet.add_argument('--add_velvetg', help='additional parameters to velvetg', default=None)

# ab
parser_ab = subparsers.add_parser('ab', help='Perform alignment and bam-processing', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py ab [options]')
parser_ab.add_argument('--se', help='input single end fqfiles', nargs='+', action=set_abspath(), default=[])
parser_ab.add_argument('--pe1', help='input paired end fqfiles PE1', nargs='+', action=set_abspath(), default=[])
parser_ab.add_argument('--pe2', help='input paired end fqfiles PE2', nargs='+', action=set_abspath(), default=[])
parser_ab.add_argument('--fa', help='bwa index to map against', action=set_abspath())
parser_ab.add_argument('--libfile', help='input parameter file', action=set_abspath())
parser_ab.add_argument('--libs', help='if --libfile is not given, library names for each input fq', nargs='+', default=['lib'])
parser_ab.add_argument('--mapq', help='if --libfile is not given, mapping quality threshold', type=int, nargs='+', default=[30])
parser_ab.add_argument('--pl', help='if --libfile is not given, platform for each input fq (order: se, pe1, pe2)', nargs='+', default=['ILLUMINA'], action=required_choices(['CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'PACBIO']))
parser_ab.add_argument('--a', help='maximum insert size for bwa sampe (-a) [500]', default=500, type=int)
parser_ab.add_argument('--N', help='maximum number of alignments to output in the XA tag [3]', default=3, type=int)
parser_ab.add_argument('--add_aln', help='additional parameters to bwa aln', default=None)
parser_ab.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/', action=set_abspath())
parser_ab.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int, action=required_interval(0,1000))
parser_ab.add_argument('--outbam', help='output file [alignment/final.flt.sort.rmdup.bam]', default='alignment/final.flt.sort.rmdup.bam')

# abgv
parser_abgv = subparsers.add_parser('abgv', help='Perform alignment, bam-processing, genotyping and filtering', parents=[parent_parser], formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=140), usage='genobox.py abgv [options]')
parser_abgv.add_argument('--se', help='input single end fqfiles', nargs='+', action=set_abspath(), default=[])
parser_abgv.add_argument('--pe1', help='input paired end fqfiles PE1', nargs='+', action=set_abspath(), default=[])
parser_abgv.add_argument('--pe2', help='input paired end fqfiles PE2', nargs='+', action=set_abspath(), default=[])
parser_abgv.add_argument('--fa', help='bwa index to map against', action=set_abspath())
parser_abgv.add_argument('--libfile', help='input parameter file', action=set_abspath())
parser_abgv.add_argument('--libs', help='if --libfile is not given, library names for each input fq', nargs='+', default=['lib'])
parser_abgv.add_argument('--mapq', help='mapping quality threshold', type=int, nargs='+', default=[30])
parser_abgv.add_argument('--pl', help='if --libfile is not given, platform for each input fq (order: se, pe1, pe2)', nargs='+', default=['ILLUMINA'], action=required_choices(['CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 'HELICOS', 'IONTORRENT', 'PACBIO']))
parser_abgv.add_argument('--a', help='maximum insert size for bwa sampe (-a) [500]', default=500, type=int)
parser_abgv.add_argument('--N', help='maximum number of alignments to output in the XA tag [3]', default=3, type=int)
parser_abgv.add_argument('--add_aln', help='additional parameters to bwa aln', default=None)
parser_abgv.add_argument('--bwa6', help='use bwa-6.1 (False)', default=False, action='store_true')
parser_abgv.add_argument('--genome', help='file containing genome to analyse, format: chrom\tchrom_len\tchrom_short_name\tploidy\tmin_depth\tmax_depth\n', default=None, action=set_abspath(), required=True)
parser_abgv.add_argument('--tmpdir', help='temporary dir for rmdup [/panvol1/simon/tmp/]', default='/panvol1/simon/tmp/', action=set_abspath())
parser_abgv.add_argument('--qtrim', help='quality threshold to trim 3\'', default=0, type=int, action=required_interval(0,1000))
parser_abgv.add_argument('--prior', help='type of prior to use for bayesian model (full, flat, cond2) [flat]', default='flat')
parser_abgv.add_argument('--pp', help='posterior probability cutoff [0.001]', default=0.001, type=float, action=required_interval(0,1))
parser_abgv.add_argument('--caller', help='samtools or gatk [samtools]', default='samtools')
parser_abgv.add_argument('--Q', help='minimum quality score', type=float, default=20.0, action=required_interval(0,10000))
parser_abgv.add_argument('--rmsk', help='rmsk to use', default=None, action=set_abspath())
parser_abgv.add_argument('--ab', help='allelic balance threshold [0.5] (no filter)', type=float, default=0.50, action=required_interval(0,1))
parser_abgv.add_argument('--prune', help='distance (nt) to prune within [0] (no filter)', type=int, default=0, action=required_interval(0,100000000))
parser_abgv.add_argument('--ex', help='exhange chromosome names using file [None]', default=None, action=set_abspath())
parser_abgv.add_argument('--dbsnp', help='dbsnp file to use (vcf.gz format)', default=None, action=set_abspath())
parser_abgv.add_argument('--ovar', help='output variant vcf-file [genotyping/\"name\".var.flt.vcf.gz]', default='genotyping/var.flt.vcf.gz')
parser_abgv.add_argument('--oref', help='output reference vcf-file [genotyping/\"name\".ref.flt.vcf.gz]', default='genotyping/ref.flt.vcf.gz')
parser_abgv.add_argument('--denovo', help='perform denovo assembly of unmapped reads [False]', default=False, action='store_true')


# parse args
args = parser.parse_args()
#args = parser.parse_args('alignment --se SRR002081se.recal.fastq SRR002082se.recal.fastq --pe1 SRR002137pe_1.recal.fastq SRR002138pe_1.recal.fastq --pe2 SRR002137pe_2.recal.fastq SRR002138pe_2.recal.fastq --n 16 --sample NA12891 --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_rCRS.fa'.split(' '))
#args = parser.parse_args('alignment --se Kleb-10-213361_2_1_sequence.trim.fq --fa kleb_pneu.fa --sample kleb_10_213361se --n 16'.split(' '))
#args = parser.parse_args('bamprocess --bam alignment/SRR002081se.recal.fastq.bam alignment/SRR002082se.recal.fastq.bam alignment/SRR002137pe_1.recal.fastq.bam alignment/SRR002138pe_1.recal.fastq.bam --mapq 30 --libs libA'.split(' '))
#args = parser.parse_args('bamprocess --bam alignment/Kleb-10-213361_2_1_sequence.trim.fq.bam --mapq 30 --outbam alignment/kleb_10_213361.flt.sort.rmdup.bam --n 16'.split(' '))
#args = parser.parse_args('bamstats --bam alignment/kleb_10_213361.flt.sort.rmdup.bam --genome ../kleb_pneu.genome'.split(' '))
#args = parser.parse_args('genotyping --bam alignment/kleb_10_213361.flt.sort.rmdup.bam --genome kleb_pneu.genome --fa ../kleb_pneu.fa --prior full --pp 0.01 --var --o genotyping/kleb_10_213361.bcf'.split(' '))
#args = parser.parse_args('vcffilter --bcf genotyping/kleb_10_213361.bcf --Q 40 --prune 5 --ab 0.20 --genome kleb_pneu.genome --o genotyping/kleb_10_213361.vcf'.split(' '))
#args = parser.parse_args('abgv --pe1 Kleb-10-213361_2_1_sequence.txt --pe2 Kleb-10-213361_2_2_sequence.txt --fa kleb_pneu.fa --genome kleb_pneu.genome --ab 0.2 --prune 5 --sample kleb_10_213361'.split(' '))
#args = parser.parse_args(' abgv --se data/*.trim --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa  --libfile libs.AusAboriginal.txt --genome build37_rCRS.genome --prior flat --pp 0.0001 --Q 40 --rmsk rmsk_build37_rCRS.number.sort.genome --ex gi2number.build37_rCRS --ab 0.2 --prune 5 --n 4 --ovar genotyping/AusAboriginal.var.flt.vcf.gz --oref genotyping/AusAboriginal.ref.flt.vcf.gz --sample AusAboriginal'.split())
#args = parser.parse_args('trim --pe1 Kleb-10-213361_2_1_sequence.txt --pe2 Kleb-10-213361_2_2_sequence.txt --l 15'.split())
#args = parser.parse_args('alignment --pe1 SRR002137pe_1.recal.fastq SRR002138pe_1.recal.fastq --pe2 SRR002137pe_2.recal.fastq SRR002138pe_2.recal.fastq --se SRR002081se.recal.fastq SRR002082se.recal.fastq --sample NA12891 --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa --libfile libs.NA12891.txt'.split())
#args = parser.parse_args('alignment --pe1 SRR002138pe_1.recal.fastq --pe2 SRR002138pe_2.recal.fastq --sample NA12891 --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa --libfile libs.NA12891.txt'.split())
#args = parser.parse_args('alignment --pe1 SRR002138pe_1.recal.fastq --pe2 SRR002138pe_2.recal.fastq --sample NA12891 --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa --libs A A --pl ILLUMINA ILLUMINA'.split())
#args = parser.parse_args('abgv --pe1 SRR002138pe_1.recal.fastq --pe2 SRR002138pe_2.recal.fastq --sample NA12891 --fa /panvol1/simon/databases/hs_ref37_rCRS/hs_ref_GRCh37_all.fa --mapq 30 30 --libs A A --pl ILLUMINA ILLUMINA'.split())
#args = parser.parse_args('alignment --se SRR075153.fastq SRR075154.fastq --sample Vcholerae_C6 --fa /panvol1/simon/databases/bacteria/vibrio_cholerae_O1_N16961.fa --libfile lib.C6.txt.2'.split())
#args = parser.parse_args('abgv --no_trim --se /panvol1/simon/projects/cge/haiti/data/ebi/SRP004712/SRR075153.fastq.trim.fq /panvol1/simon/projects/cge/haiti/data/ebi/SRP004712/SRR075154.fastq.trim.fq --sample Vcholerae_C6 --fa /panvol1/simon/databases/bacteria/vibrio_cholerae_O1_N16961.fa --genome /panvol1/simon/databases/bacteria/vibrio_cholerae_O1_N16961.genome --libfile libs.C6.txt.3'.split())
#args = parser.parse_args('ab --se /panvol1/simon/projects/cge/test/1000G/tmp_se.fastq.g --pe1 /panvol1/simon/projects/cge/test/1000G/tmp_pe_1.fastq.gz --pe2 /panvol1/simon/projects/cge/test/1000G/tmp_pe_2.fastq.gz --fa /panvol1/simon/databases/hs_ref37_gatk/human_g1k_v37.fasta --libfile libs.tmp.txt --sample tmp_testing --gz'.split())
#args = parser.parse_args('ab --se /panfs1/ponten_tom/trimmed/00020.s_1_1.fq_filtered.fq.gz.trim.fq.single.gz /panfs1/ponten_tom/trimmed/00020.s_1_2.fq_filtered.fq.gz.trim.fq.single.gz --pe1 /panfs1/ponten_tom/trimmed/00020.s_1_1.fq_filtered.fq.gz.trim.fq.gz --pe2 /panfs1/ponten_tom/trimmed/00020.s_1_2.fq_filtered.fq.gz.trim.fq.gz --libfile libs.tommerup.txt --fa /panvol1/simon/databases/hs_ref37_gatk/human_g1k_v37.fasta --a 3000 --n 4 --gz --sample test --partition cge-cluster'.split())

#args = parser.parse_args('trim --pe1 giantSquid01A_170bp.fq  --pe2 giantSquid01B_170bp.fq --min_length 60 --min_baseq 20 --min_avgq 20 --adaptors GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG ACACTCTTTCCCTACACGACGCTCTTCCGATC CTCGGCATTCCTGCTGAACCGCTCTTCCGATC'.split())


# If working dir is given, create and move to working directory else run where program is invoked

if args.sample:
   if not os.path.exists(args.sample):
      os.makedirs(args.sample)
   os.chmod(args.sample, 0777)
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
if args.module == 'trim':
   from genobox_trim import *
   trimfiles = start_trim(args, logger)

if args.module == 'alignment':
   from genobox_alignment import *
   bamfiles = start_alignment(args, logger)

elif args.module == 'bamprocess':
   from genobox_bamprocess import *
   final_bam = start_bamprocess(args.libfile, args.bam, args.mapq, args.libs, args.tmpdir, args.queue, args.outbam, args.sample, args.partition, logger)

elif args.module == 'bamstats':
   from genobox_bamstats import *
   start_bamstats(args, args.bam, args.partition, logger)

elif args.module == 'genotyping':
   from genobox_genotyping import *
   final_bcf = start_genotyping(args.bam, args.genome, args.fa, args.prior, args.pp, args.queue, args.o, args.sample, args.partition, logger)

elif args.module == 'vcffilter':
   from genobox_vcffilter import *
   final_vcf = start_vcffilter(args.bcf, args.genome, args.caller, args.Q, args.ex, args.rmsk, args.ab, args.prune, args.o, args.queue, args.sample, args.partition, logger)

elif args.module == 'dbsnp':
   from genobox_dbsnp import *
   start_dbsnp(args.vcf, args.ex, args.dbsnp, args.o, args.queue, args.partition, logger)

elif args.module == 'bcf2ref':
   from genobox_bcf2ref import *
   start_bcf2ref(args.bcf, args.genome, args.Q, args.ex, args.dbsnp, args.rmsk, args.indels, args.o, args.queue, args.sample, args.partition, logger)

elif args.module == 'ab':
   from genobox_ab import *
   start_ab(args, logger)
   
elif args.module == 'abgv':
   from genobox_abgv import *
   start_abgv(args, logger)

elif args.module == 'velvet':
   from genobox_denovo_velvet import *
   start_assembly(args, logger)
else:
   pass


