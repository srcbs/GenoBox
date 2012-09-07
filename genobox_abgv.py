#!/panvol1/simon/bin/python2.7

from genobox_trim import *
from genobox_alignment import *
from genobox_bamprocess import *
from genobox_bamstats import *
from genobox_genotyping import *
from genobox_vcffilter import *
from genobox_genotyping_gatk import *
from genobox_vcffilter_gatk import *
from genobox_dbsnp import *
from genobox_bcf2ref import *
from genobox_classes import Library

import genobox_modules         

# test
#import pickle
#bamfiles = pickle.load(open('bamfiles', 'r'))
#libfile = pickle.load(open('libfile', 'r'))

def start_abgv(args, logger):
   '''Start alignment, bam processing, genotyping, vcffiltering, dbsnp annotation, bcf2ref'''
   
   import os
   import subprocess
   
   # check genome file
   genobox_modules.check_genome(args.genome)
   
   final_bam = 'alignment/%s.flt.sort.rmdup.bam' % args.sample
   final_bcf = 'genotyping/%s.all.bcf' % args.sample
   
   # initialize library file from given arguments
   library = genobox_modules.initialize_library(args.libfile, args.se, args.pe1, args.pe2, args.sample, args.mapq, args.libs, args.pl)
   
   # start run
   if args.sample:
      print "--------------------------------------"
      print "Processing sample: %s" % args.sample
   print "--------------------------------------"
   
   # toggle start trimming
   #if args.no_trim == False:
   #   print "Starting trimming"
   #   (se_files, pe1_files, pe2_files) = start_trim(args, logger)
   #   library.update(Trim=se_files+pe1_files+pe2_files)
   
   print "Starting alignment"
   (bamfiles, library) = start_alignment(args, logger)
   
   print "Starting bam processing"
   final_bam = start_bamprocess(library, genobox_modules.unique(bamfiles.values()), args.mapq, args.libs, args.tmpdir, args.queue, final_bam, args.realignment, args.known, args.fa, args.sample, args.partition, logger)
   
   print "Starting bam stats"
   start_bamstats(args, final_bam, args.partition, logger, wait=False)
   
   print "Starting genotyping"
   if args.caller == 'samtools':
      final_bcf = start_genotyping(final_bam, args.genome, args.fa, args.prior, args.pp, args.queue, final_bcf, args.sample, args.partition, logger)
      print "Starting vcffiltering"
      final_vcf = start_vcffilter(final_bcf, args.genome, args.caller, args.Q, args.ex, args.rmsk, args.ab, args.prune, args.ovar, args.queue, args.sample, args.partition, logger)
      print "Start dbsnp"
      final_dbsnp_vcf = start_dbsnp(final_vcf, args.ex, args.dbsnp, args.ovar, args.queue, args.partition, logger)
      print "Start bcf2ref"
      start_bcf2ref(final_bcf, args.genome, args.Q, args.ex, args.dbsnp, args.rmsk, 'genotyping/indels_for_filtering.vcf', args.oref, args.queue, args.sample, args.partition, logger)
   elif args.caller ==  'gatk':
      print "Start genotyping (gatk)"
      vcffiles = start_genotyping_gatk(final_bam, args.genome, args.fa, args.dbsnp, args.call_conf, args.args.call_emit, args.output_mode, args.queue, args.sample, args.partition, logger)
      print "Start vcffiltering (gatk)"
      final_vcfs = start_vcffilter_gatk(vcffiles, args.genome, args.fa, args.Q, args.rmsk, args.ab, args.prune, args.queue, args.sample, args.partition, args.logger)
   
   # remove queuing system outfiles
   genobox_modules.rm_files(['run_genobox_*', 'semaphores.*'])
   
   print "Done"
   print "--------------------------------------"
   print "Raw genotyping is written in genotyping/all.bcf"
   print "High confidence variants: %s" % args.ovar
   print "High confidence reference: %s" % args.oref
   print "--------------------------------------"
   
   