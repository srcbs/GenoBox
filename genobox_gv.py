#!/panvol1/simon/bin/python2.7

from genobox_genotyping import *
from genobox_vcffilter import *
from genobox_genotyping_gatk import *
from genobox_vcffilter_gatk import *
import genobox_modules


def start_gv(args, logger):
   '''Perform alignment and bam processing'''
   
   import os
   import subprocess
   
   genobox_modules.check_genome(args.genome)
   final_bcf = 'genotyping/%s.all.bcf' % args.sample
   
   # start run
   if args.sample:
      print "--------------------------------------"
      print "Processing sample: %s" % args.sample
   print "--------------------------------------"
     
   if args.caller == 'samtools':
      print "Starting genotyping (samtools)"
      final_bcf = start_genotyping(args.bam, args.genome, args.fa, args.prior, args.pp, args.queue, final_bcf, args.sample, args.partition, logger)
      
      print "Starting vcffiltering"
      final_vcf = start_vcffilter(final_bcf, args.genome, args.caller, args.Q, args.ex, args.rmsk, args.ab, args.prune, args.ovar, args.queue, args.sample, args.partition, logger)
      
      print "Start dbsnp"
      final_dbsnp_vcf = start_dbsnp(final_vcf, args.ex, args.dbsnp, args.ovar, args.queue, args.partition, logger)
      
      print "Start bcf2ref"
      start_bcf2ref(final_bcf, args.genome, args.Q, args.ex, args.dbsnp, args.rmsk, 'genotyping/indels_for_filtering.vcf', args.oref, args.queue, args.sample, args.partition, logger)
   
   elif args.caller ==  'gatk':
      print "Start genotyping (gatk)"
      vcffiles = start_genotyping_gatk(args.bam, args.genome, args.fa, args.dbsnp, args.call_conf, args.call_emit, args.output_mode, args.queue, args.sample, args.partition, logger)
      
      print "Start vcffiltering (gatk)"
      final_vcfs = start_vcffilter_gatk(vcffiles, args.genome, args.fa, args.Q, args.rmsk, args.ab, args.prune, args.queue, args.sample, args.partition, logger)
   
   # remove queuing system outfiles
   genobox_modules.rm_files(['run_genobox_*', 'semaphores.*'])
   
   print "Done"
   print "--------------------------------------"
