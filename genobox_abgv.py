#!/usr/bin/python

from genobox_alignment import *
from genobox_bamprocess import *
from genobox_genotyping import *
from genobox_vcffilter import *

def start_abgv(args, logger):
   '''Start alignment, bam processing, genotyping and vcf filtering'''
      
   print "Starting alignment"
   bamfiles = start_alignment(args.se, args.pe1, args.pe2, args.bwaindex, 'alignment/', args.qtrim, args.a, args.r, args.n, args.queue, logger)
   print "Starting bamprocessing"
   final_bam = start_bamprocess(args.lib_file, bamfiles, args.mapq, args.libs, args.tmpdir, args.queue, 'alignment/final.flt.sort.rmdup.bam', logger)
   print "Starting genotyping"
   final_bcf = start_genotyping(final_bam, args.genome, args.bwaindex, args.prior, args.pp, args.queue, 'genotyping/snpcalls.bcf', logger)
   print "Starting vcffiltering"
   final_vcf = start_vcffilter(final_bcf, args.genome, args.caller, args.Q, args.rmsk, args.ab, args.prune, args.o, args.queue, logger)
   print "Done, genotyping is written in %s" % args.o
   