#!/usr/bin/python

from genobox_alignment import *
from genobox_bamprocess import *
from genobox_bamstats import *
from genobox_genotyping import *
from genobox_vcffilter import *
from genobox_dbsnp import *
from genobox_bcf2ref import *
import genobox_moab
import genobox_modules


def start_abgv(args, logger):
   '''Start alignment, bam processing, genotyping, vcffiltering, dbsnp annotation, bcf2ref'''
   
   final_bam = 'alignment/%s.flt.sort.rmdup.bam' % args.sample
   final_bcf = 'genotyping/%s.all.bcf' % args.sample
   
   print "--------------------------------------"
   print "Starting alignment"
   (bamfiles, libfile) = start_alignment(args, logger)
   print "Starting bam processing"
   final_bam = start_bamprocess(libfile, genobox_modules.unique(bamfiles.values()), args.mapq, args.libs, args.tmpdir, args.queue, final_bam, args.sample, logger)
   print "Starting bam stats"
   start_bamstats(args, final_bam, logger, wait=False)
   print "Starting genotyping"
   final_bcf = start_genotyping(final_bam, args.genome, args.fa, args.prior, args.pp, args.queue, final_bcf, logger)
   print "Starting vcffiltering"
   final_vcf = start_vcffilter(final_bcf, args.genome, args.caller, args.Q, args.rmsk, args.ab, args.prune, args.ovar, args.queue, args.sample, logger)
   print "Start dbsnp"
   final_dbsnp_vcf = start_dbsnp(final_vcf, args.ex, args.dbsnp, args.ovar, args.queue, logger)
   print "Start bcf2ref"
   start_bcf2ref(final_bcf, args.genome, args.Q, args.ex, args.dbsnp, args.rmsk, 'genotyping/indels_for_filtering.vcf', args.oref, args.queue, args.sample, logger)
   
   # remove queuing system outfiles
   genobox_moab.rm_files(['run_genobox_*', 'semaphores.*'])
   
   print "Done"
   print "--------------------------------------"
   print "Raw genotyping is written in genotyping/all.bcf"
   print "High confidence variants: %s" % args.ovar
   print "High confidence reference: %s" % args.oref
   print "--------------------------------------"
   
   