#!/panvol1/simon/bin/python2.7

from genobox_alignment import *
from genobox_bamprocess import *
from genobox_bamstats import *

import genobox_modules         

# test
#import pickle
#bamfiles = pickle.load(open('bamfiles', 'r'))
#libfile = pickle.load(open('libfile', 'r'))

def start_ab(args, logger):
   '''Perform alignment and bam processing'''
   
   import os
   import subprocess
   
   final_bam = args.outbam
   
   # initialize library file from given arguments
   library = genobox_modules.initialize_library(args.libfile, args.se, args.pe1, args.pe2, args.sample, args.mapq, args.libs, args.pl)
   
   # start run
   if args.sample:
      print "--------------------------------------"
      print "Processing sample: %s" % args.sample
   print "--------------------------------------"
      
   print "Starting alignment"
   (bamfiles, library) = start_alignment(args, logger)
   print "Starting bam processing"
   final_bam = start_bamprocess(library, genobox_modules.unique(bamfiles.values()), args.mapq, args.libs, args.tmpdir, args.queue, final_bam, args.realignment, args.known, args.fa, args.sample, args.partition, logger)
   
   # remove queuing system outfiles
   genobox_modules.rm_files(['run_genobox_*', 'semaphores.*'])
   
   print "Done"
   print "--------------------------------------"
