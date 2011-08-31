#!/panvol1/simon/bin/python2.7

import argparse

def extract_unmapped_reads(bamfiles):
   '''Generate calls to extract unmapped reads from bamfiles'''
   
   import genobox_modules
      
   paths = genobox_modules.setSystem()
   cmd = '%ssamtools view -h -b -f 4' % (paths['samtools_home'])
   
   
   calls = []
   unmapped = {}
   for id,bam in bamfiles.items():
      unmap_bam = bam+'.unmapped.bam'
      unmapped[id] = unmap_bam
      arg = ' %s > %s' % (bam, unmap_bam)
      calls.append(cmd+arg)
   
   return (calls, unmapped)

def start_unmapped_assembly():
   '''Start denovo assembly of unmapped reads given by args'''
   
   # args is different when from genobox.py
   from genobox_classes import Library
   import genobox_modules
   
   #library = Library(args.libfile)
   library = Library("libs.NA12891.txt.RISO17V9VS")
   library.read()
   #bamfiles = library.getValues('ID', 'BAM')
   libs = library.getValues('LB', 'BAM')
   
   (unmapped_calls, unmapped) = extract_unmapped_reads(bamfiles)
   # update library
   library.update_with_tag('ID', 'UNM', unmapped, force=True)
   
   # combine unmapped files to library level
   
   
   # need to get all bamfiles
   # create calls to extract unmapped reads from bamfiles: samtools view -h -b -f 4 Vchol-016_7_1_sequence.trim.fq.bam > unmapped.bam
   # create calls to run denovo assembly: genobox.py velvet --shortPaired bam unmapped.bam --sample unmapped_assembly --ksizes 33 85 4 --add_velvetg "-very_clean yes"
   


