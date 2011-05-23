#!/usr/bin/python

def start_vcffilter(bcf, genome, caller, Q, ex, rmsk, ab, prune, o, queue, dir, logger):
   '''Start variant vcf-filter
   
   Genome file must be given, format is a line for each chromosome:
   chrom\tchrom_len\tchrom_short_name\haploid/diploid\tlow_depth\thigh_depth
   
   Filtering steps:
   vcfutils.pl varFilter
   annotated repeats using rmsk
   heterozygote variants on haploid chromosomes
   allelic balance
   pruning of variants within N nt of each other
   '''
   
   import genobox_modules
   import subprocess
   import os
   
   if not os.path.exists('genotyping'):
      os.makedirs('genotyping')
   
   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   cpuE = 'nodes=1:ppn=1,mem=5gb'
   cpuF = 'nodes=1:ppn=2,mem=2gb'
   cpuB = 'nodes=1:ppn=16,mem=10gb'
   
   # create command
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_vcffilter_h.py'
   if dir and dir != 'None':
      outfile = '%s/%s.%s' % (os.path.split(o)[0], dir, os.path.split(o)[1])
   else:
      outfile = o
   arg = ' --bcf %s --genome %s --caller %s --Q %f --rmsk %s --ab %f --prune %i --o %s' % (bcf, genome, caller, Q, rmsk, ab, prune, outfile)
   vcffilter_calls = [cmd+arg]
   
   # submit jobs
   print "Submitting jobs"
   vcffilter_ids = genobox_modules.submitjob(vcffilter_calls, home, paths, logger, 'run_genobox_vcffilter', queue, cpuE, False)
   
   # release jobs #
   allids = []
   allids.extend(vcffilter_ids)
   releasemsg = genobox_modules.releasejobs(allids)
   
   # semaphore
   print "Waiting for jobs to finish ..."
   genobox_modules.wait_semaphore(vcffilter_ids, home, 'vcffilter', queue, 20, 2*86400)
   print "--------------------------------------"
   
   # return filename of final vcf
   return o