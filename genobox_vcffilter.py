#!/usr/bin/python

def start_vcffilter(bcf, genome, caller, Q, ex, rmsk, ab, prune, o, queue, dir, partition, logger):
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
   from genobox_classes import Moab
   from genobox_classes import Semaphore   
   import subprocess
   import os
   
   if not os.path.exists('genotyping'):
      os.makedirs('genotyping')
   
   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
   cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   
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
   vcffilter_moab = Moab(vcffilter_calls, logfile=logger, runname='run_genobox_vcffilter', queue=queue, cpu=cpuE, partition=partition)
   
   # release jobs #
   print "Releasing jobs"
   vcffilter_moab.release()
   
   # semaphore
   print "Waiting for jobs to finish ..."
   s = Semaphore(vcffilter_moab.ids, home, 'vcffilter', queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   
   # return filename of final vcf
   return o