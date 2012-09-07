#!/panvol1/simon/bin/python2.7

def start_vcffilter_gatk(vcfs, genome, fa, Q, rmsk, ab, prune, queue, dir, partition, logger):
   '''Start variant vcf-filter using gatk'''
   
   import genobox_modules
   from genobox_classes import Moab
   from genobox_classes import Semaphore   
   import subprocess
   import os
   
   if not os.path.exists('genotyping'):
      os.makedirs('genotyping')
   
   if not os.path.exists('tmp'):
      os.makedirs('tmp')

   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
   cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=7gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   
   vcffilter_calls = []
   cmd = paths['genobox_home'] + 'genobox_vcffilter_gatk_h.py'
   
   # for each chromosome
   for v in vcfs:
      arg = ' --vcf %s --fa %s --genome %s --Q %f' % (v, fa, genome, Q)
      if rmsk: arg = arg + ' --rmsk %s' % rmsk
      if ab != 0.5: arg = arg + ' --ab %f' % ab
      if prune != 0: arg = arg + ' --prune %i' % prune
      vcffilter_calls.append(cmd+arg)
   
   # submit jobs
   print "Submitting jobs"
   vcffilter_moab = Moab(vcffilter_calls, logfile=logger, runname='run_genobox_vcffilter_gatk', queue=queue, cpu=cpuF, partition=partition)
   
   # release jobs #
   print "Releasing jobs"
   vcffilter_moab.release()
   
   # semaphore
   print "Waiting for jobs to finish ..."
   s = Semaphore(vcffilter_moab.ids, home, 'vcffilter_gatk', queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   
   # return filename of final vcf
   