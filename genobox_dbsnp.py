#!/usr/bin/python

def start_dbsnp(vcf, ex, dbsnp, o, queue, logger):
   '''Annotate vcf.gz file with dbSNP,
   exchanging chromsome names to dbSNP version
   sort vcf and the input to dbSNP
   '''
   
   import genobox_modules
   from genobox_classes import Moab
   from genobox_classes import Semaphore   
   import subprocess
   import os
   
   if not dbsnp or dbsnp == 'None':
      print "No dbsnp file given - skipping"
      print "--------------------------------------"
      return vcf
   
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
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_dbsnp_h.py'
   arg = ' --vcf %s --ex %s --dbsnp %s --o %s' % (vcf, ex, dbsnp, o)
   dbsnp_calls = [cmd+arg]
   
   # submit jobs
   print "Submitting jobs"
   dbsnp_moab = Moab(dbsnp_calls, logfile=logger, runname='run_genobox_dbsnp', queue=queue, cpu=cpuC)
   
   # release jobs #
   print "Releasing jobs"
   dbsnp_moab.release()
   
   # semaphore
   print "Waiting for jobs to finish ..."
   s = Semaphore(dbsnp_moab.ids, home, 'dbsnp', queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   
   return o
