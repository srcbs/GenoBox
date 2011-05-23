#!/usr/bin/python

def start_dbsnp(vcf, ex, dbsnp, o, queue, logger):
   '''Annotate vcf.gz file with dbSNP,
   exchanging chromsome names to dbSNP version
   sort vcf and the input to dbSNP
   '''
   
   import genobox_modules
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
   cpuA = 'nodes=1:ppn=1,mem=512mb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   cpuE = 'nodes=1:ppn=1,mem=5gb'
   cpuF = 'nodes=1:ppn=2,mem=2gb'
   cpuB = 'nodes=1:ppn=16,mem=10gb'
   
   # create command
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_dbsnp_h.py'
   arg = ' --vcf %s --ex %s --dbsnp %s --o %s' % (vcf, ex, dbsnp, o)
   dbsnp_calls = [cmd+arg]
   
   # submit jobs
   print "Submitting jobs"
   dbsnp_ids = genobox_modules.submitjob(dbsnp_calls, home, paths, logger, 'run_genobox_dbsnp', queue, cpuC, False)
   
   # release jobs #
   allids = []
   allids.extend(dbsnp_ids)
   releasemsg = genobox_modules.releasejobs(allids)
   
   # semaphore
   print "Waiting for jobs to finish ..."
   genobox_modules.wait_semaphore(dbsnp_ids, home, 'dbsnp', queue, 20, 2*86400)
   print "--------------------------------------"
   
   return o
