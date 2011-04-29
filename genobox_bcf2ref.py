#!/usr/bin/python


def get_genome(chr_file):
   '''Read genome file into list of list'''
   fh = open(chr_file, 'r')
   L = []
   for line in fh:
      line = line.rstrip()
      s = line.split('\t')
      L.append(s)
   return L

def start_bcf2ref(bcf, genome_file, Q, ex, dbsnp, rmsk, indels, o, queue, dir, logger):
   '''Extract high confidence same-as-reference bases from bcf, options are to:
   
   exchange ids
   annotate using dbsnp
   filter rmsk
   filter ambiguous indel positions
   '''
   
   import genobox_moab
   import subprocess
   import os
   
   if not os.path.exists('genotyping'):
      os.makedirs('genotyping')
   
   # set queueing
   paths = genobox_moab.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   cpuE = 'nodes=1:ppn=1,mem=5gb'
   cpuF = 'nodes=1:ppn=2,mem=2gb'
   cpuB = 'nodes=1:ppn=16,mem=10gb'
   
   # read genome file
   genome = get_genome(genome_file)
   
   # create commands
   bcf2ref_calls = []
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_bcf2ref_h.py'
   for chr in genome:
      # set outfile name
      if len(genome) == 1:
         if dir and dir != 'None':
            outfile = '%s/%s.%s' % (os.path.split(o)[0], dir, os.path.split(o)[1])
         else:
            outfile = o
      else:
         if dir and dir != 'None':
            outfile = '%s/%s.%s.%s' % (os.path.split(o)[0], dir, chr[2], os.path.split(o)[1])
         else:
            outfile = '%s/%s.%s' % (os.path.split(o)[0], chr[2], os.path.split(o)[1])
      
      arg = ' --bcf %s --chr_id \"%s\" --chr %s --d %s --D %s --Q %f --ex %s --dbsnp %s --rmsk %s --indels %s --o %s' % (bcf, chr[0], chr[2], chr[4], chr[5], Q, ex, dbsnp, rmsk, indels, outfile)
      bcf2ref_calls.append(cmd+arg)
   
   # submit jobs
   print "Submitting jobs"
   bcf2ref_ids = genobox_moab.submitjob(bcf2ref_calls, home, paths, logger, 'run_genobox_bcf2ref', queue, cpuE, False)
   
   # release jobs #
   allids = []
   allids.extend(bcf2ref_ids)
   releasemsg = genobox_moab.releasejobs(allids)
   
   # semaphore
   print "Waiting for jobs to finish ..."
   genobox_moab.wait_semaphore(bcf2ref_ids, home, 'bcf2ref', queue, 20, 2*86400)
   print "--------------------------------------"
