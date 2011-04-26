#!/usr/bin/python

def start_vcffilter(bcf, genome, caller, Q, rmsk, ab, prune, o, queue, logger):
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
   
   import pipelinemod
   import subprocess
   import os
   
   if not os.path.exists('genotyping'):
      os.makedirs('genotyping')
   
   # set queueing
   paths = pipelinemod.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   cpuE = 'nodes=1:ppn=1,mem=5gb'
   cpuF = 'nodes=1:ppn=2,mem=2gb'
   cpuB = 'nodes=1:ppn=16,mem=10gb'
   
   # create command
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_vcffilter_h.py'
   arg = ' --bcf %s --genome %s --caller %s --Q %f --rmsk %s --ab %f --prune %i --o %s' % (bcf, genome, caller, Q, rmsk, ab, prune, o)
   vcffilter_calls = [cmd+arg]
   
   # submit jobs
   print "Submitting jobs"
   vcffilter_ids = pipelinemod.submitjob(vcffilter_calls, home, paths, logger, 'run_genobox_vcffilter', queue, cpuE, False)
   
   # release jobs #
   allids = []
   allids.extend(vcffilter_ids)
   releasemsg = pipelinemod.releasejobs(allids)
   
   # semaphore
   print "Waiting for jobs to finish ..."
   pipelinemod.wait_semaphore(vcffilter_ids, home, 'vcffilter', queue, 20, 2*86400)
