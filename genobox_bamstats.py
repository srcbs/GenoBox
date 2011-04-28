#!/panvol1/simon/bin/python2.7

def sam_flagstat(bam):
   '''Start samtools flagstat'''
   
   import os
   import pipelinemod
   paths = pipelinemod.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = paths['samtools_home'] + 'samtools flagstat %s > stats/%s.flagstat' % (bam, bamf)
   return [call]

def bed_genomeCov(bam, genome):
   '''Start bedtools genomeCoverageBed'''
   
   import os
   import pipelinemod
   paths = pipelinemod.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = paths['bedtools_home'] + 'genomeCoverageBed -ibam %s -g %s > stats/%s.coverage' % (bam, genome, bamf)
   return [call]

def plot_coverage(bam):
   '''Use output from genomeCoverageBed to plot coverage plots'''
   
   import os
   import pipelinemod
   paths = pipelinemod.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = 'R-2.12 --vanilla stats/%s.coverage %s stats/%s.coverage.pdf < %sgenobox_plotcov.R' % (bamf, bamf, bamf, paths['genobox_home'])
   return [call]

def python_avgdepth(bam):
   '''Start genobox_bam2avgdepth.py'''
   
   import os
   import pipelinemod
   paths = pipelinemod.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = 'python2.7' + paths['genobox_home'] + 'genobox_bam2avgdepth.py --i %s > stats/%s.avgdepth' % (bam, bamf)
   return [call]
   
def start_bamstats(args, logger, wait=True):
   '''Starts calculation of bam statistics'''
   
   # samtools flagstat
   # bedtools genomeCoverageBed
   # python avgdepth
   
   import subprocess
   import pipelinemod
   import os
   
   if not os.path.exists('stats'):
      os.makedirs('stats')
   
   # set queueing
   paths = pipelinemod.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   cpuE = 'nodes=1:ppn=1,mem=5gb'
   cpuF = 'nodes=1:ppn=2,mem=2gb'
   cpuB = 'nodes=1:ppn=16,mem=10gb'
      
   # create calls
   flagstat_calls = sam_flagstat(args.bam)
   coverage_calls = bed_genomeCov(args.bam, args.genome)
   plotcoverage_calls = plot_coverage(args.bam)
   avgdepth_calls = python_avgdepth(args.bam)
   
   # submit jobs
   print "Submitting jobs"
   flagstat_ids = pipelinemod.submitjob(flagstat_calls, home, paths, logger, 'run_genobox_flagstat', args.queue, cpuC, False)
   coverage_ids = pipelinemod.submitjob(coverage_calls, home, paths, logger, 'run_genobox_coverage', args.queue, cpuC, False)
   plotcoverage_ids = pipelinemod.submitjob(plotcoverage_calls, home, paths, logger, 'run_genobox_plotcoverage', args.queue, cpuA, True, 'one2one', 1, True, *coverage_ids)
   avgdepth_ids = pipelinemod.submitjob(avgdepth_calls, home, paths, logger, 'run_genobox_avgdepth', args.queue, cpuE, False)
   
   # release jobs
   allids = []
   allids.extend(flagstat_ids) ; allids.extend(coverage_ids) ; allids.extend(plotcoverage_ids) ; allids.extend(avgdepth_ids)
   releasemsg = pipelinemod.releasejobs(allids)
   semaphore_ids = allids
   
   # wait for jobs to finish
   if wait:
      print "Waiting for jobs to finish ..." 
      pipelinemod.wait_semaphore(semaphore_ids, home, 'bam_stats', args.queue, 20, 86400)
      print "--------------------------------------"
   else:
      print "Jobs running, continuing"
      print "--------------------------------------"
   
   