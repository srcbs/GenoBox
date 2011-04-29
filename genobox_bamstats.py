#!/panvol1/simon/bin/python2.7

def sam_flagstat(bam):
   '''Start samtools flagstat'''
   
   import os
   import genobox_moab
   paths = genobox_moab.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = paths['samtools_home'] + 'samtools flagstat %s > stats/%s.flagstat' % (bam, bamf)
   return [call]

def bed_genomeCov(bam, genome):
   '''Start bedtools genomeCoverageBed'''
   
   import os
   import genobox_moab
   paths = genobox_moab.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = paths['bedtools_home'] + 'genomeCoverageBed -ibam %s -g %s > stats/%s.coverage' % (bam, genome, bamf)
   return [call]

def plot_coverage(bam):
   '''Use output from genomeCoverageBed to plot coverage plots'''
   
   import os
   import genobox_moab
   paths = genobox_moab.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = 'R-2.12 --vanilla stats/%s.coverage %s stats/%s.coverage.pdf < %sgenobox_plotcov.R' % (bamf, bamf, bamf, paths['genobox_home'])
   return [call]

def python_avgdepth(bam):
   '''Start genobox_bam2avgdepth.py'''
   
   import os
   import genobox_moab
   paths = genobox_moab.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = 'python2.7 ' + paths['genobox_home'] + 'genobox_bam2avgdepth.py %s > stats/%s.avgdepth' % (bam, bamf)
   return [call]
   
def start_bamstats(args, bam, logger, wait=True):
   '''Starts calculation of bam statistics'''
   
   # samtools flagstat
   # bedtools genomeCoverageBed
   # python avgdepth
   
   import subprocess
   import genobox_moab
   import os
   
   if not os.path.exists('stats'):
      os.makedirs('stats')
   
   # set queueing
   paths = genobox_moab.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb'
   cpuC = 'nodes=1:ppn=1,mem=2gb'
   cpuE = 'nodes=1:ppn=1,mem=5gb'
   cpuF = 'nodes=1:ppn=2,mem=2gb'
   cpuB = 'nodes=1:ppn=16,mem=10gb'
      
   # create calls
   flagstat_calls = sam_flagstat(bam)
   coverage_calls = bed_genomeCov(bam, args.genome)
   plotcoverage_calls = plot_coverage(bam)
   avgdepth_calls = python_avgdepth(bam)
   
   # submit jobs
   print "Submitting jobs"
   flagstat_ids = genobox_moab.submitjob(flagstat_calls, home, paths, logger, 'run_genobox_flagstat', args.queue, cpuC, False)
   coverage_ids = genobox_moab.submitjob(coverage_calls, home, paths, logger, 'run_genobox_coverage', args.queue, cpuC, False)
   plotcoverage_ids = genobox_moab.submitjob(plotcoverage_calls, home, paths, logger, 'run_genobox_plotcoverage', args.queue, cpuA, True, 'one2one', 1, True, *coverage_ids)
   avgdepth_ids = genobox_moab.submitjob(avgdepth_calls, home, paths, logger, 'run_genobox_avgdepth', args.queue, cpuE, False)
   
   # release jobs
   allids = []
   allids.extend(flagstat_ids) ; allids.extend(coverage_ids) ; allids.extend(plotcoverage_ids) ; allids.extend(avgdepth_ids)
   releasemsg = genobox_moab.releasejobs(allids)
   semaphore_ids = allids
   
   # wait for jobs to finish
   if wait:
      print "Waiting for jobs to finish ..." 
      genobox_moab.wait_semaphore(semaphore_ids, home, 'bam_stats', args.queue, 20, 86400)
      print "--------------------------------------"
   else:
      print "Jobs running, continuing"
      print "--------------------------------------"
   
   