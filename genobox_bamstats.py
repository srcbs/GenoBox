#!/panvol1/simon/bin/python2.7

def sam_flagstat(bam):
   '''Start samtools flagstat'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = paths['samtools_home'] + 'samtools flagstat %s > stats/%s.flagstat' % (bam, bamf)
   return [call]

def bed_genomeCov(bam, genome):
   '''Start bedtools genomeCoverageBed'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = paths['bedtools_home'] + 'genomeCoverageBed -ibam %s -g %s > stats/%s.coverage' % (bam, genome, bamf)
   return [call]

def plot_coverage(bam):
   '''Use output from genomeCoverageBed to plot coverage plots'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = 'R-2.12 --vanilla stats/%s.coverage %s stats/%s.coverage.pdf < %sgenobox_plotcov.R' % (bamf, bamf, bamf, paths['genobox_home'])
   return [call]

def python_avgdepth(bam):
   '''Start genobox_bam2avgdepth.py'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
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
   import genobox_modules
   from genobox_classes import Moab
   from genobox_classes import Semaphore   
   import os
   
   if not os.path.exists('stats'):
      os.makedirs('stats')
   
   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
   cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
      
   # create calls
   flagstat_calls = sam_flagstat(bam)
   coverage_calls = bed_genomeCov(bam, args.genome)
   plotcoverage_calls = plot_coverage(bam)
   avgdepth_calls = python_avgdepth(bam)
   
   # submit jobs
   print "Submitting jobs"
   flagstat_moab = Moab(flagstat_calls, logfile=logger, runname='run_genobox_flagstat', queue=args.queue, cpu=cpuC)
   coverage_moab = Moab(coverage_calls, logfile=logger, runname='run_genobox_coverage', queue=args.queue, cpu=cpuC)
   plotcoverage_moab = Moab(plotcoverage_calls, logfile=logger, runname='run_genobox_plotcoverage', queue=args.queue, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=coverage_moab.ids)
   avgdepth_moab = Moab(avgdepth_calls, logfile=logger, runname='run_genobox_avgdepth', queue=args.queue, cpu=cpuE)
   
   # release jobs
   print "Releasing jobs"
   flagstat_moab.release()
   coverage_moab.release()
   plotcoverage_moab.release()
   avgdepth_moab.release()
      
   # wait for jobs to finish
   if wait:
      print "Waiting for jobs to finish ..."
      semaphore_ids = flagstat_moab.ids + coverage_moab.ids + plotcoverage_moab.ids + avgdepth_moab.ids
      s = Semaphore(semaphore_ids, home, 'bam_stats', args.queue, 20, 86400) 
      s.wait()
      print "--------------------------------------"
   else:
      print "Jobs running, continuing"
      print "--------------------------------------"
   
   