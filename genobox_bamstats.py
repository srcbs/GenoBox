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

def bed_genomeCov(bam):
   '''Start bedtools genomeCoverageBed'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = os.path.split(bam)[1]   
   
   call = paths['bedtools_bin'] + 'bedtools genomecov -ibam %s > stats/%s.coverage' % (bam, bamf)
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
   
   call = paths['genobox_home'] + 'genobox_bam2avgdepth1.py %s > stats/%s.avgdepth' % (bam, bamf)
   return [call]

def mapdamage(bam, fa):
   '''Run mapdamage on input bam'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   # set bam-file sans paths (input is abspath(bam))
   bamf = 'mapdamage_' + os.path.split(bam)[1]
   c1 = paths['mapdamage_home'] + 'mapdamage-0.3.6.pl map -i %s -r %s -c' % (bam, fa)
   
   # create call to move results file to stats dir
   c2 = 'mv %s stats/%s' % (bam, bamf)
   return [c1, c2]
   

def get_saturation(bam):
   '''Perform saturation calculations'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   bamf = os.path.split(bam)[1] + '.saturation'
   c1 = paths['genobox_home'] + 'genobox_bamsaturation.py --bams %s --subsample --sample stats --blocks 20 | cat -' % (bam)
   c2 = 'R-2.12 --vanilla stats/stats_Map.txt stats/stats_Map.txt stats/%s < %sgenobox_bamsaturation_plot.R' % (bamf, paths['genobox_home'])
   return [c1, c2]


def start_bamstats(args, bam, partition, logger, wait=True):
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
   cpuE = 'nodes=1:ppn=1,mem=7gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   cpuUV = 'procs=1,mem=%i,walltime=172800,flags=sharedmem'
   
   # create calls
   if args.mapdamage:
      mapdamage_calls = mapdamapge(bam, args.fa)
   else:
      flagstat_calls = sam_flagstat(bam)
      coverage_calls = bed_genomeCov(bam)
      plotcoverage_calls = plot_coverage(bam)
      avgdepth_calls = python_avgdepth(bam)
      saturation_calls = get_saturation(bam)
   
   
   # submit jobs
   print "Submitting jobs"
   if args.mapdamage:
      mapdamage_moab = Moab(mapdamage_calls, logfile=logger, runname='run_genobox_mapdamage', queue=args.queue, cpu=cpuA, partition=partition)
   else:
      flagstat_moab = Moab(flagstat_calls, logfile=logger, runname='run_genobox_flagstat', queue=args.queue, cpu=cpuC, partition=partition)
      coverage_moab = Moab(coverage_calls, logfile=logger, runname='run_genobox_coverage', queue=args.queue, cpu=cpuC, partition=partition)
      plotcoverage_moab = Moab(plotcoverage_calls, logfile=logger, runname='run_genobox_plotcoverage', queue=args.queue, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=coverage_moab.ids, partition=partition)
      avgdepth_moab = Moab(avgdepth_calls, logfile=logger, runname='run_genobox_avgdepth', queue=args.queue, cpu=cpuE, partition=partition)
      #saturation_moab = Moab(saturation_calls, logfile=logger, runname='run_genobox_saturation', queue=args.queue, cpu=cpuE, partition=partition)
   
   # release jobs
   print "Releasing jobs"
   if args.mapdamage:
      mapdamage_moab.release()
   else:
      flagstat_moab.release()
      coverage_moab.release()
      plotcoverage_moab.release()
      avgdepth_moab.release()
      #saturation_moab.release()
   
   # wait for jobs to finish
   if wait:
      print "Waiting for jobs to finish ..."
      if args.mapdamage:
         semaphore_ids = mapdamage_moab.ids
      else:
         semaphore_ids = flagstat_moab.ids + coverage_moab.ids + plotcoverage_moab.ids + avgdepth_moab.ids
      
      s = Semaphore(semaphore_ids, home, 'bam_stats', args.queue, 20, 86400) 
      s.wait()
      print "--------------------------------------"
   else:
      print "Jobs running, continuing"
      print "--------------------------------------"
   
   