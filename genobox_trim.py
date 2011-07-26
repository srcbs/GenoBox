#!/panvol1/simon/bin/python2.7

def single_trim(args):
   '''Create single end trim calls'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   cmd = '%sgenobox_trim_h.py' % (paths['genobox_home'])
   calls = []
   outfiles_se = []
   for i,f in enumerate(args.se):
      outfile_se = 'trimmed/' + os.path.split(f)[1] + '.trim.fq'
      outfiles_se.append(outfile_se)
      arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i --adaptors %s  --min_adaptor_match %i --o %s ' % (f, args.min_length,
            args.min_baseq, args.min_avgq, ' '.join(args.adaptors), args.min_adaptor_match, outfile_se)
      if args.keep_n: arg = arg + ' --keep_n'
      calls.append(cmd+arg)
   return (calls, outfiles_se)

def paired_trim(args):
   '''Create paired end trim calls'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   cmd = '%sgenobox_trim_h.py' % (paths['genobox_home'])
   calls = []
   outfiles_pe1 = []
   outfiles_pe2 = []
   for i,f in enumerate(args.pe1):
      outfile_pe1 = 'trimmed/' + os.path.split(args.pe1[i])[1] + '.trim.fq'
      outfile_pe2 = 'trimmed/' + os.path.split(args.pe2[i])[1] + '.trim.fq'
      outfiles_pe1.append(outfile_pe1)
      outfiles_pe2.append(outfile_pe2)
      arg = ' --i %s %s --min_length %i --min_baseq %i --min_avgq %i --adaptors %s --min_adaptor_match %i --o %s %s' % (args.pe1[i], args.pe2[i], args.min_length,
            args.min_baseq, args.min_avgq, ' '.join(args.adaptors), args.min_adaptor_match, outfile_pe1, outfile_pe2)
      if args.keep_n: arg = arg + ' --keep_n'
      calls.append(cmd+arg)
   return (calls, outfiles_pe1, outfiles_pe2)

def start_trim(args, logger):
   '''Start trimming from genobox.py'''
   
   import genobox_modules
   from genobox_classes import Moab, Semaphore
   import subprocess
   import os
   
   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
   cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   
   # create path
   if not os.path.exists('trimmed'):
      os.makedirs('trimmed')
   
   # create calls
   (single_calls, se_files) = single_trim(args)
   (paired_calls, pe1_files, pe2_files) = paired_trim(args)
   
   # submit jobs
   print "Submitting jobs"
   if args.se:
      single_moab = Moab(single_calls, logfile=logger, runname='run_genobox_trimse', queue=args.queue, cpu=cpuA)
   if args.pe1 and args.pe2:
      paired_moab = Moab(paired_calls, logfile=logger, runname='run_genobox_trimpe', queue=args.queue, cpu=cpuF)
      
   # release jobs
   print "Releasing jobs"
   if args.se:
      single_moab.release()
   if args.pe1 and args.pe2:
      paired_moab.release()
   
   # wait for jobs to finish
   print "Waiting for jobs to finish ..."
   semaphore_ids = []
   if args.se:
      semaphore_ids = semaphore_ids + single_moab.ids
   if args.pe1 and args.pe2:
      semaphore_ids = semaphore_ids + paired_moab.ids
   s = Semaphore(semaphore_ids, home, 'read_trimming', args.queue, 60, 86400)
   s.wait()
   print "--------------------------------------"
   
   # return trimmed files
   return (se_files, pe1_files, pe2_files)
