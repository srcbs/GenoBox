#!/usr/bin/python

def bam_index(bam):
   '''Index bam-file'''
   
   import genobox_modules
   import os.path
   paths = genobox_modules.setSystem()
   
   # skip index creation if it already exists
   calls = []
   if not os.path.isfile(bam+'.bai'):
      cmd = paths['samtools_home'] + 'samtools'
      arg = ' index %s' % (bam)
      calls.append(cmd+arg)
   else:
      calls.append('sleep 0.01')
   return calls

def get_genome(chr_file):
   '''Read chromosome file into list of list'''
   fh = open(chr_file, 'r')
   L = []
   for line in fh:
      line = line.rstrip()
      s = line.split('\t')
      L.append(s)
   return L

def mpileup(bam, chr_file, fa, prior, pp):
   '''Perform SNP calling on bam-file using samtools'''
   
   import genobox_modules
   import os
   
   paths = genobox_modules.setSystem()
   cmd = 'python2.7 ' + paths['genobox_home'] + 'genobox_mpileup.py'
   calls = []
   outfiles = []
   
   # if chromosome file is given
   if chr_file:
      chrs = get_genome(chr_file)
      for c in chrs:
         outfile = 'genotyping/tmp.' + c[2] + '.all.bcf'
         outfiles.append(outfile)         
         arg = ' --bam %s --chr \"%s\" --fa %s --prior %s --pp %f --o %s' % (bam, c[0], fa, prior, pp, outfile)
         calls.append(cmd+arg)
   else:
      tmpfile_name = os.path.split(bam)[1]
      outfile = 'genotyping/tmp.' + tmpfile_name + '.all.bcf'
      outfiles.append(outfile)
      arg = ' --bam %s --fa %s --prior %s --pp %f --o %s' % (bam, fa, prior, pp, outfile)
      calls.append(cmd+arg)
   return (calls, outfiles)

def bcf_combine(bcfs, outfile):
   '''Concatenate bcfs to a single bcf '''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   calls = []
   cmd = paths['samtools_home'] + 'bcftools'
   arg = ' cat %s > %s' % (' '.join(bcfs), outfile)
   calls.append(cmd+arg)
   return calls

def bcf_index(bcf):
   '''Index bcf file'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   calls = []
   cmd = paths['samtools_home'] + 'bcftools'
   arg = ' index %s' % (bcf)
   calls.append(cmd+arg)
   return calls

def start_genotyping(bam, chr, fa, prior, pp, queue, o, logger):
   '''Starts genotyping using samtools of input bam file'''
   
   import subprocess
   import genobox_modules
   from genobox_classes import Moab
   from genobox_classes import Semaphore   
   import os
   
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
   
   # create calls
   bamindex_calls = bam_index(bam)
   (mpileup_calls, bcffiles) = mpileup(bam, chr, fa, prior, pp)
   bcfcombine_calls = bcf_combine(bcffiles, o)
   bcfindex_calls = bcf_index(o)
   
   # submit jobs #
   print "Submitting jobs"   
   bamindex_moab = Moab(bamindex_calls, logfile=logger, runname='run_genobox_bamindex', queue=queue, cpu=cpuC)
   mpileup_moab = Moab(mpileup_calls, logfile=logger, runname='run_genobox_mpileup', queue=queue, cpu=cpuF, depend=True, depend_type='expand', depend_val=[len(mpileup_calls)], depend_ids=bamindex_moab.ids)
   bcfcombine_moab = Moab(bcfcombine_calls, logfile=logger, runname='run_genobox_bcfcombine', queue=queue, cpu=cpuC, depend=True, depend_type='conc', depend_val=[len(mpileup_calls)], depend_ids=mpileup_moab.ids)
   bcfindex_moab = Moab(bcfindex_calls, logfile=logger, runname='run_genobox_bcfindex', queue=queue, cpu=cpuC, depend=True, depend_type='one2one', depend_val=[1], depend_ids=bcfcombine_moab.ids)
      
   # release jobs #
   print "Releasing jobs"
   bamindex_moab.release()
   mpileup_moab.release()
   bcfcombine_moab.release()
   bcfindex_moab.release()
      
   # semaphore
   print "Waiting for jobs to finish ..."
   s = Semaphore(bcfindex_moab.ids, home, 'genotyping', queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   
   # remove temporary files
   genobox_modules.rm_files(bcffiles)
   
   # return output bcf
   return o
