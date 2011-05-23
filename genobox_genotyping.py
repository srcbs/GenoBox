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
   import os
   
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
   
   # get bamindex calls
   bamindex_calls = bam_index(bam)
   
   # get mpileup calls
   (mpileup_calls, bcffiles) = mpileup(bam, chr, fa, prior, pp)
   
   # get bcf combine calls
   bcfcombine_calls = bcf_combine(bcffiles, o)
   
   # get bcf index call
   bcfindex_calls = bcf_index(o)
   
   # submit jobs #
   print "Submitting jobs"
   bamindex_ids = genobox_modules.submitjob(bamindex_calls, home, paths, logger, 'run_genobox_bamindex', queue, cpuC, False)
   mpileup_ids = genobox_modules.submitjob(mpileup_calls, home, paths, logger, 'run_genobox_mpileup', queue, cpuF, True, 'expand', len(mpileup_calls), True, *bamindex_ids)
   bcfcombine_ids = genobox_modules.submitjob(bcfcombine_calls, home, paths, logger, 'run_genobox_bcfcombine', queue, cpuC, True, 'conc', len(mpileup_calls), True, *mpileup_ids)
   bcfindex_ids = genobox_modules.submitjob(bcfindex_calls, home, paths, logger, 'run_genobox_bcfindex', queue, cpuC, True, 'one2one', 1, True, *bcfcombine_ids)
   
   # release jobs #
   allids = []
   allids.extend(bamindex_ids) ; allids.extend(mpileup_ids) ; allids.extend(bcfcombine_ids) ; allids.extend(bcfindex_ids)
   releasemsg = genobox_modules.releasejobs(allids)
   
   # semaphore
   print "Waiting for jobs to finish ..."
   genobox_modules.wait_semaphore(bcfindex_ids, home, 'genotyping', queue, 20, 2*86400)
   print "--------------------------------------"
   
   # remove temporary files
   genobox_modules.rm_files(bcffiles)
   
   # return output bcf
   return o
