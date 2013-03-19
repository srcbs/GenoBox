#!/panvol1/simon/bin/python2.7

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

def unified_genotyper(bam, genome, fa, dbsnp,  call_conf, call_emit, output_mode):
   '''Perform genotyping on bam-file using GATK unified genotyper'''
   
   import genobox_modules
   import os
   
   paths = genobox_modules.setSystem()
   gatk_cmd = paths['GATK_home'] + 'GenomeAnalysisTK.jar'
   java_cmd = 'java -Djava.io.tmpdir=/panvol1/simon/tmp/ -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms3000m -Xmx3000m -jar '
   cmd = java_cmd + gatk_cmd
   
   calls = []
   outfiles = []
   basename = os.path.split(bam)[1]
   
   chrs = get_genome(genome)
   for c in chrs:
      outfile = 'genotyping/%s.%s.raw.vcf.gz' % (basename.replace('.bam', ''), c[2])
      logfile = 'log/run_unified_genotyper.%s.%s.log' % (basename.replace('.bam', ''), c[2])
      outfiles.append(outfile)
      arg = ' -T UnifiedGenotyper -R %s -I %s -o /dev/stdout -log %s -stand_call_conf %f -stand_emit_conf %f -L %s -baq CALCULATE_AS_NECESSARY --num_threads 1 -glm BOTH --output_mode %s ' % (fa, bam, logfile, call_conf, call_emit, c[2], output_mode)
      if dbsnp: arg = arg + '--dbsnp %s ' % dbsnp
      arg = arg + ''' | perl -ne 'if ($_ =~ m/^INFO/ or $_ =~ m/^WARN/) {} else {print $_}' | gzip -c - > %s''' % outfile
      calls.append(cmd+arg)
   
   return (calls, outfiles)


def start_genotyping_gatk(bam, genome, fa, dbsnp, call_conf, call_emit, output_mode, queue, sample, partition, logger):
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
   cpuE = 'nodes=1:ppn=1,mem=3gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   
   # create calls
   bamindex_calls = bam_index(bam)
   (gatk_calls, vcffiles) = unified_genotyper(bam, genome, fa, dbsnp, call_conf, call_emit, output_mode)
   
   # submit jobs #
   print "Submitting jobs"   
   bamindex_moab = Moab(bamindex_calls, logfile=logger, runname='run_genobox_bamindex', queue=queue, cpu=cpuC, partition=partition)
   gatk_moab = Moab(gatk_calls, logfile=logger, runname='run_genobox_genotyping_gatk', queue=queue, cpu=cpuE, depend=True, depend_type='expand', depend_val=[len(gatk_calls)], depend_ids=bamindex_moab.ids, partition=partition)
   
   # release jobs #
   print "Releasing jobs"
   #bamindex_moab.release()
   #gatk_moab.release()
      
   # semaphore (consensus is currently not waited for)
   print "Waiting for jobs to finish ..."
   s = Semaphore(gatk_moab.ids, home, 'genotyping', queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   
   # remove temporary files
   #genobox_modules.rm_files(bcffiles)
   
   # return output variant files
   return vcffiles
