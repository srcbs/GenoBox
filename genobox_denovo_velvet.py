#!/panvol1/simon/bin/python

import argparse
import os
import logging

def required_nargs_abspath(min,max):
   '''Enforces input to nargs to be between min and max long'''
   class RequiredInterval(argparse.Action):
      def __call__(self, parser, args, value, option_string=None):
         if not min<=len(value)<=max:
            msg='argument "{f}" requires between {min} and {max} arguments'.format(
               f=self.dest,min=min,max=max)
            raise argparse.ArgumentTypeError(msg)
         setattr(args, self.dest, value)
         import os
         if type(value) == str:
            f_abs = os.path.abspath(value)
            setattr(args, self.dest, f_abs)
         elif type(value) == list:
            new_list = [value[0]]
            for f in value[1:]:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, value)
   return RequiredInterval

def required_nargs(min,max):
   '''Enforces input to nargs to be between min and max long'''
   class RequiredInterval(argparse.Action):
      def __call__(self, parser, args, value, option_string=None):
         if not min<=len(value)<=max:
            msg='argument "{f}" requires between {min} and {max} arguments'.format(
               f=self.dest,min=min,max=max)
            raise argparse.ArgumentTypeError(msg)
         setattr(args, self.dest, value)
   return RequiredInterval

def set_kmersizes(args, no_reads=10000, step=4):
   '''Automatically sets kmersizes from avg read lengths'''
   
   from Bio.SeqIO.QualityIO import FastqGeneralIterator
   from Bio import SeqIO
   import sys
   
   def average(values):
      '''Computes the arithmetic mean of a list of numbers.'''
      return sum(values, 0.0) / len(values)
   
   def get_readlengths(f, no_reads):
      '''Get avg. readlengths of no_reads'''
            
      def get_filetype(f):
         '''Return filetype (fasta/fastq)'''
         
         inhandle = open(f, "r")
         line = inhandle.readline()
         if line.startswith(">"):
            out = 'fasta'
         elif line.startswith("@"):
            out = 'fastq'
         else:
            raise ValueError('Input must be fasta or fastq')
         return out
      
      # main
      fh = open(f, 'r')
      count = 0
      L = []
      if get_filetype(f) == 'fastq':
         for (title, sequence, quality) in FastqGeneralIterator(fh):
            if count < no_reads:
               L.append(len(sequence))
               count = count + 1
            else:
               return average(L)
      elif get_filetype(f) == 'fasta':
         for rec in SeqIO.parse(fh, 'fasta'):
            if count < no_reads:
               L.append(len(rec.seq))
               count = count + 1
            else:
               return average(L)
      else:
         raise ValueError('Input must be fasta or fastq')
   
   def floor_to_odd(number):
      number = int(number)
      if number%2==0:
         return number - 1
      else:
         return number
   
   # main
   files = []
   if args.short: files.extend(args.short)
   if args.short2: files.extend(args.short2)
   if args.shortPaired: files.extend(args.shortPaired)
   if args.shortPaired2: files.extend(args.shortPaired2)
   
   # remove any occurences of fasta, fastq, fastq.gz, fasta.gz
   files = filter (lambda a: a != 'fasta', files)
   files = filter (lambda a: a != 'fastq', files)
   files = filter (lambda a: a != 'fasta.gz', files)
   files = filter (lambda a: a != 'fastq.gz', files)
   
   # get aerage lengths
   avg_lengths = []
   for f in files:
      avg_f = get_readlengths(f, no_reads)
      avg_lengths.append(avg_f)
   
   # get average of files and set ksizes from this
   avg = average(avg_lengths)
   ksizes = [str(floor_to_odd(avg/3)), str(floor_to_odd(avg/3*2.5)), '4']
   # change min ksizes if it is smaller than 15. Memory requirements becomes very large at k < 15
   if int(ksizes[0]) < 15: ksizes[0] = '15'
   sys.stderr.write('Ksizes set to (min, max, step) %s\n' % ' '.join(ksizes))
   return ksizes


def illumina_trim(args, min_length, min_baseq, min_avgq, min_adaptor_match, keep_n):
   '''Create single end trim calls'''
   
   import os
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   cmd = '%sgenobox_trim_h.py' % (paths['genobox_home'])
   calls = []
   if args.short:
      if args.short[0] == 'fastq':
         outfiles_short = []   
         for i,f in enumerate(args.short):
            if i == 0: continue
            outfile_short = 'trimmed/' + os.path.split(f)[1] + '.trim.fq'
            outfiles_short.append(outfile_short)
            arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s ' % (f, min_length,
               min_baseq, min_avgq, min_adaptor_match, outfile_short)
            if keep_n: arg = arg + ' --keep_n'
            calls.append(cmd+arg)
         args.short[1:] = outfiles_short
   
   if args.short2:
      if args.short2[0] == 'fastq':
         outfiles_short2 = []   
         for i,f in enumerate(args.short2):
            if i == 0: continue
            outfile_short2 = 'trimmed/' + os.path.split(f)[1] + '.trim.fq'
            outfiles_short2.append(outfile_short2)
            arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i  --min_adaptor_match %i --o %s ' % (f, min_length,
               min_baseq, min_avgq, min_adaptor_match, outfile_short2)
            if keep_n: arg = arg + ' --keep_n'
            calls.append(cmd+arg)
         args.short2[1:] = outfiles_short2
   
   if args.shortPaired and args.shortPaired[0] == 'fastq':
      outfiles_shortPaired = []
      if len(args.shortPaired) == 3:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired[1])[1] + '.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired[2])[1] + '.trim.fq'
         outfiles_shortPaired.append(outfile_pe1)
         outfiles_shortPaired.append(outfile_pe2)
         arg = ' --i %s %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired[1], args.shortPaired[2], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      elif len(args.shortPaired) == 2:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired[1])[1] + '_1.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired[1])[1] + '_2.trim.fq'
         outfiles_shortPaired.append(outfile_pe1)
         outfiles_shortPaired.append(outfile_pe2)
         arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired[1], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      else:
         raise ValueError('Length of input to shortPaired is not correct')
      
      if keep_n: arg = arg + ' --keep_n'
      calls.append(cmd+arg)
      args.shortPaired[1:] = outfiles_shortPaired
   
   if args.shortPaired2 and args.shortPaired2[0] == 'fastq':
      outfiles_shortPaired2 = []
      if len(args.shortPaired2) == 3:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired2[1])[1] + '.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired2[2])[1] + '.trim.fq'
         outfiles_shortPaired2.append(outfile_pe1)
         outfiles_shortPaired2.append(outfile_pe2)
         arg = ' --i %s %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired2[1], args.shortPaired2[2], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      elif len(args.shortPaired2) == 2:
         outfile_pe1 = 'trimmed/' + os.path.split(args.shortPaired2[1])[1] + '_1.trim.fq'
         outfile_pe2 = 'trimmed/' + os.path.split(args.shortPaired2[1])[1] + '_2.trim.fq'
         outfiles_shortPaired2.append(outfile_pe1)
         outfiles_shortPaired2.append(outfile_pe2)
         arg = ' --i %s --min_length %i --min_baseq %i --min_avgq %i --min_adaptor_match %i --o %s %s' % (args.shortPaired2[1], min_length, min_baseq, min_avgq,  min_adaptor_match, outfile_pe1, outfile_pe2)
      else:
         raise ValueError('Length of input to shortPaired2 is not correct')
      if keep_n: arg = arg + ' --keep_n'
      calls.append(cmd+arg)
      args.shortPaired2[1:] = outfiles_shortPaired2
   
   if len(calls) > 0:
      if not os.path.exists('trimmed'):
         os.makedirs('trimmed')
   return calls

def interleave(i, sample):
   '''Identify paired end read files that should be interleaved'''
   
   import os
   
   def merge(reads, format, interleaved):
      '''Perform merging'''
      
      import genobox_modules
      paths = genobox_modules.setSystem()
      
      # shuffle <file1> <file2> <out>
      if format.find('fastq') > -1: cmd = '%sshuffleSequences_fastq.pl %s %s %s' % (paths['velvet_home'], reads[0], reads[1], interleaved)
      if format.find('fasta') > -1: cmd = '%sshuffleSequences_fasta.pl %s %s %s' % (paths['velvet_home'], reads[0], reads[1], interleaved)
      return cmd
   
   # check if format can be merged (fasta, fastq)
   if i:
      if i[0].find('fasta') > -1 or i[0].find('fastq') > -1:
         # check if merging is needed
         if len(i) > 2:
            format = i[0]
            reads = i[1:3]
            interleaved =  os.path.split(reads[0])[1] + '.interleaved'
            #interleaved = reads[0] + '.interleaved'
            # set new format
            call = merge(reads, format, interleaved)
            return (call, [format, interleaved])
         else:
            return (None, i)
      else:
         return (None, i)
   else:
      return (None, i)
        

def create_velveth_calls(args):
   '''Return velveth calls'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   cmd = '%svelveth' % paths['velvet_home']
   
   # create calls, outpath, ksizes, format, readtypes, reads
   velveth_calls = []
   if len(args.ksizes) == 1:
      arg = ' %s %s ' % (args.outpath, args.ksizes[0])
      if args.short: arg = arg + ' -short -%s %s' % (args.short[0], ' '.join(args.short[1:]))
      if args.short2: arg = arg + ' -short2 -%s %s' % (args.short2[0], ' '.join(args.short2[1:]))
      if args.shortPaired: arg = arg + ' -shortPaired -%s %s' % (args.shortPaired[0], args.shortPaired[1])
      if args.shortPaired2: arg = arg + ' -shortPaired2 -%s %s' % (args.shortPaired2[0], args.shortPaired2[1])
      if args.long: arg = arg + ' -long -%s %s' % (args.long[0], ' '.join(args.long[1:]))
      if args.longPaired: arg = arg + ' -longPaired -%s %s' % (args.longPaired[0], args.longPaired[1])
      if args.add_velveth: arg = arg + ' %s' % args.add_velveth
      call = cmd + arg
      velveth_calls.append(call)
   elif len(args.ksizes) >= 2 and len(args.ksizes) <= 3:
      if len(args.ksizes) == 2:
         step = 2
      elif len(args.ksizes) == 3:
         step = args.ksizes[2]
       
      # create calls, outpath, ksizes, format, readtypes, reads
      for k in range(int(args.ksizes[0]), int(args.ksizes[1]), int(step)):
         arg = ' %s_%s %s ' % (args.outpath, k, k)
         if args.short: arg = arg + ' -short -%s %s' % (args.short[0], ' '.join(args.short[1:]))
         if args.short2: arg = arg + ' -short2 -%s %s' % (args.short2[0], ' '.join(args.short2[1:]))
         if args.shortPaired: arg = arg + ' -shortPaired -%s %s' % (args.shortPaired[0], args.shortPaired[1])
         if args.shortPaired2: arg = arg + ' -shortPaired2 -%s %s' % (args.shortPaired2[0], args.shortPaired2[1])
         if args.long: arg = arg + ' -long -%s %s' % (args.long[0], ' '.join(args.long[1:]))
         if args.longPaired: arg = arg + ' -longPaired -%s %s' % (args.longPaired[0], args.longPaired[1])
         if args.add_velveth: arg = arg + ' %s' % args.add_velveth
         call = cmd + arg
         velveth_calls.append(call)
   else:
      raise ValueError('ksizes must be one value giving ksize, two values giving lower and upper limit (step will be 2) or three values giving lower limit, upper limit and step')
   return velveth_calls   

def create_velvetg_calls(args):
   '''Return velvetg calls'''
   
   import genobox_modules   
   paths = genobox_modules.setSystem()
   
   # create cmd
   cmds = []
   if len(args.ksizes) == 1:
      cmd = '%svelvetg %s' % (paths['velvet_home'], args.outpath)
      cmds.append(cmd)
   elif len(args.ksizes) >= 2 and len(args.ksizes) <= 3:
      if len(args.ksizes) == 2:
         step = 2
      elif len(args.ksizes) == 3:
         step = args.ksizes[2]
      
      for k in range(int(args.ksizes[0]), int(args.ksizes[1]), int(step)):
         cmd = '%svelvetg %s_%s' % (paths['velvet_home'], args.outpath, k)
         cmds.append(cmd)
   
   # create arg: cov_cut, exp_cov, ins_length, add_velvetg
   velvetg_calls = []
   # add other parameters
   for i in range(len(cmds)):
      arg = ' -min_contig_lgth %i' % args.min_contig_lgth
      if args.cov_cut: arg = arg + ' -cov_cut %f' % args.cov_cut
      if args.exp_cov != "None": arg = arg + ' -exp_cov %s' % args.exp_cov
      if args.ins_length: arg = arg + ' -ins_length %i' % args.ins_length
      if args.add_velvetg: arg = arg + ' %s' % args.add_velvetg
      velvetg_calls.append(cmds[i]+arg)
    
   return velvetg_calls

def get_best_assembly(args):
   '''Identify the best assembly from several k-mers'''
   
   # read in stats.txt files for each assembly. Calc sum of contigs and N50.
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   cmd = '%sR-2.12 --vanilla ' % paths['R_home']
   
   # set argument
   if len(args.ksizes) == 1:
      arg = ' %s %s' % (args.outpath, args.ksizes[0])
   elif len(args.ksizes) >= 2:
      if len(args.ksizes) == 2:
         step = 2
      elif len(args.ksizes) == 3:
         step = args.ksizes[2]
      
      arg_list = []
      for k in range(int(args.ksizes[0]), int(args.ksizes[1]), int(step)):
         out = '%s_%s/stats.txt %s' % (args.outpath, k, k)
         arg_list.append(out)
      arg = ' '.join(arg_list)
   
   call = [cmd + arg + ' < %sgenobox_denovo_velvet_parse.R' % (paths['genobox_home'])]
   return call

def accept_assembly(args):
   '''Parse best assembly and remove other assemblies'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   call = '%sgenobox_denovo_velvet_accept.py %s' % (paths['genobox_home'], args.outpath)
   return [call]

def clean():
   '''Clean sample directory'''
   
   import genobox_modules
   paths = genobox_modules.setSystem()
   
   call = '%sgenobox_denovo_velvet_clean.py' % (paths['genobox_home'])
   return [call]

def start_assembly(args, logger):
   '''Start assembly'''
   
   import genobox_modules
   from genobox_classes import Moab
   from genobox_classes import Semaphore   
   import os
   
   # set queueing
   paths = genobox_modules.setSystem()
   home = os.getcwd()
   cpuV = 'nodes=1:ppn=%i,mem=%s,walltime=172800' % (args.n, args.m)
   cpuA = 'nodes=1:ppn=1,mem=512mb,walltime=172800'
   cpuC = 'nodes=1:ppn=1,mem=2gb,walltime=172800'
   cpuE = 'nodes=1:ppn=1,mem=5gb,walltime=172800'
   cpuF = 'nodes=1:ppn=2,mem=2gb,walltime=172800'
   cpuB = 'nodes=1:ppn=16,mem=10gb,walltime=172800'
   
   # set kmersizes (if auto)
   if args.ksizes == ['auto']:
      args.ksizes = set_kmersizes(args)   
   
   # trimming calls
   if args.trim: illuminatrim_calls = illumina_trim(args, int(args.ksizes[0]), 15, 20, 15, False)
   
   # checking if files needs to be interleaved
   interleave_dict = {}    
   interleave_dict['shortPaired'] = interleave(args.shortPaired, args.sample)[0] ; args.shortPaired = interleave(args.shortPaired, args.sample)[1]
   interleave_dict['shortPaired2'] = interleave(args.shortPaired2, args.sample)[0] ; args.shortPaired2 = interleave(args.shortPaired2, args.sample)[1]
   interleave_dict['longPaired'] = interleave(args.longPaired, args.sample)[0] ; args.longPaired = interleave(args.longPaired, args.sample)[1]
   
   # interleave calls
   interleave_calls = []
   for key,value in interleave_dict.items():
      if value:
         interleave_calls.append(value)
   
   # velvet calls
   velveth_calls = create_velveth_calls(args)
   velvetg_calls = create_velvetg_calls(args)
   
   # velvet parse calls
   velvetparse_calls = get_best_assembly(args)
   velvetaccept_calls = accept_assembly(args)
   velvetclean_calls = clean()
   
   # set environment variable:
   env_var = 'OMP_NUM_THREADS=%i' % int(args.n - 1)
   
   # submit and release jobs
   print "Submitting jobs"
   # if trimming is needed
   if args.trim:
      illuminatrim_moab = Moab(illuminatrim_calls, logfile=logger, runname='run_genobox_trim', queue=args.queue, cpu=cpuF)
      # if no interleaving is needed
      if len(interleave_calls) == 0:
         velveth_moab = Moab(velveth_calls, logfile=logger, runname='run_genobox_velveth', queue=args.queue, cpu=cpuV, depend=True, depend_type='all', depend_val=[1], depend_ids=illuminatrim_moab.ids, env=env_var)
         velvetg_moab = Moab(velvetg_calls, logfile=logger, runname='run_genobox_velvetg', queue=args.queue, cpu=cpuV, depend=True, depend_type='one2one', depend_val=[1], depend_ids=velveth_moab.ids)
      # if interleaving is needed
      else:
         interleave_moab = Moab(interleave_calls, logfile=logger, runname='run_genobox_interleave', queue=args.queue, cpu=cpuF, depend=True, depend_type='all', depend_val=[1], depend_ids=illuminatrim_moab.ids)
         velveth_moab = Moab(velveth_calls, logfile=logger, runname='run_genobox_velveth', queue=args.queue, cpu=cpuV, depend=True, depend_type='all', depend_val=[1], depend_ids=interleave_moab.ids, env=env_var)
         velvetg_moab = Moab(velvetg_calls, logfile=logger, runname='run_genobox_velvetg', queue=args.queue, cpu=cpuV, depend=True, depend_type='one2one', depend_val=[1], depend_ids=velveth_moab.ids)
   # if no trimming
   else:
      # if no interleaving is needed
      if len(interleave_calls) == 0:
         velveth_moab = Moab(velveth_calls, logfile=logger, runname='run_genobox_velveth', queue=args.queue, cpu=cpuV, env=env_var)
         velvetg_moab = Moab(velvetg_calls, logfile=logger, runname='run_genobox_velvetg', queue=args.queue, cpu=cpuV, depend=True, depend_type='one2one', depend_val=[1], depend_ids=velveth_moab.ids)
      # if interleaving is needed
      else:
         interleave_moab = Moab(interleave_calls, logfile=logger, runname='run_genobox_interleave', queue=args.queue, cpu=cpuF)
         velveth_moab = Moab(velveth_calls, logfile=logger, runname='run_genobox_velveth', queue=args.queue, cpu=cpuV, depend=True, depend_type='all', depend_val=[1], depend_ids=interleave_moab.ids, env=env_var)
         velvetg_moab = Moab(velvetg_calls, logfile=logger, runname='run_genobox_velvetg', queue=args.queue, cpu=cpuV, depend=True, depend_type='one2one', depend_val=[1], depend_ids=velveth_moab.ids)
   
   # submit job for velvetparse if more than one ksize was chosen
   if len(args.ksizes) > 1:
      velvetparse_moab = Moab(velvetparse_calls, logfile=logger, runname='run_genobox_velvetparse', queue=args.queue, cpu=cpuA, depend=True, depend_type='conc', depend_val=[len(velvetg_calls)], depend_ids=velvetg_moab.ids)
      velvetaccept_moab = Moab(velvetaccept_calls, logfile=logger, runname='run_genobox_velvetaccept', queue=args.queue, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=velvetparse_moab.ids) 
      velvetclean_moab = Moab(velvetclean_calls, logfile=logger, runname='run_genobox_velvetclean', queue=args.queue, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=velvetaccept_moab.ids)
   else:
      velvetclean_moab = Moab(velvetclean_calls, logfile=logger, runname='run_genobox_velvetclean', queue=args.queue, cpu=cpuA, depend=True, depend_type='one2one', depend_val=[1], depend_ids=velvetg_moab.ids)
   
   # release jobs
   print "Releasing jobs"
   if args.trim and len(illuminatrim_calls) > 0: illuminatrim_moab.release()
   if len(interleave_calls) > 0: interleave_moab.release()
   velveth_moab.release()
   velvetg_moab.release()
   if len(args.ksizes) > 1: 
      velvetparse_moab.release()
      velvetaccept_moab.release()
   velvetclean_moab.release()
   
   # semaphore (consensus is currently not waited for)
   print "Waiting for jobs to finish ..."
   s = Semaphore(velvetclean_moab.ids, home, 'velvet', args.queue, 20, 2*86400)
   s.wait()
   print "--------------------------------------"
   

if __name__ == '__main__':
   # create the parser
   parser = argparse.ArgumentParser(prog='genobox_denovo_velvet.py', formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog,max_help_position=50, width=110), usage='%(prog)s [options]', description='''
   Run Velvet denovo assembly. Read input is given by eg. --short <format> <reads>
   Format can be: fasta, fastq, raw, fasta.gz, fastq.gz, raw.gz, sam, bam
   add_velvetg/add_velveth has to be added with quotes, eg: add_velvetg "-very_clean yes
   if mate_pairs are used (in eg. shortPaired2) one should add: add_velvetg "-shortMatePaired2 yes -ins_length2 <MP size>"\n''')

   # add the arguments
   parser.add_argument('--short', help='input read format and short reads', nargs='+', action=required_nargs_abspath(0,3))
   parser.add_argument('--shortPaired', help='input read format and short paired reads', nargs='+', action=required_nargs_abspath(0,3))
   parser.add_argument('--short2', help='input read format and short2 reads', nargs='+', action=required_nargs_abspath(0,3))
   parser.add_argument('--shortPaired2', help='input read format and short paired2 reads', nargs='+', action=required_nargs_abspath(0,3))
   parser.add_argument('--long', help='input read format and long reads', nargs='+', action=required_nargs_abspath(0,3))
   parser.add_argument('--longPaired', help='input read format and long paired reads', nargs='+', action=required_nargs_abspath(0,3))
   parser.add_argument('--ksizes', help='kmers to run assemblies for (single (m) or m M s (min, max, step)) [auto]', nargs='+', default=['auto'])
   parser.add_argument('--sample', help='name of run and output directory [velvet_assembly]', default='velvet_assembly')
   parser.add_argument('--outpath', help='name of assembly output dir [assembly]', default='assembly')
   parser.add_argument('--trim', help='should input files be trimmed (illumina only) [False]', default=False, action='store_true')
   parser.add_argument('--min_contig_lgth', help='mininum length to report contig [100]', default=100, type=int)
   parser.add_argument('--cov_cut', help='coverage cutoff for removal of low coverage (float) [None]', default=None, type=float)
   parser.add_argument('--exp_cov', help='Expected mean coverage (None, float, auto) [auto]', default='auto')
   parser.add_argument('--ins_length', help='insert size (reads included) [None]', default=None, type=int)
   parser.add_argument('--add_velveth', help='additional parameters to velveth', default=None)
   parser.add_argument('--add_velvetg', help='additional parameters to velvetg', default=None)
   parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
   parser.add_argument('--m', help='memory needed for assembly [2gb]', default='2gb')
   parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cbs]', default='cbs')
   parser.add_argument('--log', help='log level [info]', default='info')

   args = parser.parse_args()
   #args = parser.parse_args('--short fastq test_1.fq test_2.fq --ksizes 33 49 4 --outpath test'.split())
   #args = parser.parse_args('--short fastq.gz Kleb-10-213361_2.interleaved.fastq.test.gz --ksizes 41 55 4 --sample Kleb'.split())

   # add_velveth and add_velvetg works from commandline, eg:
   # genobox_denovo_velvet.py --short fastq.gz interleaved.fastq.gz --ksizes 33 --sample Kleb --add_velvetg "-very_clean yes"
   
   # change to sample dir if set
   if args.sample:
      if not os.path.exists(args.sample):
         os.makedirs(args.sample)
      os.chmod(args.sample, 0777)
      os.chdir(args.sample)
   else:
      pass
   
   # create log dir
   if not os.path.exists('log'):
      os.makedirs('log')
   
   # start logging
   logger = logging.getLogger('genobox_denovo_velvet.py')
   hdlr = logging.FileHandler('genobox_denovo_velvet.log')
   formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
   hdlr.setFormatter(formatter)
   logger.addHandler(hdlr) 
   if args.log == 'info':
      logger.setLevel(logging.INFO)
   
   # start assembly
   start_assembly(args, logger)
