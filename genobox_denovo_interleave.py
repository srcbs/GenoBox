#!/panvol1/simon/bin/python

import argparse

def readline(fh, no_lines):
    '''Readline from open filehandle'''
    lines = []
    for i in range(no_lines):
        lines.append(fh.readline())
    try:
        lines.index('')
    except:
        return lines
    else:
        return None

def interleave(pe1, pe2, out, format, gzip):
   '''Interleave two paired end files to one file'''
   
   import os
   import subprocess
   
   # open handles
   if format.find("gz") > -1:
      fh_pe1 = os.popen('gzip -d %s' % pe1, 'r')
      fh_pe2 = os.popen('gzip -d %s' % pe2, 'r')
   else:
      fh_pe1 = open(pe1, 'r')
      fh_pe2 = open(pe2, 'r')
    
   
   # check format
   if format.find('fasta') > -1:
      no_lines = 2
   elif format.find('fastq') > -1:
      no_lines = 4
   else:
      raise ValueError('fasta or fastq not specified in format')
   
   # start write out
   if gzip:
      process = subprocess.Popen('gzip -1 > %s' % out, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
      while 1:
         lines1 = readline(fh_pe1, no_lines)
         lines2 = readline(fh_pe2, no_lines)
         if not lines1 or not lines2:
            break
         for line in lines1:
            process.stdin.write(line)
         for line in lines2:
            process.stdin.write(line)
      
      # end Popen
      stdout, stderr = process.communicate()
   else:
      fh_out = open(out, 'w')
      while 1:
         lines1 = readline(fh_pe1, no_lines)
         lines2 = readline(fh_pe2, no_lines)
         if not lines1 or not lines2:
            break
         for line in lines1:
            fh_out.write(line)
         for line in lines2:
            fh_out.write(line)
      

if __name__ == '__main__':

   parser = argparse.ArgumentParser(prog='genobox_denovo_interleave.py',
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=25, width=140),
                                 description='''Mix reads from two pe-files and create an interleaved zipped file''', 
                                 usage='%(prog)s [options]')
   
   # add the arguments
   parser.add_argument('i', help='paired end files', nargs=2)
   parser.add_argument('format', help='input format of file (fasta, fastq, fasta.gz, fastq.gz)')
   parser.add_argument('out', help='output interleaved file')
   parser.add_argument('--gzip', help='should output be gzipped [False]', default=False, action='store_true')
   
   args = parser.parse_args()
   #args = parser.parse_args('Kleb-10-213361_2_1_sequence.trim.fq Kleb-10-213361_2_2_sequence.trim.fq fastq Kleb-10-213361_2.interleaved.fastq.test.gz'.split())
   interleave(args.i[0], args.i[1], args.out, args.format, args.gzip)
