#!/usr/bin/python

# functions for use in the genobox pipeline #


################################

###      LIBRARY FILES       ###

################################

# format is tab separated:
# ID, Data, Mapq, Read group information
# eg: ID Data	MAPQ	LB	PL	SM	CN

# test
#se = ['SRR002081se.recal.fastq', 'SRR002082se.recal.fastq']
#pe1 = ['SRR002137pe_1.recal.fastq', 'SRR002137pe_1.recal.fastq']
#pe2 = ['SRR002137pe_2.recal.fastq', 'SRR002137pe_2.recal.fastq']
#input = se + pe1 + pe2
#sample = 'NA12891'
#mapq = [30]
#libs = ['libSe']

def read_library(library_file):
   '''Reads in library file and returns dict with ID as keys'''
   
   fh = open(library_file, 'r')
   header = fh.readline().rstrip().split('\t')
   
   # fill dict with data
   id_dict = dict()
   for line in fh:
      fields = line.rstrip().split('\t')
      if len(header) == len(fields):
         id_dict[fields[0]] = fields[0:]
      else:
         raise ValueError('Line \"%s\" do not have same number of fields as header')
   
   id_dict['header'] = header
   return id_dict

def library_from_input(input, sample='sample', mapq=[30], libs=['lib']):
   '''Create library dict from inputs'''
   
   import sys
   
   # check lengths of mapq and libs
   if len(mapq) < len(input):
      sys.stderr.write("warning: length of mapq is shorter than length of inputfiles, reusing %i\n" % mapq[0])   
      while len(mapq) < len(input):
         mapq.append(mapq[0])
   if len(libs) < len(input):
      sys.stderr.write("warning: length of libs is shorter than length of inputfiles, reusing %s\n" % libs[0])   
      while len(libs) < len(input):
         libs.append(libs[0])
   
   # create dict
   id_dict = dict()
   c = 0
   for i,f in enumerate(input):
      id = '%s_%i' % (sample, c + 1)
      id_dict[id] = [id, f, str(mapq[c]), libs[c], 'ILLUMINA', sample]
      c = c + 1
   
   id_dict['header'] = ['ID', 'Data', 'MAPQ', 'LB', 'PL', 'SM']
   
   # write to disk
   fh = open('libs.%s.txt' % sample, 'w')
   fh.write('%s\n' % '\t'.join(id_dict['header']))
   for key in sorted(id_dict.keys()):
      if key == 'header':
         continue
      fh.write('%s\n' % '\t'.join(id_dict[key]))
   return (id_dict, 'libs.%s.txt' % sample)

def update_libfile(libfile, key_col, new_col, val_dict, force=False):
   '''Update libfile with a new column of data'''
   
   import sys
   
   # read in current file
   fh = open(libfile, 'r')
   header = fh.readline().rstrip().split('\t')
   data = fh.read().split('\n')
   if data[-1] == '':
      data = data[:-1]
   fh.close()
   
   # check length of data and new val_dict
   if len(data) == len(val_dict):
      pass
   else:
      raise ValueError('Rows of libfile (%i) does not match length of new val_dict (%i)' % (len(data), len(val_dict)))
   
   # check if column already exist, if forced then remove old column
   if new_col in header:
      if force == False:
         sys.stderr.write('column %s already exists, %s not updated\n' % (new_col, libfile))
         return 
      else:
         n = header.index(new_col)
         header.remove(new_col)
         new_data = []
         for line in data:
            fields = line.split('\t')
            new_fields = fields[:n] + fields[(n+1):]
            new_data.append('\t'.join(new_fields))
         data = new_data
   
   # append data
   fh = open(libfile, 'w')
   header.append(new_col)
   fh.write('%s\n' % '\t'.join(header))
   for line in data:
      fields = line.split('\t')
      for key in val_dict.keys():
         if fields[header.index(key_col)] == key:
            new_line = '%s\t%s\n' % (line, val_dict[key])
            fh.write(new_line)
            break
   
   fh.close()
   return

def read_groups_from_libfile(index, libfile):
   '''Return read group from libfile with index as key'''
   
   allowedRGs = ['ID', 'CN', 'DS', 'DT', 'FO', 'KS', 'LB', 'PG', 'PI', 'PL', 'PU', 'SM']
   
   RG = dict()
   #RG = '\@RG\tID:foo\tSM:bar'
   header = libfile['header']
   for key,value in libfile.items():
      if key == 'header':
         continue
      
      currRG = ['@RG']
      for h in header:
         if h in allowedRGs:
            currRG.append('%s:%s' % (h, value[header.index(h)]))
      RG[value[header.index(index)]] = currRG
   return RG

# read library file 
def read_bam_libs(libfile):
   '''Read library file and return dict'''
   
   from collections import defaultdict
   
   fh = open(libfile, 'r')
   header = fh.readline().rstrip().split('\t')
   
   # check if bamfiles are written in libfile
   if 'BAM' in header:
      pass
   else:
      raise IndexError('There is no column \"BAM\" in the libfile (%s)' % libfile)
   
   bam2lib = {}
   lib2bam = defaultdict(list)
   for line in fh:
      fields = line.rstrip().split('\t')
      bam2lib[fields[header.index('BAM')]] = (fields[header.index('MAPQ')], fields[header.index('LB')])
      
      lib2bam[fields[header.index('LB')]].append(fields[header.index('BAM')])
   
   # unique on lib2bam
   for key,values in lib2bam.items():
      lib2bam[key] = unique(values)
   
   return (bam2lib, lib2bam)


##################################

###      FASTQ DETECTION       ###

##################################

def set_filetype(f):
   '''Detects filetype from fa, fq and sff input file'''
   
   inhandle = open(f, "r")
   line = inhandle.readline()
   if line.startswith(">"):
      out = 'fasta'
   elif line.startswith("@"):
      out = 'fastq'
   else:
      inhandle = open(f, "rb")
      line = inhandle.readline()
      if line.startswith(".sff"):
         out = 'sff'
      else:
         raise ValueError('Input must be fasta, fastq or sff')
   
   return out

def set_fqtype(f):
   '''Detects sanger or illumina format from fastq'''
   
   # Open fastq, convert ASCII to number, check if number is above or below certain thresholds to determine format
   
   from Bio.SeqIO.QualityIO import FastqGeneralIterator
   
   type = 'not determined'
   inhandle = open(f, 'r')
   for (title, sequence, quality) in FastqGeneralIterator(inhandle):
      qs = map(ord, quality)
      for q in qs:
         if q > 73:
            type = 'Illumina'
            break
         elif q < 59:
            type = 'Sanger'
            break
      if type != 'not determined':
         break
   
   if type == 'not determined':
      raise ValueError('Fastq format not identified, are you sure it is sanger/illumina?')
   else:
      return type


###################################

###      GENERAL FUNCTIONS      ###

###################################


def unique(seq):
   '''Return unique and ordered list'''
   
   # order preserving 
   checked = [] 
   for e in seq: 
      if e not in checked: 
         checked.append(e) 
   return checked

