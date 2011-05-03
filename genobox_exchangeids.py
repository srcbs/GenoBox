#!/panvol1/simon/bin/python2.7

import argparse
import sys

def exchange_id(handle, outhandle, x, b, d, f, s):
   ''' Exchange id of input line, column x '''
   for line in handle:
      if line.startswith('#'):
         outhandle.write(line)
         continue
      
      # filter on f (if set)
      if f:
         if not line.startswith(f):
            outhandle.write(line)
            continue
      
      # filter on s (if set)
      if s:
         if line.startswith(s):
            outhandle.write(line)
            continue
      
      # do the exchange
      line = line.rstrip()
      if line == '':
         continue
      
      fields = line.split('\t')
      fields[x] = d[fields[x]]
      
      newline = ('\t').join(fields) + '\n'
      outhandle.write(newline)
   


# create the parser
parser = argparse.ArgumentParser(description='''
   Exchange ids in file -a column X (0-based) based on mapping provided in file b.
   File b must have common ids in column 0 and new id in column 1
   ''')

# add the arguments
parser.add_argument('--a', help='file to exchange ids [stdin]', default=None)
parser.add_argument('--x', help='column to swap ids in [0]', default=0, type=int)
parser.add_argument('--b', help='file containing annotation')
parser.add_argument('--f', help='only if line starts with []', default=None)
parser.add_argument('--s', help='skip line if starts with []', default=None)
parser.add_argument('--o', help='output file [stdout]', default=None)
parser.add_argument('--log', help='log level [INFO]', default='info')

# parse the command line
args = parser.parse_args()
#args = parser.parse_args('--a rmsk_build37.sort.genome --x 0 --b chr2gi.build37 --o rmsk_build37.gi.genome'.split())
#args = parser.parse_args(' --a hs_ref_GRCh37_all_gatk_tmp.fa --x 0 --b gi_fa2number.build37_rCRS --o hs_ref_GRCh37_all_gatk_number.fa --f ">"'.split())
#args = parser.parse_args(' --x 0 --b gi_fa2number.build37_rCRS --f ">"'.split())

# read annotation to dict
d = {}
bhandle = open(args.b, 'r')
for line in bhandle:
   line = line.rstrip()
   fields = line.split('\t')
   d[fields[0]] = fields[1]

bhandle.close()

# set outhandle
if args.o:
   outhandle = open(args.o, 'w')
else:
   outhandle = sys.stdout

# set inhandle
if args.a:
   inhandle = open(args.a, 'r')
else:
   inhandle = sys.stdin

# parse and exchange ids
exchange_id(inhandle, outhandle, args.x, args.b, d, args.f, args.s)

outhandle.close()
inhandle.close()
