#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description='''
   Create library file from tab sep file of "filename\tlibname"
   ''')

# add the arguments
parser.add_argument('--i', help='input file')
parser.add_argument('--o', help='output file')
parser.add_argument('--log', help='log level [INFO]', default='info')

args = parser.parse_args()

lib = {}
fh = open(args.i, 'r')
outhandle = open(args.o, 'w')

for line in fh:
   line = line.rstrip()
   fields = line.split('\t')
   if lib.has_key(fields[1]):
      lib[fields[1]].append(fields[0])
   else:
      lib[fields[1]] = [fields[0]]

for key in lib.keys():
   outhandle.write('lib%s\t%s\n' % (key, ' '.join(lib[key])))
