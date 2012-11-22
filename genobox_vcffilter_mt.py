#!/panvol1/simon/bin/python2.7

import argparse
import sys

def read_rmsk(rmsk):
   '''Read rmsk file and return dictionary of positions to filter
   only keep positions that matches chromosome being analyzed'''
   
   # because this is only used for MT dna only keep chromosomal positions from rmsk file
   # that matches MT (to minimize memory)
   fh = open(rmsk, 'r')
   dict = {}
   for line in fh:
      line = line.rstrip()
      fields = line.split('\t')
      if fields[0] == 'MT':
         for pos in xrange(int(fields[1]), int(fields[2])+1, 1):
            dict[fields[0], str(pos)] = 1
   
   return dict

def manual_rmsk_filter(args):
   '''Manually filter vcf using rmsk file'''
   
   # read in rmsk
   rmsk_dict = read_rmsk(args.rmsk)
   
   # read in vcf
   if args.vcf == '-':
      fh = sys.stdin
   else:
      fh = open(args.vcf, 'r')
   
   # filter vcf
   for line in fh:
      if line.startswith('#'):
         sys.stdout.write(line)
         continue
      
      fields = line.split('\t')
      if rmsk_dict.has_key((fields[0], fields[1])):
         pass
      else:
         sys.stdout.write(line)
   
   

if __name__ == '__main__':
   
   # create the parser
   parser = argparse.ArgumentParser(description=
   '''Filter rmsk for MT dna (because it fails for bedtools)''')
   
   # add the arguments
   parser.add_argument('--vcf', help='input vcf [-]', default='-')
   parser.add_argument('--rmsk', help='input rmsk', required=True)
   parser.add_argument('--log', help='log level [info]', default='info')
   
   args = parser.parse_args()
   #args = parser.parse_args('--vcf tmp/tmp.BI16.flt.sort.rmdup.realign.MT.raw.vcf.gz.ref.vcf --rmsk /panvol1/simon/projects/arctic/references/build37/rmsk.build37'.split())
   
   manual_rmsk_filter(args)
