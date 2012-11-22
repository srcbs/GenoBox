#!/panvol1/simon/bin/python2.7

from __future__ import division

import argparse
import logging
import pysam
from collections import defaultdict

def average(values):
    '''Computes the arithmetic mean of a list of numbers.'''
    return sum(values, 0.0) / len(values)


parser = argparse.ArgumentParser(description='''Calculate average depth from bam file.''')

parser.add_argument('i', help='input bamfile')
parser.add_argument('--log', help='log level [info]', default='info')

args = parser.parse_args()
#args = parser.parse_args('libA_hg19.sort.q30.bam.rmdup.bam'.split())
#args = parser.parse_args('alignment/kleb_10_213361.flt.sort.rmdup.bam'.split())
#args = parser.parse_args('A.flt.sort.bam'.split())


# set logging
logger = logging.getLogger('genobox.py')
hdlr = logging.FileHandler('genobox.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# read in from bam
samfile = pysam.Samfile(args.i, 'rb')

# get length of references from header
refl = defaultdict(lambda:0) 
header = samfile.header['SQ']
lengths = {}
for ref in header:
   refl[ref['SN']] = ref['LN']
   lengths[ref['SN']] = []

# parse bam
total_count = 0
not_aligned = 0

d_counts = defaultdict(lambda:0)
d_sums = defaultdict(lambda:0)
s = 0
with samfile as bam_file:
   for alignedRead in bam_file:
      total_count += 1
      if alignedRead.flag == 4L:
         not_aligned += 1
         continue
      refname = samfile.getrname( alignedRead.rname )
      d_counts[refname] += 1
      d_sums[refname] += alignedRead.rlen
      s += alignedRead.rlen


print 'Counted %i, with %i not aligned' % (total_count, not_aligned)
logger.info('Counted %i, with %i not aligned' % (total_count, not_aligned))

print 'Average length of mapped reads %g' % (s/total_count)
logger.info('Average length of mapped reads %g' % (s/total_count))

# output
print 'Average depth of %s:' % args.i
print 'Chrom\tAvgReadDepth\tSumReadLength\tChromLength\tAverageReadLength'
logger.info('Average depth of %s:' % args.i)
logger.info('Chrom\tAvgReadDepth\tSumReadLength\tChromLength\tAverageReadLength')
totalLength = 0
totalReadlength = 0
for key in d_counts.keys():
   totalLength += refl[key]
   totalReadlength += d_sums[key]
   avg = d_sums[key] / refl[key]
   avgReadL = d_sums[key]/d_counts[key]
   print '%s: %.1f\t%i\t%i\t%1.f' % (key, avg, d_sums[key], refl[key], avgReadL)
   logger.info('%s: %.1f\t%i\t%i\t%1.f' % (key, avg, d_sums[key], refl[key], avgReadL))

print 'Total: %.1f\t%i\t%i' % (totalReadlength / totalLength, totalReadlength, totalLength)
logger.info('Total: %.1f\t%i\t%i' % (totalReadlength / totalLength, totalReadlength, totalLength))


