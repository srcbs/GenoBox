#!/panvol1/simon/bin/python2.7

from __future__ import division

import argparse
import sys

def calc_TiTv_ratio(i):
    '''Calculate transition / transversion ratio of vcf'''
    
    if args.i == 'stdin':
        fh = sys.stdin
    else:
        fh = open(i, 'r')
    
    # parse
    transversion = 0
    transition = 0
    total = 0
    for line in fh:
        if line.startswith('#'):
            continue
        
        total += 1
        fields = line.rstrip().split('\t')
        ref = fields[3]
        alt = fields[4]
        genotype = ref+alt
        
        # purines: A, G
        # pyrimidines: C, T, U
        
        # Transition: C->T, T->C, A->G, G->A
        # Transversion: C->A, T->A, T->G, C->G, A->C, A->T, G->T, C->A
        if "A" in genotype and "T" in genotype: transversion += 1
        if "C" in genotype and "A" in genotype: transversion += 1 
        if "T" in genotype and "G" in genotype: transversion += 1
        if "C" in genotype and "G" in genotype: transversion += 1 
        if "C" in genotype and "T" in genotype: transition += 1 
        if "A" in genotype and "G" in genotype: transition += 1 
    
    titv = transition / transversion
    print "Transitions: %i" % transition
    print "Transversions: %i" % transversion
    print "Total: %i" % total
    print "Ratio: %.2f" % titv

# create the parser
parser = argparse.ArgumentParser(description='''
   Calculate transition and transversion ratio from input vcf or stdin
   ''')

# add the arguments
parser.add_argument('--i', help='input vcf [stdin]', default='stdin')

# parse the command line
args = parser.parse_args()
#args = parser.parse_args('--i tmp.vcf'.split())

if __name__ == '__main__':
    calc_TiTv_ratio(args.i)
