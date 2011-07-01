#!/panvol1/simon/bin/python2.7

import argparse
import sys

def read_bed(bed):
    '''Read rmsk file and return dictionary of positions to filter
        only keep positions that matches chromosome being analyzed'''
    
    fh = open(bed, 'r')
    dict = {}
    for line in fh:
        line = line.rstrip()
        fields = line.split('\t')
        for pos in xrange(int(fields[1])+1, int(fields[2])+1, 1):
            dict[fields[0], str(pos)] = 1
    
    return dict

def manual_bed_filter(vcf, bed, vcf_out):
    '''Manually filter vcfgz using bed file'''
    
    # get rmsk as dict
    rmsk_dict = read_bed(bed)
    
    # open in/outhandles
    if vcf == 'stdin':
        fh = sys.stdin
    else:
        fh = open(vcf, 'r')
    
    if vcf_out == 'stdout':
       fh_out = sys.stdout
    else:
        fh_out = open(vcf_out, 'w')
    
    #parse file
    for line in fh:
        if line.startswith('#'):
            fh_out.write(line)
            continue
        
        fields = line.split('\t')
        if rmsk_dict.has_key((fields[0], fields[1])):
            fh_out.write(line)
        else:
            pass


# create the parser
parser = argparse.ArgumentParser(description='''
    Filter vcf based on input bed
    ''')

# add the arguments
parser.add_argument('--vcf', help='input vcf [stdin]', default='stdin')
parser.add_argument('--bed', help='input bed', required=True)
parser.add_argument('--o', help='output vcf [stdout]', default='stdout')

# parse the command line
args = parser.parse_args()

manual_bed_filter(args.vcf, args.bed, args.o)
