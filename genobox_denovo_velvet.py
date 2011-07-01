#!/panvol1/simon/bin/python

import argparse

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


def which(L, st):
    out = []
    for i,item in enumerate(L):
        if item == st:
            out.append(i)
    return out

def interleave(i):
    '''Identify paired end read files that should be interleaved'''
    
    def merge(reads, format, interleaved):
        '''Perform merging'''
        
        import genobox_modules
        paths = genobox_modules.setSystem()
        
        # change to genobox_denovo_interleave.py <file1> <file2> <format> <out.gz>
        cmd = '%sgenobox_denovo_interleave.py %s %s %s %s' % (paths['genobox_home'], reads[0], reads[1], format, interleaved)
        return cmd
    
    # check if format can be merged (fasta, fastq)
    if i:
        if i[0].find('fasta') > -1 or i[0].find('fastq') > -1:
            # check if merging is needed
            if len(i) > 2:
                format = i[0]
                reads = i[1:3]
                interleaved = reads[0] + '.interleaved.gz'
                # set new format
                if format == 'fasta': new_format = 'fasta.gz'
                elif format == 'fastq': new_format = 'fastq.gz'
                else: new_format = format
                call = merge(reads, format, interleaved)
                return (call, [format, interleaved])
            else:
                return (None, i)
        else:
            return (None, i)
    else:
        return (None, i)
        

def velveth_calls(args):
    '''Return velveth calls'''
    
    import genobox_modules
    
    paths = genobox_modules.setSystem()
    cmd = '%svelveth' % paths['velvet_home']
    
    velveth_calls = []
    if len(args.ksizes) == 1:
        arg = ' %s %s ' % (args.outpath, args.ksizes[0])
        if args.short: arg = arg + ' '.join(args.short) + ' '
        if args.short2: arg = arg + ' '.join(args.short2) + ' '
        if args.shortPaired: arg = arg + ' '.join(args.shortPaired) + ' '
        if args.shortPaired2: arg = arg + ' '.join(args.shortPaired2) + ' '
        if args.long: arg = arg + ' '.join(args.long) + ' '
        if args.longPaired: arg = arg + ' '.join(args.longPaired) + ' '
        call = cmd + arg
        velveth_calls.append(call)
    
    # repeat above for range(args.ksizes[0], args.ksizes[1], 2)


def start_assembly(args):
    '''Start assembly'''
    
    # checking if files needs to be interleaved
    interleave_dict = {}    
    interleave_dict['short'] = interleave(args.short)[0] ; args.short = interleave(args.short)[1]
    interleave_dict['short2'] = interleave(args.short2)[0] ; args.short2 = interleave(args.short2)[1]
    interleave_dict['shortPaired'] = interleave(args.shortPaired)[0] ; args.shortPaired = interleave(args.shortPaired)[1]
    interleave_dict['shortPaired2'] = interleave(args.shortPaired2)[0] ; args.shortPaired2 = interleave(args.shortPaired2)[1]
    interleave_dict['long'] = interleave(args.long)[0] ; args.long = interleave(args.long)[1]
    interleave_dict['longPaired'] = interleave(args.longPaired)[0] ; args.longPaired = interleave(args.longPaired)[1]
    
    # creating calls
    interleave_calls = []
    for key,value in interleave_dict.items():
        if value:
            interleave_calls.append(value)
    
    # create velveth calls
    



# create the parser
parser = argparse.ArgumentParser(description='''Run Velvet denovo assembly\nRead input is given by eg. --short <format> <reads>\n
    Format can be: fasta, fastq, raw, fasta.gz, fastq.gz, raw.gz, sam, bam\n''')

# add the arguments
parser.add_argument('--short', help='input read format and short reads', nargs='+', action=required_nargs(0,3))
parser.add_argument('--shortPaired', help='input read format and short paired reads', nargs='+', action=required_nargs(0,3))
parser.add_argument('--short2', help='input read format and short2 reads', nargs='+', action=required_nargs(0,3))
parser.add_argument('--shortPaired2', help='input read format and short paired2 reads', nargs='+', action=required_nargs(0,3))
parser.add_argument('--long', help='input read format and long reads', nargs='+', action=required_nargs(0,3))
parser.add_argument('--longPaired', help='input read format and long paired reads', nargs='+', action=required_nargs(0,3))
parser.add_argument('--ksizes', help='kmers to run assemblies for (single no or range) [33]', nargs='+', default=[33])
parser.add_argument('--outpath', help='name of run, also output dir [velvet_assembly]', default='velvet_assembly')
parser.add_argument('--min_length', help='mininum length to report contig [100]', default=100, type=int)
parser.add_argument('--cov_cut', help='removal of low coverage nodes after Tour bus (float) [None]', default=None, type=float)
parser.add_argument('--exp_cov', help='Expected mean coverage [auto]', default='auto')
parser.add_argument('--ins_length', help='insert size (reads included) [None]', default=None)
parser.add_argument('--add_velveth', help='additional parameters to velveth', default=None)
parser.add_argument('--add_velvetg', help='additional parameters to velvetg', default=None)
parser.add_argument('--n', help='number of threads for parallel run [4]', default=4, type=int)
parser.add_argument('--queue', help='queue to submit jobs to (idle, cbs, cpr, cge, urgent) [cbs]', default='cbs')
parser.add_argument('--log', help='log level [info]', default='info')

args = parser.parse_args()
#args = parser.parse_args('--short fastq test_1.fq test_2.fq --ksizes 33 --outpath test'.split())
#args = parser.parse_args('--short fastq.gz interleaved.fastq.gz --ksizes 33 --outpath test'.split())

start_assembly(args)




