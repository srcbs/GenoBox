#!/panvol1/simon/bin/python

# Creator: Anna Sapfo
# Modified: Simon Rasmussen

import argparse
import subprocess
import pysam
import sys
import os 
import random
import re
try:
    import pylab as py
    import numpy as np
    figure=True
except:
    figure=False

def set_abspath():
   '''Returns absolute path of file to argparse'''
   class SetAbspath(argparse.Action):
      def __call__(self, parser, args, filenames, option_string=None):
         import os
         if type(filenames) == str:
            f_abs = os.path.abspath(filenames)
            setattr(args, self.dest, f_abs)
         elif type(filenames) == list:
            new_list = []
            for f in filenames:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, filenames)
   return SetAbspath

def mapping(se, pe1, pe2, dico_c):
   '''bwa mapping of fastq files'''
   
   # mapping se files
   samfiles = []
   
   if args.se:
      for s in args.se:
         out = os.path.split(s)[1]+'.bam'
         samfiles.append(out)
         call = "/panvol1/simon/bin/bwa-0.6.2/bwa aln -t %i -l 10000 %s %s | /panvol1/simon/bin/bwa-0.6.2/bwa samse %s -  %s | samtools view -Sb - > %s" % (args.n, args.fa, s, args.fa, s, out)
         ec = subprocess.call(call, shell=True)
   
   # map pe files
   if args.pe1:
      for i in range(len(args.pe1)):
         p1 = args.pe1[i]
         p2 = args.pe2[i]
         sai1 = os.path.split(p1)[1]+'.sai'
         sai2 = os.path.split(p2)[1]+'.sai'
         out = os.path.split(p1)[1]+'.bam'
         samfiles.append(out)
         aln1_call = "/panvol1/simon/bin/bwa-0.6.2/bwa aln -t %i -l 10000 %s %s > %s" % (args.n, args.fa, p1, sai1)
         aln2_call = "/panvol1/simon/bin/bwa-0.6.2/bwa aln -t %i -l 10000 %s %s > %s" % (args.n, args.fa, p2, sai2)
         sampe_call = "/panvol1/simon/bin/bwa-0.6.2/bwa sampe %s %s %s %s %s | samtools view -Sb - > %s" % (args.fa, sai1, sai2, p1, p2, out)
         ec = subprocess.call(aln1_call, shell=True)
         ec = subprocess.call(aln2_call, shell=True)
         ec = subprocess.call(sampe_call, shell=True)
         
         # remove saifiles
         call = 'rm %s %s' % (sai1, sai2)
         ec = subprocess.call(call, shell=True)
   
   # merge sams to one file
   if len(samfiles) > 1:
      call = 'samtools merge -f %s %s' % (dico_c['SAMfile_all'], ' '.join(samfiles))
      ec = subprocess.call(call, shell=True)
   else:
      if samfiles[0] == dico_c['SAMfile_all']:
         pass
      else:
         call = 'mv %s %s' % (samfiles[0], dico_c['SAMfile_all'])
         ec = subprocess.call(call, shell=True)
   
   # remove samfiles
   for s in samfiles:
      call = 'rm %s' % s
      ec = subprocess.call(call, shell=True)


def input_bam(bams, dico_c):
   '''Copy bam to SAMfiles_all or merge to Samfiles_all'''
   
   if len(bams) > 1:
      call = 'samtools merge -f %s %s' % (dico_c['SAMfile_all'], ' '.join(bams))
      ec = subprocess.call(call, shell=True)
   else:
      call = 'cp %s %s' % (bams[0], dico_c['SAMfile_all'])
      ec = subprocess.call(call, shell=True)


def split_bam(bam, nlines_block, prefix, subset):
   '''Split bam into nlines_pr_block bams'''
      
   samfile = pysam.Samfile(bam, 'rb')
   k = 'a'
   c = 0
   fh_out = pysam.Samfile('%s.%s.bam' % (prefix, k), 'wb', template=samfile)
   kmers = {'ACG':1, 'CTA':1, 'GAA':1, 'TGT': 1}
   
   for read in samfile:
      c += 1
      if subset:
         if kmers.has_key(read.seq[0:3]):
            pass
         else:
            continue
      
      if c > nlines_block:
         c = 0
         k = chr(ord(k)+1)
         fh_out.close()
         fh_out = pysam.Samfile('%s.%s.bam' % (prefix, k), 'wb', template=samfile)
         d = fh_out.write(read)
      else:
         d = fh_out.write(read)

def sort_bam(bams):
   '''Sorting input bams (list)'''
   
   sorted = []
   for b in bams:
      s = os.path.split(b)[1].replace('.bam', '.sort')
      sorted.append(s+'.bam')
      call = 'samtools sort %s %s' % (b, s)
      ec = subprocess.call(call, shell=True)
   
   return sorted


if __name__ == '__main__':
   
   # create the parser
   parser = argparse.ArgumentParser(prog='bam_saturation.py', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=130), usage='%(prog)s [options]', description='''Calculate saturation plots from data. Input either fastq or bams''')
   
   parser.add_argument('--se', help='input SE fastq', nargs='+', action=set_abspath())
   parser.add_argument('--pe1', help='input PE fastq pair 1', nargs='+', action=set_abspath())
   parser.add_argument('--pe2', help='input PE fastq pair 2', nargs='+', action=set_abspath())
   parser.add_argument('--bams', help='input bamfiles', nargs='+', action=set_abspath())
   parser.add_argument('--fa', help='reference genome', action=set_abspath())
   parser.add_argument('--subsample', help='subsample reads [False]', default=False, action='store_true')
   parser.add_argument('--sample', help='name of sample/output dir', required=True)
   parser.add_argument('--n', help='number of cores to use [1]', default=1, type=int)
   parser.add_argument('--blocks', help='number of blocks to use [10]', default=10, type=int)
   parser.add_argument('--replicates', help='number of replicates [1]', default=1, type=int)
   
   
   args = parser.parse_args()
   #args = parser.parse_args('--se test.se.fastq.gz --pe1 test.pe1.fastq.gz --pe2 test.pe2.fastq.gz --sample Aleutian_1_2_2.01082012 --fa /panvol1/simon/projects/arctic/references/build37/hs.build37.1.fa --n 8 --blocks 10 --replicates 1'.split())


dico_c={}

# maintaining Annas dict
dico_c['dirSaturation']=args.sample
dico_c['basefile']=args.sample
dico_c['hg']=args.fa
dico_c['ncpus'] = args.n
nblocks=args.blocks
ncpus=args.n
nrepeats=args.replicates

dico_c['SAMfile_all']="%(basefile)s_all.bam"%dico_c
dico_c['HEADERfile']="Header_%(basefile)s.txt"%dico_c
dico_c['index']=dico_c['hg']

# change dir
if not os.path.exists(dico_c['dirSaturation']):
   os.makedirs(dico_c['dirSaturation'])

os.chdir(dico_c['dirSaturation'])

# map or use bams
if args.se or args.pe1:
   mapping(args.se, args.pe1, args.pe2, dico_c)
elif args.bams:
   input_bam(args.bams, dico_c)
else:
   raise ValueError('Input must be given as either fastq or bam')


print "#3. get header "
commandSEP1="samtools view -H %(SAMfile_all)s > %(HEADERfile)s;"%dico_c
print commandSEP1
os.system(commandSEP1)
headerC=int(os.popen("wc -l %(HEADERfile)s"%dico_c).readlines()[0].split()[0])
#headerC=26
#print headerC
dico_c['headerC']=headerC
print "nlines in headerC: ",headerC

print "#4. count the number of lines in sam file"
commandCount="samtools view -c %(SAMfile_all)s"%dico_c
print commandCount
nlines=int(os.popen(commandCount).readlines()[0].rstrip())
#nlines=3992689
print "nlines ",nlines
nreads=int(1.*nlines/nblocks)+1
dico_c['nlines_block']=nreads
print "nlines_block: ",dico_c['nlines_block']


print "#6. split the sam file"
split_bam(dico_c['SAMfile_all'], dico_c['nlines_block'], dico_c['basefile'], args.subsample)

# create list of split bam-files
DirList=os.listdir('.')
ListSamFiles = []
pattern = re.compile(dico_c['basefile']+'\.\w+\.bam')
for f in DirList:
   if pattern.match(f): ListSamFiles.append(f)

# sort bams
sortedBamFiles = sort_bam(ListSamFiles)

# mv sorted_bams to org bams
for i in range(len(ListSamFiles)):
   cmd = 'mv %s %s' % (sortedBamFiles[i], ListSamFiles[i])
   ec = subprocess.call(cmd, shell=True)


print "#9. every repeat: shuffle files differently: needto keep each part"
STRRm=''
STRMap=''

for repeat in range(nrepeats):
    print "REPEAT: ",repeat
    random.shuffle(ListSamFiles)
    print ListSamFiles
    
    ColleFiles=[ListSamFiles[:xx+1] for xx in range(len(ListSamFiles))]
    print len(ColleFiles)
    #ColleFiles[xx+1 for xx in range(len(ListSamFiles))]
    print ColleFiles
    
    for colI,col in enumerate(ColleFiles):
        if colI==0:
            commandMerge='cp %s %s.%s.bam'%(col[0],dico_c['basefile'],colI)
        else:
            commandMerge="samtools merge -f %(basefile)s.%(count)s.bam %(col)s"
        
        dico_c['count']=colI
        dico_c['col']='\t'.join(col)
        
        print commandMerge%dico_c
        os.system(commandMerge%dico_c)
        
        
    print "#10. rmdup each replicate "
    
    for colI in range(len(ColleFiles)):
        dico_c['count']=colI
        
        commandRMDUP="java -XX:ParallelGCThreads=8 -XX:+UseParallelGC -XX:-UsePerfData -Xms4500m -Xmx4500m -jar /panvol1/simon/bin/picard-tools-1.56/MarkDuplicates.jar INPUT=%(basefile)s.%(count)s.bam OUTPUT=%(basefile)s.%(count)s.srt.rmdup.bam METRICS_FILE=/dev/null REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=/panvol1/simon/tmp VALIDATION_STRINGENCY=SILENT"%dico_c
        
        print commandRMDUP
        os.system(commandRMDUP)
        
        
    print "#11. count each for each case 1. the number or reads 2. the number of reads that map to the human "
    
    X=[];YRm=[];YMap=[];
    
    for colI in range(len(ColleFiles)):
        dico_c['count']=colI
        
        commandCountRm="samtools view -c %(basefile)s.%(count)s.srt.rmdup.bam"%dico_c
        commandCountMap="samtools view -F4 -c %(basefile)s.%(count)s.srt.rmdup.bam"%dico_c
        print commandCount
        
        nlinesRmdup=int(os.popen(commandCountRm).readlines()[0].rstrip())
        nlinesMap=int(os.popen(commandCountMap).readlines()[0].rstrip())
        
        print nlinesRmdup
        print nlinesMap
        
        X.append(colI)
        
        YRm.append(nlinesRmdup)
        YMap.append(nlinesMap)
        
        STRRm+="%i\t%i\t%i\n"%(repeat,colI,nlinesMap)
        STRMap+="%i\t%i\t%i\n"%(repeat,colI,nlinesRmdup)
        
    
    print STRRm
    print STRMap
    
    print "#12. remove the collection of samfiles: rmdup "
    
    commandRM="rm %(basefile)s.*.bam  "%dico_c
    print commandRM
    os.system(commandRM)
    
    if figure:
        py.subplot(211)
        py.title("Rmdup, clonality, all trimmed reads")
        py.plot(nreads*(1+np.array(X)),YRm,'o-',label='observed #%s'%repeat)
        
        py.subplot(212)
        py.title("Rmdup, clonality, all trimmed reads that map to human")
        py.plot(nreads*(1+np.array(X)),YMap,'o-',label='observed #%s'%repeat)


if figure:
    py.subplot(211)
    py.plot(nreads*(1+np.array(X)),nreads*(1+np.array(X)),'ro-',label='expected')
    py.subplot(212)    
    py.plot(nreads*(1+np.array(X)),YMap[0]*(1+np.array(X)),'ro-',label='expected')
    py.legend(loc='best')
    py.savefig("%(basefile)s.pdf"%dico_c)

outfile=open("%(basefile)s_Rm.txt"%dico_c,'w')
outfile.write(STRRm)
outfile.close()
outfile=open("%(basefile)s_Map.txt"%dico_c,'w')
outfile.write(STRMap)
outfile.close()


#print "#13. remove each part "
commandRM="rm %(basefile)s.?.bam"%dico_c
print commandRM
os.system(commandRM)

# remove original bam
commandRM="rm %s*" % dico_c['SAMfile_all']
print commandRM
os.system(commandRM)



# Map=read.table("../Aleutian_1_2_2.01082012_Map.txt", header=FALSE)
# Rm=read.table("../Aleutian_1_2_2.01082012_Rm.txt", header=FALSE)
#
# Map_s=read.table("Aleutian_1_2_2.01082012_Map.txt", header=FALSE)
# Rm_s=read.table("Aleutian_1_2_2.01082012_Rm.txt", header=FALSE)
#
# par(mfrow=c(1,2))
# plot(Map[,3], Rm[,3], pch=16, xlab="All mapped reads", ylab="Mapped rmduped reads", main="All reads")
# abline(0,1, col="grey")
# plot(Map_s[,3], Rm_s[,3], pch=16, xlab="All mapped reads", ylab="Mapped rmduped reads", main="Read subset")
# abline(0,1, col="grey")
# dev.print("test.saturation.full_vs_subset.pdf", device=pdf)
#