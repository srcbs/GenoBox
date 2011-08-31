#!/panvol1/simon/bin/python2.7

import argparse

def parse_accept(outpath):
   '''Take best assembly and remove worst. Parses velvet_parse.txt'''
   
   import re
   import subprocess
   
   # parse velvet_parse.txt file
   reg = re.compile('(\d+)')
   fh = open('velvet_parse.txt', 'r')
   for line in fh:
      line = line.rstrip()
      if line.find('Best assembly') > -1:
         match = reg.findall(line)
         best_assembly = match[0]
      if line.find('Rank') > -1:
         match = reg.findall(line)
         rm_assemblies = match
         rm_assemblies.remove(best_assembly)
         
   calls = []
   # move that assembly to final args.outpath
   cmd = 'mv'
   arg = ' %s_%s %s' % (outpath, best_assembly, outpath)
   calls.append(cmd+arg)
   
   # remove other assemblies
   cmd = 'rm'
   for k in rm_assemblies:
      arg = ' -r %s_%s' % (outpath, k)
      calls.append(cmd+arg)
   
   # run processes
   for call in calls:
      subprocess.call(call, shell=True)


if __name__ == '__main__':
   parser = argparse.ArgumentParser(prog='genobox_denovo_accept.py', description='''Keep best assembly and remove others''')
   
   parser.add_argument('outpath', help='outpath given to assemblies')
   
   args = parser.parse_args()
   
   # start 
   parse_accept(args.outpath)