#!/panvol1/simon/bin/python2.7

from genobox_modules import rm_files
import subprocess
import os

rm_files(['run_genobox_velveth.*', 'run_genobox_velvetg.*', 'run_genobox_interleave.*', '*.interleaved', 'pbsjob.tmp*', 'run_genobox_velvetaccept.*', 'run_mlst_trim.*'])

if os.path.exists('trimmed'):
   subprocess.call('rm -r trimmed/', shell=True)
