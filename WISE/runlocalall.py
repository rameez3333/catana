from glob import glob
import os, sys

flist = sorted(glob('wise-allsky-cat-part??'))

for fname in flist:
    os.system('python processor_AllSky.py 16 '+fname)
