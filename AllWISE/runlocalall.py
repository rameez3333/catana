from glob import glob
import os, sys

flist = glob('wise-allwise-cat-part??')

for fname in flist:
    os.system('python processor_AllWISE.py 16 '+fname)
