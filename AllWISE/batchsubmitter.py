from glob import glob
import os, sys
import numpy as np

nside = sys.argv[1]


flines = open('submittemplate.sh').readlines()


flist = sorted(glob('wise-allwise-cat-part??'))

fdonelist = sorted(glob('*result.txt'))

for fname in flist:
    for glcut in np.arange(10,11):
        for jcut in ['0']:
            jobname = str(glcut)+str(jcut)+'allwise'+fname[-2:]+'job'
            fout = open(jobname+str(glcut)+'.slurm', "w")
            jobline = 'python processor_AllWISE.py '+nside+ ' '+ fname + ' '+ str(glcut)+' '+str(jcut)
            for line in flines:
                fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
            fout.close()
            os.system('chmod +x ' + jobname+str(glcut)+'.slurm')
            os.system('sbatch -p icecube '+ jobname+str(glcut)+'.slurm')
