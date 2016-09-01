from glob import glob
import os, sys
import numpy as np

nside = sys.argv[1]


flines = open('submittemplate.sh').readlines()


#speclist = ['05', '06', '13', '25', '39', '43', '46' ,'48']

flist = sorted(glob('wise-allsky-cat-part??'))


for fname in flist:
    for glcut in np.arange(10,11):
        for jcut in ['0']:
    #fname = 'wise-allsky-cat-part'+spec
            jobname = str(glcut)+str(jcut)+'allsky'+fname[-2:]+'job'
            fout = open(jobname+str(glcut)+'.slurm', "w")
            jobline = 'python processor_AllSky.py '+nside+ ' '+ fname + ' '+ str(glcut)+' '+str(jcut)
            for line in flines:
                fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
            fout.close()
            os.system('chmod +x ' + jobname+str(glcut)+'.slurm')
            os.system('sbatch -p icecube '+ jobname+str(glcut)+'.slurm')
