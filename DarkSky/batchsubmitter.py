from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)


flines = open('submittemplatehighmem.sh').readlines()


#speclist = ['05', '06', '13', '25', '39', '43', '46' ,'48']
#flist = sorted(glob('WISE/wise-allsky-cat-part??'))
flist = sorted(glob('AllWISE/wise-allwise-cat-part??'))
#flist = ['AllWISE/wise-allwise-cat-part46', 'AllWISE/wise-allwise-cat-part47', 'AllWISE/wise-allwise-cat-part48']
jname = "DarkSkyQuery"

njobs = int(sys.argv[1])

h = 0.688062

edgeclearance = 2100.

totallength = 8000./h




for i in range(0, njobs):
    CX = np.random.uniform(2100., totallength-2100.)
    CY = np.random.uniform(2100., totallength-2100.)
    CZ = np.random.uniform(2100., totallength-2100.)
    jobname = jname+str(i)+str(CX)+str(CY)+str(CZ)
    fout = open(jobname+'.slurm', "w")
    jobline = 'python AnalyzeBox.py '+' -x '+str(CX)+' -y '+str(CY)+' -z '+str(CZ)

    for line in flines:
        fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
    fout.close()
    os.system('chmod +x ' + jobname+'.slurm')
    os.system('sbatch -p icecube '+ jobname+'.slurm')
    time.sleep(0.01)
                        #raw_input("Press Enter to Continue")
