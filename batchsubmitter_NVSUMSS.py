from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)


flines = open('submittemplate.sh').readlines()


#speclist = ['05', '06', '13', '25', '39', '43', '46' ,'48']
#flist = sorted(glob('WISE/wise-allsky-cat-part??'))
flist = sorted(glob('AllWISE/wise-allwise-cat-part??'))
#flist = ['AllWISE/wise-allwise-cat-part46', 'AllWISE/wise-allwise-cat-part47', 'AllWISE/wise-allwise-cat-part48']
jname = "NVSUMSS"


for gcut in [10, 15, 20]:
    for sgcut in [0]:
        for tmrscut in [0.00028, 0.0028, 0.01]:
            for tcut in [10, 15, 20, 25, 30, 35, 40, 45, 50]:
                for patch in [-30, -40]:
                    #for jcut in [1, 0]:
                        #for fgo in flist:
                            jobname = jname+str(gcut)+str(sgcut)+str(tmrscut)+str(tcut)+str(patch)+'job'
                            fout = open(jobname+'.slurm', "w")
                            jobline = 'python processor_NVSUMSSminus2MRS.py -f NVSS/'+' -g '+str(gcut)+' -s '+str(sgcut)+' -r 32 -d -c -m '+str(tmrscut)+' -t '+str(tcut)+' -p '+str(patch)
                            #if diffcuts:
                                #jobline = jobline + ' -d '
                            #if jcut:
                                #jobline = jobline + ' -j '
                            for line in flines:
                                fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
                            fout.close()
                            os.system('chmod +x ' + jobname+'.slurm')
                            os.system('sbatch -p icecube '+ jobname+'.slurm')
                            #time.sleep(300)
                            #raw_input("Press Enter to Continue")
