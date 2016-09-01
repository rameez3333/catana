from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)


flines = open('submittemplate.sh').readlines()


##speclist = ['05', '06', '13', '25', '39', '43', '46' ,'48']
##flist = sorted(glob('WISE/wise-allsky-cat-part??'))
#flist = sorted(glob('AllWISE/wise-allwise-cat-part??'))
##flist = ['AllWISE/wise-allwise-cat-part46', 'AllWISE/wise-allwise-cat-part47', 'AllWISE/wise-allwise-cat-part48']
#jname = "AllWISE"


#for gcut in [10, 15, 20]:
    #for sgcut in [0,5,10]:
        #for tmrscut in [0.00028, 0.0028, 0.01]:
            #for tmzcut in [0.02, 0.05, 0.09]:
                #for diffcuts in [1, 0]:
                    #for jcut in [1, 0]:
                        #for fgo in flist:
                            #jobname = jname+fgo[29:]+str(gcut)+str(sgcut)+str(tmzcut)+str(tmrscut)+str(diffcuts)+str(jcut)+'job'
                            #fout = open(jobname+'.slurm', "w")
                            #jobline = 'python processor_WISE.py -f AllWISE/'+' -i '+fgo+' -g '+str(gcut)+' -s '+str(sgcut)+' -r 32 -m '+str(tmrscut)+' -z '+str(tmzcut)
                            #if diffcuts:
                                #jobline = jobline + ' -d '
                            #if jcut:
                                #jobline = jobline + ' -j '
                            #for line in flines:
                                #fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
                            #fout.close()
                            #os.system('chmod +x ' + jobname+'.slurm')
                            #os.system('sbatch -p icecube '+ jobname+'.slurm')
                        #time.sleep(300)
                            #raw_input("Press Enter to Continue")

jname = 'Simjob_'

for ncat in [600000]:
    for vel in [369.]:# 738., 1107., 1476., 1845.]:
        for i in range(200,1000):
                            jobname = jname+str(ncat)+str(vel)+str(i)+'job'
                            fout = open(jobname+'.slurm', "w")
                            jobline = 'python DoSim.py -r 168.0 -d -7.0'+' -v '+str(float(vel))+' -n '+str(ncat)+' -s '+str(i) + ' -p 0'
                            for line in flines:
                                fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
                            fout.close()
                            os.system('chmod +x ' + jobname+'.slurm')
                            os.system('sbatch -p icecube '+ jobname+'.slurm')
                            raw_input('test')
