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
jname = "AllWISE"

#flist = [ 'AllWISE/wise-allwise-cat-part03', 'AllWISE/wise-allwise-cat-part05', 'AllWISE/wise-allwise-cat-part10', 'AllWISE/wise-allwise-cat-part18', 'AllWISE/wise-allwise-cat-part25','AllWISE/wise-allwise-cat-part34', 'AllWISE/wise-allwise-cat-part35', 'AllWISE/wise-allwise-cat-part39']

for gcut in [15]:
    for sgcut in [5]:
        for tmrscut in [0.00028]:
            for tmzcut in [0.03]:
                for diffcuts in [1]:
                    for jcut in [0]:
                        for pv in [400]:
                            for lcut in [0]:
                                for ecut in [1]:
                                    for fgo in flist:
                                        fname = fgo + "_gbcut"+str(gcut)
                                        if diffcuts:
                                                fname = fname+"_gselect"
                                        if jcut:
                                                fname = fname+"_jmcut"
                                        if sgcut:
                                                fname = fname+"_sgbcut"+str(sgcut)
                                        fname = fname+'_2mrscut'+str(tmrscut)+'_2mrszcut'+str(tmzcut)
                                        if lcut:
                                                fname = fname+"_lrem"
                                        if pv:
                                                fname = fname+"_vcut"+str(pv)
                                        if ecut:
                                                fname = fname+"_extcut"
                                        fname=fname+'_galaxyselection.txt'
                                        f = glob(fname)
                                        print f
                                        if len(f):
                                            continue
                                        jobname = jname+fgo[29:]+str(gcut)+str(sgcut)+str(tmzcut)+str(tmrscut)+str(diffcuts)+str(jcut)+str(pv)+str(jcut)+str(ecut)+'job'
                                        fout = open(jobname+'.slurm', "w")
                                        jobline = 'python processor_WISE.py -f AllWISE/'+' -i '+fgo+' -g '+str(gcut)+' -s '+str(sgcut)+' -r 32 -o -m '+str(tmrscut)+' -z '+str(tmzcut) + ' -v '+str(pv)
                                        if diffcuts:
                                            jobline = jobline + ' -d '
                                        if jcut:
                                            jobline = jobline + ' -j '
                                        if lcut:
                                            jobline = jobline + ' -l '
                                        if ecut:
                                            jobline = jobline + ' -e '                                           
                                        for line in flines:
                                            fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
                                        fout.close()
                                        os.system('chmod +x ' + jobname+'.slurm')
                                        os.system('sbatch -p icecube '+ jobname+'.slurm')
                                        time.sleep(0.01)
                                        raw_input("Press Enter to Continue")
