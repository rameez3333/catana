import numpy as np
from glob import glob
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import sys

#nside = int(sys.argv[1])
##glcut = sys.argv[2]

#npix = hp.nside2npix(nside)



fout = open("AllWISEcrwford.txt", "w")
fout.write("Angular tolerance within which 2MRS sources are looked for|Redshift cut on the 2MRS catalog |Sourced removed because of 2MRS overlap|Supergalactic latitude cut(+/-)|Galactic Latitude Cut(+/-)|Galaxy Selection Cuts?|Magnitude Cut?|Total Number of sources left after cuts|dipole DEC(deg)|dipole RA(deg)|dipole value|Galactic b (deg)|Galactic(l)(deg)|Radio eq Velocity(km/s)|Angle to CMB dipole\n")

flist = sorted(glob('AllWISE/wise-allwise-cat-part??'))

c = 299792.458

counter=0
for gcut in [10, 15, 20]:
    for sgcut in [0,5,10]:
        for tmrscut in [0.00028, 0.0028, 0.01]:
            for tmzcut in [0.02, 0.05, 0.09]:
                for diffcuts in [1, 0]:
                    for jcut in [1, 0]:
                        totsources=0.
                        cx,cy,cz = 0.,0.,0.
                        tmrssources = 0.
                        for fgo in flist:
                            fname = fgo + "_gbcut"+str(gcut)
                            if diffcuts:
                                fname = fname+"_gselect"
                            if jcut:
                                fname = fname+"_jmcut"
                            if sgcut:
                                fname = fname+"_sgbcut"+str(sgcut)
                            fname = fname+'_2mrscut'+str(tmrscut)+'_2mrszcut'+str(tmzcut)+'_cwrfrd_result.txt'
                            print fname
                            flines = open(fname).readlines()
                            print "Dec range:", flines[0]
                            totsources = totsources + float(flines[1].replace('\n','').split(",")[-1]) 
                            l2 = np.asarray([float(x) for x in flines[1].replace('\n','').split(",")])
                            tmrssources =  tmrssources + l2[-2] - l2[np.nonzero(l2)[0][-3]]
                            
                            l3 = flines[2].replace('\n','').split(",")
                            cx = cx+float(l3[0])
                            cy = cy+float(l3[1])
                            cz = cz+float(l3[2])
                    
                        ske = SkyCoord(x=cx, y=cy, z=cz, representation='cartesian').frame.represent_as('spherical')
                        ra = ske.lon.degree
                        dec = ske.lat.degree
                        resdip = np.sqrt(cx*cx + cy*cy + cz*cz)
                        vel = c*resdip*3./(totsources*(2. + 1.05*(1+0.75)))
                        gb = SkyCoord(ra = ra*u.degree, dec=dec*u.degree).galactic.b.value
                        gl = SkyCoord(ra = ra*u.degree, dec=dec*u.degree).galactic.l.value
                        v = vel
                        dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(-7.0))*np.cos(np.deg2rad(dec))*np.cos(np.deg2rad(168.0) - np.deg2rad(ra))+np.sin(np.deg2rad(-7.0))*np.sin(np.deg2rad(dec))))
                        fout.write(str(tmrscut)+"|"+str(tmzcut)+"|"+str(tmrssources)+"|"+str(sgcut)+"|"+str(gcut)+"|"+str(diffcuts)+"|"+str(jcut)+"|"+str(totsources)+"|"+str(dec)+"|"+str(ra)+"|"+str(resdip/totsources)+"|"+str(gb)+"|"+str(gl)+"|"+str(v)+'|'+str(dangle)+"\n")
                        counter+=1

print "Total ", counter
            


            
