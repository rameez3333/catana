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



fout = open("AllWISEcrwford2mrsz03lrem.txt", "w")
fout.write("Angular tolerance within which 2MRS sources are looked for|Redshift cut on the 2MRS catalog |Sourced removed because of 2MRS overlap|Supergalactic latitude cut(+/-)|Galactic Latitude Cut(+/-)|Galaxy Selection Cuts?|Magnitude Cut?|Total Number of sources left after cuts|dipole|vcut|lcut|ecut|DEC(deg)|dipole RA(deg)|dipole value|Galactic b (deg)|Galactic(l)(deg)|Radio eq Velocity(km/s)|Angle to CMB dipole\n")

flist = sorted(glob('AllWISE/wise-allwise-cat-part??'))

c = 299792.458

counter=0
for gcut in [15]:
    for sgcut in [5]:
        for tmrscut in [0.00028]:
            for tmzcut in [0.03]:
                for diffcuts in [1]:
                    for jcut in [0]:
                        for pv in [0, 200, 400, 600]:
                            for lcut in [1,0]:
                                for ecut in [0, 1]:
#                                    try:
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
                                            fname = fname+'_2mrscut'+str(tmrscut)+'_2mrszcut'+str(tmzcut)
                                            if lcut:
                                                fname = fname+"_lrem"
                                            if pv:
                                                fname = fname+"_vcut"+str(pv)
                                            if ecut:
                                                fname = fname+"_extcut"
                                            fname=fname+'_cwrfrd_result.txt'
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
                                        fout.write(str(tmrscut)+"|"+str(tmzcut)+"|"+str(tmrssources)+"|"+str(sgcut)+"|"+str(gcut)+"|"+str(diffcuts)+"|"+str(jcut)+"|"+str(totsources)+"|"+str(pv)+"|"+str(lcut)+"|"+str(ecut)+"|"+str(dec)+"|"+str(ra)+"|"+str(resdip/totsources)+"|"+str(gb)+"|"+str(gl)+"|"+str(v)+'|'+str(dangle)+"\n")
                                        counter+=1
 #                                   except:
  #                                      print "Fuck this"

print "Total ", counter
            


            
