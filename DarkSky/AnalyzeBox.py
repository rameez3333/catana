import yt
import numpy as np
from astropy.cosmology import Planck15 as cosmo
from scipy.interpolate import InterpolatedUnivariateSpline
from optparse import OptionParser


inv_comoving_distance = InterpolatedUnivariateSpline(cosmo.comoving_distance(np.linspace(0, 0.6, 1000)).value, np.linspace(0, 0.6, 1000))

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-x", "--xcenter", action="store", type="float", default=0.0, dest="XC", help="Xcoord of the Center")
parser.add_option("-y", "--ycenter", action="store", type="float", default=0.0, dest="YC", help="Ycoord of the Center")
parser.add_option("-z", "--zcenter", action="store", type="float", default=0.0, dest="ZC", help="Zcoord of the Center")
parser.add_option("-w", "--width", action="store", type="float", default=4000., dest="WIDTH", help="Width of the box, radius of the circle, in MPc")

(options, args) = parser.parse_args()

center = np.array([options.XC,  options.YC, options.ZC])*0.688062

width = options.WIDTH*0.688062

outfname = 'DarkSkyDipoleQuery_CenterX_'+str(options.XC)+'_Y_'+str(options.YC)+'_Z_'+str(options.ZC)+'_W_'+str(width)+'_Out.txt'

bbox = np.array([center-(width+100.*0.688062)/2., center+(width+100.*0.688062)/2.])


ds = yt.load('Halos/ds14_a_halos_1.0000',midx_filename='Halos/ds14_a_halos_1.0000.midx9',bounding_box=bbox)
 
ad = ds.sphere(center, (100.0, 'Mpccm/h'))


X = np.asarray(ad['x'].in_units('Mpc'))
Y = np.asarray(ad['y'].in_units('Mpc'))
Z = np.asarray(ad['z'].in_units('Mpc'))

VX = np.asarray(ad['vx'].in_units('km/s'))
VY = np.asarray(ad['vy'].in_units('km/s'))
VZ = np.asarray(ad['vz'].in_units('km/s'))

M = np.asarray(ad['m200c'])



xmed, ymed, zmed = np.median(X), np.median(Y), np.median(Z)

print 'Median:', xmed, ymed, zmed

disttomed = np.power((X-xmed), 2) + np.power((Y-ymed), 2) + np.power((Z-zmed), 2)

observerindex = np.argmin(disttomed)

xob, yob, zob = X[observerindex], Y[observerindex], Z[observerindex]

center = np.array([xob,  yob, zob])

ad = ds.sphere(center, (width/2., 'Mpccm/h'))

X = np.asarray(ad['x'].in_units('Mpc'))
Y = np.asarray(ad['y'].in_units('Mpc'))
Z = np.asarray(ad['z'].in_units('Mpc'))

VX = np.asarray(ad['vx'].in_units('km/s'))
VY = np.asarray(ad['vy'].in_units('km/s'))
VZ = np.asarray(ad['vz'].in_units('km/s'))

M = np.asarray(ad['m200c'])

disttoobsq = np.power((X-xob), 2) + np.power((Y-yob), 2) + np.power((Z-zob), 2)

observerindex = np.argmin(disttoobsq)

fluxatob = M/disttoobsq

nhalos = len(X)

disttoob = np.sqrt(disttoobsq)

print 'Maximum Distance to Observer', np.nanmax(disttoob)

redshift = inv_comoving_distance(disttoob)

halarr = np.vstack([X-xob,Y-yob,Z-zob, VX,VY,VZ, M, fluxatob, redshift, disttoob])

observer = halarr.transpose()[observerindex]

fout = open(outfname, 'w')

print 'Observer', observer
fout.write('Number of Halos: '+str(nhalos)+'\n')
fout.write('Observer: '+str(xob)+'|'+str(yob)+'|'+str(zob)+'|'+str(observer[3])+'|'+str(observer[4])+'|'+str(observer[5])+'|'+str(observer[6])+'\n' )

fout.write('Flux Limited Trials:\n')

halarr = np.delete(halarr.transpose(),observerindex, axis=0).transpose()

halarr = halarr.transpose()[halarr[6]>0].transpose()




#find center


fluxthreshes = np.power(np.linspace(np.log10(np.min(halarr[7])), np.log10(np.max(halarr[7])), 5000, ), 10.)

zhistbins = np.linspace(0, 0.6, 61)


for ft in fluxthreshes:
    halarr = halarr.transpose()[halarr[7]>ft].transpose()
    ncuthalos=len(halarr[0])
    cx = np.sum(halarr[0]/halarr[9])
    cy = np.sum(halarr[1]/halarr[9])
    cz = np.sum(halarr[2]/halarr[9])
    
    zchalarr = halarr.transpose()[halarr[8]>0.03].transpose()
    zncuthalos=len(zchalarr[0])
    
    zcx = np.sum(zchalarr[0]/zchalarr[9])
    zcy = np.sum(zchalarr[1]/zchalarr[9])
    zcz = np.sum(zchalarr[2]/zchalarr[9])
    
    zhist = np.histogram(halarr[8], bins=zhistbins)[0]
    strhist=''
    for b in zhist:
        strhist = strhist+str(b)+','
    fout.write(str(ft)+'|'+str(ncuthalos)+'|'+str(cx)+'|'+str(cy)+'|'+str(cz)+'|'+str(zncuthalos)+'|'+str(zcx)+'|'+str(zcy)+'|'+str(zcz)+'|'+strhist[:len(strhist)-1]+'\n')
    
fout.close()








