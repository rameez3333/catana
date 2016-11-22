import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import healpy as hp
import os, sys



twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
twomrsarr = np.vstack((twomrsarr[0], twomrsarr[1], twomrsarr[2]/299792.458))
starburstarr = np.genfromtxt('Starburst/20mJy_sample_181208.csv', skip_header=1, usecols=(1,2,3)).transpose()
AGNS = np.genfromtxt('AGNs/CarameteBiermanncatalog.csv', skip_header=1, usecols=(1,2,6), delimiter=",").transpose()
agncoord =SkyCoord(l = AGNS[0]*u.degree, b = AGNS[1]*u.degree, frame='galactic')
AGNS = np.vstack((agncoord.icrs.ra.value, agncoord.icrs.dec.value, AGNS[2]))

def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x = np.remainder(RA+360-org,360) # shift RA values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
    ax.scatter(np.radians(x),np.radians(Dec), alpha=0.5, s=5)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
    
    
twomrsarr = twomrsarr.transpose()[twomrsarr[2]<0.01].transpose()
starburstarr = starburstarr.transpose()[starburstarr[2]<0.02].transpose()
AGNS = AGNS.transpose()[AGNS[2]<0.02].transpose()

plot_mwd(twomrsarr[0], twomrsarr[1], title="2MRS, z<0.01")
plt.show()


def scattomap(dec,ra, nside=16):
    #hmap = np.zeros(hp.nside2npix(nside))
    hmap = np.bincount(hp.ang2pix(nside, np.deg2rad(90.-dec), np.deg2rad(ra)), minlength=hp.nside2npix(nside))
    return hmap
    
    
twomap = scattomap(twomrsarr[1], twomrsarr[0])
starburstmap = scattomap(starburstarr[1], starburstarr[0])
AGNmap = scattomap(AGNS[1], AGNS[0])

randmap = np.random.rand(hp.nside2npix(16))

a = hp.sphtfunc.anafast(map1 = twomap, map2=twomap)
b  = hp.sphtfunc.anafast(map1 = twomap, map2=starburstmap)
c = hp.sphtfunc.anafast(map1 = twomap, map2=AGNmap)
d = hp.sphtfunc.anafast(map1 = twomap, map2=randmap)

print a

print b

print c

print d
hp.mollview(twomap, title="2MRS , z<0.01")
plt.show()
    
