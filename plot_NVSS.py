import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
from glob import glob
import astropy.io.fits as fits


pdec = float(sys.argv[1])



inarr = np.genfromtxt('NVSS/NVSSNoBullshit.txt', usecols=(0,1,2,3,4,5,7), skip_footer=1).transpose()
ra = (inarr[0] + inarr[1]/60. + inarr[2]/3600.)/24.*360.
dec = (inarr[3] + np.sign(inarr[3])*inarr[4]/60. + np.sign(inarr[3])*inarr[5]/3600.)
NVarr = np.vstack((ra, dec, inarr[6]))

#NVarr = np.vstack((din['RA(2000)'], din['DEC(2000)']))
SUin = np.genfromtxt('NVSS/SUMSS.txt', usecols = (0,1,2,3,4,5,8)).transpose()
SUra = (SUin[0] + SUin[1]/60. + SUin[2]/3600.)/24.*360.
SUdec = (SUin[3] - SUin[4]/60. - SUin[5]/3600.)
SUarr = np.vstack((SUra, SUdec, SUin[6]))



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
    ax.scatter(np.radians(x),np.radians(Dec), alpha=0.1, s=0.001)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)


NVcompdeccount = len(NVarr.transpose()[NVarr[1]> -1.*pdec].transpose()[0])

print NVcompdeccount, "sources found in NVSS in complementary declination region"

SUdeccount = len(SUarr.transpose()[SUarr[1]< pdec].transpose()[0])

print SUdeccount, "sources found in SUMSS in declination region"



totarr = np.append(NVarr, SUarr, axis=1)

#arr = arr.transpose()[(arr[2]<-20.) + (arr[2]>20.)].transpose()

print len(totarr[0])

plot_mwd(totarr[0], totarr[1], org=180., title = "NVSS+SUMSS Raw")

#hp.projscatter(np.deg2rad(90.-arr[1]), np.deg2rad(arr[0]))
plt.show()
plt.savefig(tagname+'scatter.png')
