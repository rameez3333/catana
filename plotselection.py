import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from glob import glob

flist = sorted(glob("AllWISE/wise-allwise-cat-part??_gbcut15_gselect_sgbcut5_2mrscut0.00028_2mrszcut0.02_lrem_galaxyselection.txt"))

arr = np.genfromtxt(flist[0], delimiter="|")


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
    ax.scatter(np.radians(x),np.radians(Dec), alpha=0.1, s=0.003)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True, which='both')


for fname in flist[1:]:
    print fname
    arr = np.append(arr, np.genfromtxt(fname, delimiter="|"), axis=1)





#data = np.log10(totalmotion)

allarr=arr

lra = [202.5,  241.3, 247.05, 194.95, 80.893, 246.3, 165.3, 285.4, 10.7]
ldec = [-31.0, 17.75, 40.93, 27.98, -69.756, -24.5, -78.2, -37.47, 40.48]
lwidth = [4., 3.,     3, 3, 5, 5, 2, 2, 2]

print 'Removing Around Objects'
for ra, dec, width in zip(lra, ldec, lwidth):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(allarr[1]))*np.cos(dec)*np.cos(np.deg2rad(allarr[0]) - ra)+np.sin(np.deg2rad(allarr[1]))*np.sin(dec)))
    allarr = allarr.transpose()[dangle>width].transpose()
print 'Removing diametrically opposite around objects'
for ra, dec, width in zip(lra, ldec, lwidth):
    ra, dec = np.deg2rad(ra-180.), np.deg2rad(-1.*dec)
    if ra<0:
        ra = ra+2.*np.pi
    dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(allarr[1]))*np.cos(dec)*np.cos(np.deg2rad(allarr[0]) - ra)+np.sin(np.deg2rad(allarr[1]))*np.sin(dec)))
    allarr = allarr.transpose()[dangle>width].transpose()   
countlrem = len(allarr[0])
print countlrem, 'objects survive specific superclusters removal'

arr=allarr

totalmotion = np.rad2deg(np.power((np.power((np.cos(np.deg2rad(arr[1]))*np.deg2rad(arr[17]/3.6e6)), 2) + np.power(np.deg2rad(arr[18]/3.6e6), 2)), 0.5))*3.6e6

totalmotion[np.isnan(totalmotion)]=0

#arr = arr.transpose()[totalmotion<200].transpose()

#np.savetxt('wise-allwise-galaxyselection_jcut_merged.txt',arr.transpose(), delimiter="|")

#arr = arr.transpose()[(arr[2]<-20.) + (arr[2]>20.)].transpose()

print len(arr[0])

#plot_mwd(arr[0], arr[1], org=180., title = "AllWISE")

#hp.projscatter(np.deg2rad(90.-arr[1]), np.deg2rad(arr[0]))

binwidth=10
bins = range(int(min(totalmotion-binwidth)), 1000, binwidth)
plt.hist(arr[3], bins=bins)
plt.xlabel('Proper Motion [mas/yr]', fontsize=15.)
plt.ylabel('Number of Sources', fontsize=15.)
plt.show()
#plt.savefig('AllWISEbestmincrwfrdlrem.png')
