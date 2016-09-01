import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
from glob import glob

foldname = sys.argv[1]
tagname = sys.argv[2]

flist = sorted(glob(foldname+'*'+tagname+"_galaxyselectionnew.txt"))
print len(flist), 'files found'

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
    ax.scatter(np.radians(x),np.radians(Dec), alpha=0.5, s=0.005)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)


#for fname in flist[1:]:
    #print fname
    #arr = np.append(arr, np.genfromtxt(fname, delimiter="|", skip_footer=4), axis=1)

#arr = arr.transpose()[(arr[2]<-20.) + (arr[2]>20.)].transpose()

print len(arr[0])

plot_mwd(arr[0], arr[1], org=180., title = "NVSS+SUMSS")

#hp.projscatter(np.deg2rad(90.-arr[1]), np.deg2rad(arr[0]))
plt.show()
#plt.savefig(tagname+'scatter.png')
