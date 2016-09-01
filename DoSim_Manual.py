import numpy as np
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.modeling.rotations as rot
from optparse import OptionParser
import healpy as hp
import sys, os

usage = 'usage: %prog [options]'
parser = OptionParser(usage)


parser.add_option("-r", "--ra", action="store", type="float", default=120.0, dest="RA", help="R.A of the velocity direction in decimal degrees")                  
parser.add_option("-d", "--dec", action = "store", type = "float", default=-30.0, dest="DEC", help = "Declination of the velocity direction in decimal degrees")
parser.add_option("-v", "--velocity", action = "store", type = "float", default=1000.0, dest="VEL", help = "Velocity in km/s")
parser.add_option("-n", "--number", action = "store", type = "int", default=100000, dest="NS", help = "Number of sources in catalog")
parser.add_option("-s", "--serial", action = "store", type = "int", default=0, dest="SER", help = "Unique Identifier")

(options, args) = parser.parse_args()


raab = options.RA
decab = options.DEC
velocab = options.VEL
ncat = options.NS
serno = options.SER

fhead = str(serno)+'_Sim_'+str(raab)+'_'+str(decab)+'_'+str(velocab)+'_'+str(ncat)+'_'

c = 299792.458

fout = open(fhead+'summary.txt',"w")

fout.write("Input RA(deg)|Input Dec(deg)|Input velocity(km/s)|Size of random catalog|Trial number\n")
fout.write(str(raab) + '|'+str(decab)+'|'+str(velocab)+'|'+str(ncat)+'|'+str(serno)+'\n')

fout.write("Cut around galactic plane(+/-deg)|Cut around supergalactic plane(+/-deg)|Result RA(deg)|Result Dec(deg)|Result velocity(km/s)|Result dipole strength\n")
fout.write("Before Aberrations\n")

def GenerateDataSet(size):
    ra = np.random.uniform(0, 360, size=size)
    dec = np.rad2deg(np.arcsin(np.random.uniform(0, 1, size=size))*np.random.choice([1.,-1.], size=size))
    return np.vstack((ra, dec))
    
    
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
    
    
allarr = GenerateDataSet(ncat)


def aberrator(arr, ra, dec, v):
    gamma = 1./np.sqrt(1-(v/c)**2.)
    #print gamma, v/c
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    dangle = np.arccos(np.cos(np.deg2rad(arr[1]))*np.cos(dec)*np.cos(np.deg2rad(arr[0]) - ra)+np.sin(np.deg2rad(arr[1]))*np.sin(dec))
    #print np.rad2deg(np.min(dangle)), np.rad2deg(np.max(dangle))
    tanphi = np.sin(dangle)/(gamma*(np.cos(dangle)-v/c))
    rotor = rot.RotateCelestial2Native(lon = np.rad2deg(ra), lat=np.rad2deg(dec), lon_pole=5.*u.degree)
    rotorback = rot.RotateNative2Celestial(lon = np.rad2deg(ra), lat=np.rad2deg(dec), lon_pole=5.*u.degree)
    lon, lat = rotor(arr[0], arr[1])
    #print np.rad2deg(dangle) + lat
    #print np.rad2deg(dangle)
    #print np.rad2deg(np.arctan(tanphi))
    phi = np.rad2deg(np.arctan(tanphi))
    phi = phi+(-1.*np.sign(phi) + 1.)*90.
    #print phi
    #print dangle - np.arctan(tanphi)
    diff = np.rad2deg(dangle) - phi
    print diff
    #latab = np.rad2deg(np.arctan(tanphi))
    latab = lat - diff
    
    #print (np.rad2deg(dangle) - np.rad2deg(np.arctan(tanphi)))#[0], lon[0], lat[0]
    
    #print lat - latab
    #print lat
    cellon, cellat = rotorback(lon, latab)
    return np.vstack((cellon, cellat))
    





def Processor(allarr,nside, gcut=0., sgcut=0.,fhead='wtf'):
    npix = hp.nside2npix(nside)
    LH = np.ones(npix)
    UH = np.ones(npix)
    
    #if gcut:
        #gb = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).galactic.b.value
        #allarr = allarr.transpose()[(gb<-1.*gcut) + (gb>gcut)].transpose()
    #if sgcut:
        #sgb = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).supergalactic.sgb.value
        #allarr = allarr.transpose()[(sgb<-1.*sgcut) + (sgb>sgcut)].transpose()
    #for i in range(0, npix):
        #print i
        #dec, ra = np.deg2rad((90. - np.rad2deg(hp.pix2ang(nside, i)[0]))), hp.pix2ang(nside, i)[1]
        #print "DEC. ", (90. - np.rad2deg(hp.pix2ang(nside, i)[0])), " R.A.", np.rad2deg(hp.pix2ang(nside, i)[1])
        #dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(allarr[1]))*np.cos(dec)*np.cos(np.deg2rad(allarr[0]) - ra)+np.sin(np.deg2rad(allarr[1]))*np.sin(dec)))
        #print "UH: ", len(allarr.transpose()[dangle<90.]), "LH: ", len(allarr.transpose()[dangle>90.])
        #LH[i], UH[i] = len(allarr.transpose()[dangle>90.]), len(allarr.transpose()[dangle< 90.])
    plot_mwd(allarr[0], allarr[1], org=0., title = "Monte Carlo")
    plt.show()
    plt.savefig(fhead+'_gcut_'+str(gcut)+'_sgcut_'+str(sgcut)+'_Cat.pdf')
    return LH, UH    


def Rescompiler(LH, UH, gcut=0., sgcut=0., fhead='wtf'):
    hmap = (UH-LH)/(UH+LH)
    v = hmap[np.argmax(hmap)]*2./(2.)*299792.458
    print (90. - np.rad2deg(hp.pix2ang(32, np.argmax(hmap))[0])), np.rad2deg(hp.pix2ang(32, np.argmax(hmap))[1]), hmap[np.argmax(hmap)], v
    hp.mollview(hmap)
    plt.savefig(fhead+'_gcut_'+str(gcut)+'_sgcut_'+str(sgcut)+'_Scan.pdf')
    return (90. - np.rad2deg(hp.pix2ang(32, np.argmax(hmap))[0])), np.rad2deg(hp.pix2ang(32, np.argmax(hmap))[1]), hmap[np.argmax(hmap)], v



#for cuts in [[0,0]]:#,[5,0],[10,0],[20,0],[10,5],[10,10], [20,5], [20,10]]:
    #LH,UH = Processor(allarr, 32, cuts[0], cuts[1], fhead)
    #resra,resdec,resdip, resvel = Rescompiler(LH, UH, cuts[0],cuts[1], fhead)
    ##fout.write(str(cuts[0])+'|'+str(cuts[1])+'|'+str(resra)+'|'+str(resdec)+'|'+ str(resvel)+'|'+str(resdip)+'\n')

allarr=aberrator(allarr, raab, decab, velocab)
plot_mwd(allarr[0], allarr[1], org=0., title = "After aberration: RA 125.0, DEC -30.0 v=c/2 Monte Carlo Catalogsize: "+str(ncat))

plt.show()
               
#allarr=aberrator(allarr, raab, decab, velocab)

#fhead = fhead+'_postab'

#fout.write('After Aberration applied\n')

#for cuts in [[0,0]]:#,[5,0],[10,0],[20,0],[10,5],[10,10], [20,5], [20,10]]:
    #LH,UH = Processor(allarr, 32, cuts[0], cuts[1], fhead)
    #resra,resdec,resdip, resvel = Rescompiler(LH, UH, cuts[0],cuts[1], fhead)
    ##fout.write(str(cuts[0])+'|'+str(cuts[1])+'|'+str(resra)+'|'+str(resdec)+'|'+ str(resvel)+'|'+str(resdip)+'\n')

#fout.close()
#os.system('mv '+fhead.replace('_postab','*.*')+' Sim/')
