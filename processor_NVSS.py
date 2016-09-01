import numpy as np
import healpy as hp
import sys, os
from astropy import units as u
from astropy.coordinates import SkyCoord
from optparse import OptionParser
import astropy.io.fits as fits
import matplotlib.pyplot as plt

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

#parser.add_option("-i", "--input", action="store", type="string", default="", dest="INPUT", help="Input i3 file to process")
parser.add_option("-f", "--folder", action="store", type="string", default="", dest="FOLDER", help="Folder")
parser.add_option("-g", "--gallatcut", action="store", type="int", default=10., dest="GCUT", help="Width around the galactic plane to cut")
parser.add_option("-s", "--sgallatcut", action="store", type="int", default=0., dest="SGCUT", help="Width around the supergalactic plane to cut")
parser.add_option("-r", "--resolution", action="store", type="int", default=1, dest="RES", help="Resolution of the scan, in healpy nside") 
parser.add_option("-t", "--thresholdflux", action="store", type="int", default=1, dest="TF", help="Threshold of the flux")      
parser.add_option("-d", "--datfile", action = "store_true", default=False, dest="DATFILE", help = "Whether to use the ascii file instead of the fits table")
#parser.add_option("-j", "--jmcut", action = "store_true", default=False, dest="JMCUT", help = "Whether to apply the j2mass > 15.8 cut")
parser.add_option("-o", "--savesel", action = "store_true", default=False, dest="SOUT", help = "Whether to save a selection")

(options, args) = parser.parse_args()

scanres = options.RES
#scanfile = options.INPUT
outfolder = options.FOLDER
gallatcut = options.GCUT
sgallatcut = options.SGCUT
useascii = options.DATFILE
#jmcut = options.JMCUT
fluxcut = options.TF


npix = hp.nside2npix(scanres)

if not useascii:
    fin = fits.open('NVSS/CATALOG.FIT')
    din = fin[1].data
    NVarr = np.vstack((din['RA(2000)'], din['DEC(2000)'], din['PEAK INT']*1000.))
    allarr = NVarr
else:
    inarr = np.genfromtxt('NVSS/NVSSNoBullshit.txt', usecols=(0,1,2,3,4,5,7), skip_footer=1).transpose()
    ra = (inarr[0] + inarr[1]/60. + inarr[2]/3600.)/24.*360.
    dec = (inarr[3] + np.sign(inarr[3])*inarr[4]/60. + np.sign(inarr[3])*inarr[5]/3600.)
    allarr = np.vstack((ra, dec, inarr[6]))

#print "Loading file"

#def nulltonan(val):
    #if val=='NULL' or val=='null':
        #return np.nan
    #else:
        #return float(val)

def SpaceAngleRad(zen1, azi1, zen2, azi2):
    return acos( sin(zen1)*sin(zen2)*cos(azi1-azi2)+cos(zen1)*cos(zen2))

#if 'allwise' in scanfile:
    #allarr = np.genfromtxt(scanfile, usecols=(1,2,7,16,286,287 ), delimiter='|', converters={7:nulltonan,16:nulltonan,286:nulltonan,287:nulltonan}).transpose()
#else:
    #allarr = np.genfromtxt(scanfile, usecols=(1,2,7,16,272,273 ), delimiter='|', converters={7:nulltonan,16:nulltonan,272:nulltonan,273:nulltonan}).transpose()


decmin, decmax = np.min(allarr[1]), np.max(allarr[1])
totno = len(allarr[0])


print totno, ' objects in declination range ', decmin, ' - ', decmax, ' loaded'

outfile = 'nvssonly'

if useascii:
    outfile=outfile+'ascii'

gallatcutno = 0
n2masscutno = 0
w1mprocutno = 0
wmproj2masscutno = 0
j2masscutno = 0
sgbcutno = 0
finalno = 0

if gallatcut:
    print "Applying Gal Lat cut, +/-", gallatcut 
    gb = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).galactic.b.value
    allarr = allarr.transpose()[(gb<-1.*gallatcut) + (gb>gallatcut)].transpose()
    gallatcutno = len(allarr[0])
    print gallatcutno, ' objects survive Galactic Latitude cut'
    outfile = outfile+"_gbcut"+str(gallatcut)
    
print "Applying Flux Cut", fluxcut
allarr = allarr.transpose()[(allarr[2]>fluxcut)*(allarr[2]<1000.)].transpose()
fcno = len(allarr[0])
print fcno, ' objects survive flux cut (incl <1000mJy cut)'
outfile = outfile+"_fluxcut"+str(fluxcut)

print "Applying declination Cut <", 40.
allarr = allarr.transpose()[(allarr[1]<40.)*(allarr[1]>-40.)].transpose()
dcno = len(allarr[0])
print dcno, ' objects survive declinationcut'

#if restselectioncut:
    #print "Now applying N2Mass Cut"
    #allarr = allarr.transpose()[allarr[4]>0].transpose()
    #n2masscutno = len(allarr[0])
    #print n2masscutno, ' objects survive n_2mass cut'
    #print "Now applying w1mpro Cut"
    #allarr = allarr.transpose()[(allarr[3]>12.0) * (allarr[3]<15.2)].transpose()
    #w1mprocutno = len(allarr[0])
    #print w1mprocutno, ' objects survive w1mpro cut'
    #print "Now applying w1mpro - jm2mass Cut"
    #allarr = allarr.transpose()[(allarr[3]-allarr[5])<-1.7].transpose()
    #wmproj2masscutno = len(allarr[0])
    #print wmproj2masscutno, ' objects survive w1mpro -jm2mass cut'
    #print "Now applying jm2mass cuts"
    #allarr = allarr.transpose()[allarr[5]<16.5].transpose()
    #outfile = outfile+"_gselect"
    
#if jmcut:
    #allarr = allarr.transpose()[allarr[5]>15.8].transpose()
    #j2masscutno = len(allarr[0])
    #print j2masscutno, ' objects survive j2mass cut'
    #outfile = outfile+"_jmcut"

if sgallatcut:
    print "Now applying supergalactic latitude cut"
    sgb = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).supergalactic.sgb.value
    allarr = allarr.transpose()[(sgb<-1.*sgallatcut) + (sgb>sgallatcut)].transpose()
    sgbcutno = len(allarr[0])
    print sgbcutno, 'objects survive Supergalactic latitude cut'
    outfile = outfile+"_sgbcut"+str(sgallatcut)

if options.SOUT:
    np.savetxt(outfile+'_galaxyselection.txt',allarr, delimiter="|")

finalno = len(allarr[0])


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
    ax.scatter(np.radians(x),np.radians(Dec), alpha=0.6, s=0.001)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)


#plot_mwd(allarr[0], allarr[1], org=180., title = "NVSS Raw, +/- 5deg b cut, +/-5deg supergalactic b cut")

#plt.show()

fout = open(outfile+'_result.txt', "w")
fout.write(str(decmin)+ "," + str(decmax)+"\n")
fout.write(str(totno)+","+str(gallatcutno)+","+str(n2masscutno)+","+str(w1mprocutno)+","+str(wmproj2masscutno)+","+str(j2masscutno)+","+str(sgbcutno)+","+str(finalno)+"\n")

for i in range(0, npix):
    dec, ra = np.deg2rad((90. - np.rad2deg(hp.pix2ang(scanres, i)[0]))), hp.pix2ang(scanres, i)[1]
    print "DEC. ", (90. - np.rad2deg(hp.pix2ang(scanres, i)[0])), " R.A.", np.rad2deg(hp.pix2ang(scanres, i)[1])
    dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(allarr[1]))*np.cos(dec)*np.cos(np.deg2rad(allarr[0]) - ra)+np.sin(np.deg2rad(allarr[1]))*np.sin(dec)))
    print "UH: ", len(allarr.transpose()[dangle<90.]), "LH: ", len(allarr.transpose()[dangle>90.])
    fout.write(str(np.rad2deg(dec))+","+str(np.rad2deg(ra))+","+str(len(allarr.transpose()[dangle<90.]))+","+str(len(allarr.transpose()[dangle>90.]))+"\n")

fout.close()
os.system('mv '+ outfile+'_result.txt '+outfolder)
os.system('mv '+ outfile+'_galaxyselection.txt '+outfolder)
