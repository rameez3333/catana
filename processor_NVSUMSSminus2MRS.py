import numpy as np
import healpy as hp
import sys, os
from astropy import units as u
from astropy.coordinates import SkyCoord
from optparse import OptionParser
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import k3match

c = 299792.458

usage = 'usage: %prog [options]'
parser = OptionParser(usage)

#parser.add_option("-i", "--input", action="store", type="string", default="", dest="INPUT", help="Input i3 file to process")
parser.add_option("-f", "--folder", action="store", type="string", default="", dest="FOLDER", help="Folder")
parser.add_option("-g", "--gallatcut", action="store", type="int", default=10., dest="GCUT", help="Width around the galactic plane to cut")
parser.add_option("-s", "--sgallatcut", action="store", type="int", default=0., dest="SGCUT", help="Width around the supergalactic plane to cut")
parser.add_option("-r", "--resolution", action="store", type="int", default=1, dest="RES", help="Resolution of the scan, in healpy nside") 
parser.add_option("-t", "--thresholdflux", action="store", type="int", default=1, dest="TF", help="Threshold of the flux")      
parser.add_option("-p", "--patchdeclination", action="store", type="float", default=40.0, dest="PD", help="Declination value at which the two catalogs are to be patched together")
parser.add_option("-m", "--twomrsoverlap", action="store", type="float", default=0.0, dest="TW", help="Angular separation to look for 2MRS catalog objects within")      
parser.add_option("-d", "--datfile", action = "store_true", default=False, dest="DATFILE", help = "Whether to use the ascii file instead of the fits table")
#parser.add_option("-j", "--jmcut", action = "store_true", default=False, dest="JMCUT", help = "Whether to apply the j2mass > 15.8 cut")
parser.add_option("-o", "--savesel", action = "store_true", default=False, dest="SOUT", help = "Whether to save a selection")
parser.add_option("-c", "--crawford", action = "store_true", default=False, dest="CRWF", help = "Use the Crawford estimator instead of the hemispherical count")

(options, args) = parser.parse_args()

scanres = options.RES
#scanfile = options.INPUT
outfolder = options.FOLDER
gallatcut = options.GCUT
sgallatcut = options.SGCUT
useascii = options.DATFILE
#jmcut = options.JMCUT
fluxcut = options.TF
patchdec = options.PD
tmrsovang = options.TW
crwford = options.CRWF

if gallatcut<10.:
    print "Warning: SUMSS already includes a cut on the Galactic latitude at +/-10."

npix = hp.nside2npix(scanres)

if not useascii:
    print "Loading from FIT file"
    fin = fits.open('NVSS/CATALOG.FIT')
    din = fin[1].data
    NVarr = np.vstack((din['RA(2000)'], din['DEC(2000)'], din['PEAK INT']*1000.))

else:
    print "Loading from the Ascii file"
    inarr = np.genfromtxt('NVSS/NVSSNoBullshit.txt', usecols=(0,1,2,3,4,5,7), skip_footer=1).transpose()
    ra = (inarr[0] + inarr[1]/60. + inarr[2]/3600.)/24.*360.
    dec = (inarr[3] + np.sign(inarr[3])*inarr[4]/60. + np.sign(inarr[3])*inarr[5]/3600.)
    NVarr = np.vstack((ra, dec, inarr[6]))

#print "Loading file"

#def nulltonan(val):
    #if val=='NULL' or val=='null':
        #return np.nan
    #else:
        #return float(val)

SUin = np.genfromtxt('NVSS/SUMSS.txt', usecols = (0,1,2,3,4,5, 8)).transpose()
SUra = (SUin[0] + SUin[1]/60. + SUin[2]/3600.)/24.*360.
SUdec = (SUin[3] - SUin[4]/60. - SUin[5]/3600.)
SUarr = np.vstack((SUra, SUdec, SUin[6]))  
        

def SpaceAngleRad(zen1, azi1, zen2, azi2):
    return acos( sin(zen1)*sin(zen2)*cos(azi1-azi2)+cos(zen1)*cos(zen2))

#if 'allwise' in scanfile:
    #allarr = np.genfromtxt(scanfile, usecols=(1,2,7,16,286,287 ), delimiter='|', converters={7:nulltonan,16:nulltonan,286:nulltonan,287:nulltonan}).transpose()
#else:
    #allarr = np.genfromtxt(scanfile, usecols=(1,2,7,16,272,273 ), delimiter='|', converters={7:nulltonan,16:nulltonan,272:nulltonan,273:nulltonan}).transpose()


ndecmin, ndecmax = np.min(NVarr[1]), np.max(NVarr[1])
sdecmin, sdecmax = np.min(SUarr[1]), np.max(SUarr[1])
ntotno = len(NVarr[0])
stotno = len(SUarr[0])

print ntotno, ' NVSS: objects in declination range ', ndecmin, ' - ', ndecmax, ' loaded'

print stotno, ' SUMSS: objects in declination range ', sdecmin, ' - ', sdecmax, ' loaded'

outfile = 'nvsumss_patch'+str(patchdec)

if useascii:
    outfile=outfile+'ascii'

sgallatcutno = 0
sn2masscutno = 0
sw1mprocutno = 0
swmproj2masscutno = 0
sj2masscutno = 0
ssgbcutno = 0
sfinalno = 0


ngallatcutno = 0
nn2masscutno = 0
nw1mprocutno = 0
nwmproj2masscutno = 0
nj2masscutno = 0
nsgbcutno = 0
nfinalno = 0

sgbcutno=0.


if gallatcut:
    print "Applying Gal Lat cut, +/-", gallatcut 
    ngb = SkyCoord(ra = NVarr[0]*u.degree, dec=NVarr[1]*u.degree).galactic.b.value
    NVarr = NVarr.transpose()[(ngb<-1.*gallatcut) + (ngb>gallatcut)].transpose()
    ngallatcutno = len(NVarr[0])
    sgb = SkyCoord(ra = SUarr[0]*u.degree, dec=SUarr[1]*u.degree).galactic.b.value
    SUarr = SUarr.transpose()[(sgb<-1.*gallatcut) + (sgb>gallatcut)].transpose()
    sgallatcutno = len(SUarr[0])    
    print ngallatcutno, sgallatcutno, ' objects survive Galactic Latitude cut'
    outfile = outfile+"_gbcut"+str(gallatcut)
    
print "Applying Flux Cut", fluxcut
NVarr = NVarr.transpose()[(NVarr[2]>fluxcut)*(NVarr[2]<1000.)].transpose()
SUarr = SUarr.transpose()[(SUarr[2]>fluxcut)*(SUarr[2]<1000.)].transpose()
nfcno, sfcno = len(NVarr[0]), len(SUarr[0])
print nfcno, sfcno, ' objects survive flux cut (incl <1000mJy cut)'
outfile = outfile+"_fluxcut"+str(fluxcut)



print "Patching NVSS and SUMSS together:"



print "NVSS : Applying declination Cut >", patchdec
NVarr = NVarr.transpose()[(NVarr[1]>patchdec)].transpose()
print "SUMSS : Applying declination Cut <", patchdec
SUarr = SUarr.transpose()[(SUarr[1]<patchdec)].transpose()


ndcno = len(NVarr[0])
print 'NVSS: ',ndcno, ' objects survive declinationcut'
sdcno = len(SUarr[0])
print 'SUMSS: ',sdcno, ' objects survive declinationcut'

ndecpatch = len(NVarr.transpose()[NVarr[1]>-1.*patchdec].transpose()[0])
sdecpatch = len(SUarr.transpose()[SUarr[1]<patchdec].transpose()[0])

print "NVSS complementary dec:", ndecpatch
print "SUMSS complementary dec:", sdecpatch


print "Removing", sdecpatch-ndecpatch, 'sources from SUMSS randomly'

ones  =  np.ones(ndecpatch)
zeros = np.zeros(sdecpatch-ndecpatch)

randomcutter = np.append(ones,zeros, axis=0)



np.random.shuffle(randomcutter)


SUarr = SUarr.transpose()[randomcutter>0.5].transpose()

print "SUMSS now has :", len(SUarr[0]), " sources"

allarr = np.append(NVarr, SUarr, axis=1)
twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2)).transpose()

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
    
    print sgbcutno, 'objects survive Supergalactic latitude cut'
    outfile = outfile+"_sgbcut"+str(sgallatcut)
sgbcutno = len(allarr[0])
counttmrs=0
if tmrsovang:
    inda,indb,cc = k3match.celestial(allarr[0], allarr[1], twomrsarr[0],twomrsarr[1], tmrsovang)
    allarr = np.delete(allarr, inda, axis=1)
    outfile=outfile+'_2mrscut'+str(tmrsovang)
    counttmrs=len(allarr[0])
    print counttmrs, 'objects survive 2mrs intersection cut'

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
if not crwford:
    fout = open(outfile+'_result.txt', "w")
    fout.write(str(ndecmin)+ "," + str(ndecmax)+','+str(sdecmin)+ "," + str(sdecmax)+"\n")
    fout.write(str(stotno)+","+str(sgallatcutno)+","+str(sn2masscutno)+","+str(sw1mprocutno)+","+str(swmproj2masscutno)+","+str(sj2masscutno)+","+str(sgbcutno)+","+str(counttmrs)+","+str(finalno)+"\n")

    for i in range(0, npix):
        dec, ra = np.deg2rad((90. - np.rad2deg(hp.pix2ang(scanres, i)[0]))), hp.pix2ang(scanres, i)[1]
        print "DEC. ", (90. - np.rad2deg(hp.pix2ang(scanres, i)[0])), " R.A.", np.rad2deg(hp.pix2ang(scanres, i)[1])
        dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(allarr[1]))*np.cos(dec)*np.cos(np.deg2rad(allarr[0]) - ra)+np.sin(np.deg2rad(allarr[1]))*np.sin(dec)))
        print "UH: ", len(allarr.transpose()[dangle<90.]), "LH: ", len(allarr.transpose()[dangle>90.])
        fout.write(str(np.rad2deg(dec))+","+str(np.rad2deg(ra))+","+str(len(allarr.transpose()[dangle<90.]))+","+str(len(allarr.transpose()[dangle>90.]))+"\n")

    fout.close()
else:
    fout = open(outfile+'_cwrfrd_result.txt', "w")
    fout.write(str(ndecmin)+ "," + str(ndecmax)+','+str(sdecmin)+ "," + str(sdecmax)+"\n")
    fout.write(str(stotno)+","+str(sgallatcutno)+","+str(sn2masscutno)+","+str(sw1mprocutno)+","+str(swmproj2masscutno)+","+str(sj2masscutno)+","+str(sgbcutno)+","+str(counttmrs)+","+str(finalno)+"\n")
    skc = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).cartesian
    cx = skc.x.value
    cy = skc.y.value
    cz = skc.z.value
    xf = cx.sum()
    yf = cy.sum()
    zf = cz.sum()

    ske = SkyCoord(x=xf, y=yf, z=zf, representation='cartesian').frame.represent_as('spherical')
    ra = ske.lon.degree
    dec = ske.lat.degree
    resdip = np.sqrt(xf*xf + yf*yf + zf*zf)
    vel = c*resdip*3./(float(len(allarr[0]))*(2. + 1.05*(1+0.75)))
    print vel
    fout.write(str(ra)+','+str(dec)+','+str(resdip)+','+str(vel)+"\n")
    fout.close()
    
os.system('mv '+ outfile+'*_result.txt '+outfolder)
os.system('mv '+ outfile+'_galaxyselection.txt '+outfolder)
