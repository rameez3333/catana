import numpy as np
import healpy as hp
import sys, os
from astropy import units as u
from astropy.coordinates import SkyCoord
from optparse import OptionParser
import k3match



usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-i", "--input", action="store", type="string", default="", dest="INPUT", help="Input i3 file to process")
parser.add_option("-f", "--folder", action="store", type="string", default="", dest="FOLDER", help="Folder")
parser.add_option("-g", "--gallatcut", action="store", type="int", default=10., dest="GCUT", help="Width around the galactic plane to cut")
parser.add_option("-s", "--sgallatcut", action="store", type="int", default=0., dest="SGCUT", help="Width around the supergalactic plane to cut")
parser.add_option("-r", "--resolution", action="store", type="int", default=1, dest="RES", help="Resolution of the scan, in healpy nside")                  
parser.add_option("-d", "--diffcuts", action = "store_true", default=False, dest="DIFFCUTS", help = "Whether to apply the remaining cuts")
parser.add_option("-j", "--jmcut", action = "store_true", default=False, dest="JMCUT", help = "Whether to apply the j2mass > 15.8 cut")
parser.add_option("-o", "--savesel", action = "store_true", default=False, dest="SOUT", help = "Whether to save a selection")
parser.add_option("-m", "--twomrsoverlap", action="store", type="float", default=0.0, dest="TW", help="Angular separation to look for 2MRS catalog objects within")   
parser.add_option("-z", "--zcut2mrs", action="store", type="float", default=0.0, dest="ZC", help="Redshift cut on the 2MRS catalog before correlation")
parser.add_option("-c", "--crawford", action = "store_true", default=False, dest="CRWF", help = "Use the Crawford estimator instead of the hemispherical count")
parser.add_option("-l", "--lrem", action = "store_true", default=False, dest="LREM", help = "Remove the list of directions")
parser.add_option("-v", "--velocitcut", action="store", type="int", default=0., dest="VELCUT", help="Velocity Cut in mas/yr")
parser.add_option("-e", "--extcut", action = "store_true", default=False, dest="EXTC", help = "Cut away extended flag sources")

(options, args) = parser.parse_args()

scanres = options.RES
scanfile = options.INPUT
outfolder = options.FOLDER
gallatcut = options.GCUT
sgallatcut = options.SGCUT
restselectioncut = options.DIFFCUTS
jmcut = options.JMCUT
tmrsovang = options.TW
zcut2mrs = options.ZC
crwford = options.CRWF
lrem =options.LREM
vel = options.VELCUT
extc = options.EXTC
npix = hp.nside2npix(scanres)

print "Loading file"

def nulltonan(val):
        if val=='NULL' or val=='null':
            return np.nan
        else:
            return float(val)
 
def SpaceAngleRad(zen1, azi1, zen2, azi2):
    return acos( sin(zen1)*sin(zen2)*cos(azi1-azi2)+cos(zen1)*cos(zen2))

if 'allwise' in scanfile:
    cols = (1,2,7,16,286,287, 3,4, 18, 20, 22, 24, 26, 28, 30, 33, 34, 45, 47, 57)
    allarr = np.genfromtxt(scanfile, usecols=cols, delimiter='|', converters={7:nulltonan,16:nulltonan,286:nulltonan,287:nulltonan, 3:nulltonan, 4:nulltonan,18:nulltonan, 20:nulltonan, 22:nulltonan, 24:nulltonan, 26:nulltonan, 28:nulltonan, 30:nulltonan, 33:nulltonan, 34:nulltonan, 45:nulltonan, 47:nulltonan, 57:nulltonan}, dtype=["float" for n in range(len(cols))])
    allarr = allarr.view((float, len(allarr.dtype.names))).transpose()
else:
    allarr = np.genfromtxt(scanfile, usecols=(1,2,7,16,272,273 ), delimiter='|', converters={7:nulltonan,16:nulltonan,272:nulltonan,273:nulltonan}).transpose()


decmin, decmax = np.min(allarr[1]), np.max(allarr[1])
totno = len(allarr[0])

twomrsarr = np.genfromtxt('2MRS/catalog/2mrs_1175_done.dat', skip_header=10, usecols=(1,2,24)).transpose()
print len(twomrsarr[0]), 'objects loaded from the 2MRS catalog' 

print totno, ' objects in declination range ', decmin, ' - ', decmax, ' loaded'

outfile = scanfile

#in order: Shapley, Hercules(2), Coma, LMC,  
lra = [202.5,  241.3, 247.05, 194.95, 80.893, 246.3, 165.3, 285.4, 10.7]
ldec = [-31.0, 17.75, 40.93, 27.98, -69.756, -24.5, -78.2, -37.47, 40.48]
lwidth = [4., 3.,     3, 3, 5, 5, 2, 2, 2]

gallatcutno = 0
n2masscutno = 0
w1mprocutno = 0
wmproj2masscutno = 0
j2masscutno = 0
sgbcutno = 0
finalno = 0

if gallatcut:
    print "Applying Gal Lat cut, +/-", gallatcut 
    allarr = allarr.transpose()[(allarr[2]<-1.*gallatcut) + (allarr[2]>gallatcut)].transpose()
    gallatcutno = len(allarr[0])
    print gallatcutno, ' objects survive Galactic Latitude cut'
    outfile = outfile+"_gbcut"+str(gallatcut)

if restselectioncut:
    print "Now applying N2Mass Cut"
    allarr = allarr.transpose()[allarr[4]>0].transpose()
    n2masscutno = len(allarr[0])
    print n2masscutno, ' objects survive n_2mass cut'
    print "Now applying w1mpro Cut"
    allarr = allarr.transpose()[(allarr[3]>12.0) * (allarr[3]<15.2)].transpose()
    w1mprocutno = len(allarr[0])
    print w1mprocutno, ' objects survive w1mpro cut'
    print "Now applying w1mpro - jm2mass Cut"
    allarr = allarr.transpose()[(allarr[3]-allarr[5])<-1.7].transpose()
    wmproj2masscutno = len(allarr[0])
    print wmproj2masscutno, ' objects survive w1mpro -jm2mass cut'
    print "Now applying jm2mass cuts"
    allarr = allarr.transpose()[allarr[5]<16.5].transpose()
    outfile = outfile+"_gselect"
    
if jmcut:
    allarr = allarr.transpose()[allarr[5]>15.8].transpose()
    j2masscutno = len(allarr[0])
    print j2masscutno, ' objects survive j2mass cut'
    outfile = outfile+"_jmcut"

if sgallatcut:
    print "Now applying supergalactic latitude cut"
    sgb = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).supergalactic.sgb.value
    allarr = allarr.transpose()[(sgb<-1.*sgallatcut) + (sgb>sgallatcut)].transpose()
    sgbcutno = len(allarr[0])
    print sgbcutno, 'objects survive Supergalactic latitude cut'
    outfile = outfile+"_sgbcut"+str(sgallatcut)
    
counttmrs=0
if tmrsovang:
    z = twomrsarr[2]/299792.458
    twomrsarr=twomrsarr.transpose()[z<zcut2mrs].transpose()
    print len(twomrsarr[0]), 'objects remaining in the 2MRS catalog after redshift cut' 
    inda,indb,c = k3match.celestial(allarr[0], allarr[1], twomrsarr[0],twomrsarr[1], tmrsovang)
    allarr = np.delete(allarr, inda, axis=1)
    outfile=outfile+'_2mrscut'+str(tmrsovang)+'_2mrszcut'+str(zcut2mrs)
    counttmrs=len(allarr[0])
    print counttmrs, 'objects survive 2mrs intersection cut'

countlrem = 0
if lrem:
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
    outfile = outfile+'_lrem'
  
countlvem = 0
if vel:
    totalmotion = np.rad2deg(np.power((np.power((np.cos(np.deg2rad(allarr[1]))*np.deg2rad(allarr[17]/3.6e6)), 2) + np.power(np.deg2rad(allarr[18]/3.6e6), 2)), 0.5))*3.6e6
    totalmotion[np.isnan(totalmotion)]=0
    allarr = allarr.transpose()[totalmotion<vel].transpose()
    countlvem = len(allarr[0])
    print countlvem, 'objects survive velocity cut'
    outfile=outfile+'_vcut'+str(vel)

countextc=0
if extc:
    allarr = allarr.transpose()[allarr[19]==0].transpose()
    countextc = len(allarr[0])
    print countextc, 'objects survive extension flag cut'
    outfile=outfile+'_extcut'
    
if options.SOUT:
    np.savetxt(outfile+'_galaxyselection.txt',allarr, delimiter="|")


finalno = len(allarr[0])

if not crwford:
    fout = open(outfile+'_result.txt', "w")
    fout.write(str(decmin)+ "," + str(decmax)+"\n")
    fout.write(str(totno)+","+str(gallatcutno)+","+str(n2masscutno)+","+str(w1mprocutno)+","+str(wmproj2masscutno)+","+str(j2masscutno)+","+str(sgbcutno)+","+str(counttmrs)+","+str(countlrem)+","+str(countlvem)+","+str(countextc)+","+str(finalno)+"\n")

    for i in range(0, npix):
        dec, ra = np.deg2rad((90. - np.rad2deg(hp.pix2ang(scanres, i)[0]))), hp.pix2ang(scanres, i)[1]
        print "DEC. ", (90. - np.rad2deg(hp.pix2ang(scanres, i)[0])), " R.A.", np.rad2deg(hp.pix2ang(scanres, i)[1])
        dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(allarr[1]))*np.cos(dec)*np.cos(np.deg2rad(allarr[0]) - ra)+np.sin(np.deg2rad(allarr[1]))*np.sin(dec)))
        print "UH: ", len(allarr.transpose()[dangle<90.]), "LH: ", len(allarr.transpose()[dangle>90.])
        fout.write(str(np.rad2deg(dec))+","+str(np.rad2deg(ra))+","+str(len(allarr.transpose()[dangle<90.]))+","+str(len(allarr.transpose()[dangle>90.]))+"\n")

else:
    fout = open(outfile+'_cwrfrd_result.txt', "w")
    fout.write(str(decmin)+ "," + str(decmax)+"\n")
    fout.write(str(totno)+","+str(gallatcutno)+","+str(n2masscutno)+","+str(w1mprocutno)+","+str(wmproj2masscutno)+","+str(j2masscutno)+","+str(sgbcutno)+","+str(counttmrs)+","+str(finalno)+"\n")
    skc = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).cartesian
    cx = skc.x.value
    cy = skc.y.value
    cz = skc.z.value
    xf = cx.sum()
    yf = cy.sum()
    zf = cz.sum()
    fout.write(str(xf)+','+str(yf)+','+str(zf)+"\n")
    
fout.close()
os.system('mv '+ outfile+'*_result.txt '+outfolder)
if options.SOUT:
    os.system('mv '+ outfile+'_galaxyselection.txt '+outfolder)
