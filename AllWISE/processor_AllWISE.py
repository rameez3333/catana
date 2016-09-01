import numpy as np
import healpy as hp
import sys
import astropy import units as u
from astropy.coordinates import SkyCoord


scanres = int(sys.argv[1])
scanfile = sys.argv[2]
gallatcut = int(sys.argv[3])
jlcutyn = int(sys.argv[4])
sgallatcut  = int(sys.argv[5])

npix = hp.nside2npix(scanres)

print "Loading file"

def nulltonan(val):
    if val=='NULL' or val=='null':
        return np.nan
    else:
        return float(val)

def SpaceAngleRad(zen1, azi1, zen2, azi2):
    return acos( sin(zen1)*sin(zen2)*cos(azi1-azi2)+cos(zen1)*cos(zen2))


allarr = np.genfromtxt(scanfile, usecols=(1,2,7,16,286,287 ), delimiter='|', converters={7:nulltonan,16:nulltonan,286:nulltonan,287:nulltonan}).transpose()



decmin, decmax = np.min(allarr[1]), np.max(allarr[1])
totno = len(allarr[0])


print totno, ' objects in declination range ', decmin, ' - ', decmax, ' loaded'


print "Applying Gal Lat cut, +/-", gallatcut 

allarr = allarr.transpose()[(allarr[2]<-1.*gallatcut) + (allarr[2]>gallatcut)].transpose()

gallatcutno = len(allarr[0])

print gallatcutno, ' objects survive Galactic Latitude cut'

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
#if jlcutyn:
    #allarr = allarr.transpose()[allarr[5]>15.8].transpose()


j2masscutno = len(allarr[0])

print j2masscutno, ' objects survive j2mass cut'

print "Now applying supergalactic latitude cut"

sgb = SkyCoord(ra = allarr[0]*u.degree, dec=allarr[1]*u.degree).supergalactic.sgb.value

allarr = allarr.transpose()[(sgb<-1.*sgallatcut) + (sgb>sgallatcut)].transpose()

sgbcutno = len(allarr[0])

print sgbcutno, 'objects survive Supergalactic latitude cut'


np.savetxt(scanfile+'_galaxyselection.txt',allarr, delimiter="|")
allarr = allarr.transpose()[allarr[5]>15.8].transpose()
np.savetxt(scanfile+'_jcut_galaxyselection.txt',allarr, delimiter="|")

#fout = open(scanfile+"_"+str(scanres)+'nside_'+ str(gallatcut) +'glcut'+str(jlcutyn)+'jmcut_result.txt', "w")

#fout.write(str(decmin)+ "," + str(decmax)+"\n")
#fout.write(str(totno)+","+str(gallatcutno)+","+str(n2masscutno)+","+str(w1mprocutno)+","+str(wmproj2masscutno)+","+str(j2masscutno)+"\n")

#for i in range(0, npix):
    #dec, ra = np.deg2rad((90. - np.rad2deg(hp.pix2ang(scanres, i)[0]))), hp.pix2ang(scanres, i)[1]
    #print "DEC. ", (90. - np.rad2deg(hp.pix2ang(scanres, i)[0])), " R.A.", np.rad2deg(hp.pix2ang(scanres, i)[1])
    #dangle = np.rad2deg(np.arccos(np.cos(np.deg2rad(allarr[1]))*np.cos(dec)*np.cos(np.deg2rad(allarr[0]) - ra)+np.sin(np.deg2rad(allarr[1]))*np.sin(dec)))
    #print "UH: ", len(allarr.transpose()[dangle<90.]), "LH: ", len(allarr.transpose()[dangle>90.])
    #fout.write(str(np.rad2deg(dec))+","+str(np.rad2deg(ra))+","+str(len(allarr.transpose()[dangle<90.]))+","+str(len(allarr.transpose()[dangle>90.]))+"\n")

#fout.close()
