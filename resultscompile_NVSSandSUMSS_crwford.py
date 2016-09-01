import numpy as np
from glob import glob
import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import sys

#nside = int(sys.argv[1])
##glcut = sys.argv[2]

#npix = hp.nside2npix(nside)


def process(tagname):
        #print "NVSS/"+tagname+"_result.txt"
        flist = sorted(glob("NVSS/"+tagname+"_result?*.txt"))
        print flist
	print len(flist), ' files found'
        totsources=0
        cutsources=0
        tmrsremsources=0
        totUH = np.zeros(npix)
        totLH = np.zeros(npix)
        for f in flist:
            fline = open(f).readlines()
            totsources = totsources + float(fline[1].split(",")[0])
            cutsources = cutsources + float(fline[1].split(",")[-1])
            tmrsremsources = tmrsremsources + (float(fline[1].split(",")[-3]) - float(fline[1].split(",")[-2]))
	    if tmrsremsources < 0:
		    tmrsremsources = tmrsremsources + (float(fline[1].split(",")[-3]) - float(fline[1].split(",")[-2]))
            i=0
            mapUH = np.zeros(npix)
            mapLH = np.zeros(npix)
            for line in fline[2:]:
                #if (90. - np.rad2deg(hp.pix2ang(16, i)[0]) - float(line.split(',')[0]) + float(line.split(',')[1]) - np.rad2deg(hp.pix2ang(16, i)[1])):
                #print "wtf", i, (90. - np.rad2deg(hp.pix2ang(16, i)[0]) - float(line.split(',')[0]) + float(line.split(',')[1]) - np.rad2deg(hp.pix2ang(16, i)[1]))
                
                mapUH[i] = float(line.split(',')[2])
                mapLH[i] = float(line.split(',')[3])
                i+=1
            if (i - npix):
                print "wtf now"
#            if ((mapUH+mapLH)[0:npix].sum()/npix - float(fline[1].split(",")[-1])):
                #print "wtf 3"
                #print f
                #print (mapUH + mapLH)[0:npix].sum()/(npix)
                #print float(fline[1].split(",")[-1])
                #print ((mapUH+mapLH)[0:npix/2].sum()*2./npix - float(fline[1].split(",")[-1]))

            totUH = totUH+mapUH
            totLH = totLH+mapLH

        """
        for i in range(1, npix/2 +1):
            totUH[npix-i] = totLH[i-1]
            totLH[npix-i] = totUH[i-1]
        """
        #print totUH
        #print totLH

        map = (totUH-totLH)/(totUH+totLH)

        #print map

        #print np.min(map), np.max(map)

        #print "Total Sources: ", totsources
        #print "After Cut: ", cutsources

        #print "The minimum is at", np.argmin(map), "with a value of ",  map[np.argmin(map)]

        #print "The minimum is at ", (90. - np.rad2deg(hp.pix2ang(nside, np.argmin(map))[0])), np.rad2deg(hp.pix2ang(nside, np.argmin(map))[0])

        #print "The maximum is at", np.argmax(map), "with a value of ",  map[np.argmax(map)]

        #print "The maximum is at ", (90. - np.rad2deg(hp.pix2ang(nside, np.argmax(map))[0])), np.rad2deg(hp.pix2ang(nside, np.argmax(map))[0])

        coords = hp.pix2ang(nside,np.arange(npix))

        angs = np.rad2deg(np.arccos(np.sin(coords[0])*np.sin(hp.pix2ang(nside, np.argmax(map))[0])*np.cos(coords[1] - hp.pix2ang(nside, np.argmax(map))[1]) + np.cos(coords[0])*np.cos(hp.pix2ang(nside, np.argmax(map))[0])))

        #print np.size(angs)

        #plt.scatter(angs, map)
        #plt.xlabel("Angle")
        #plt.ylabel("Hemispheric count difference")
        if tmrsremsources<0.:
            tmrsremsources=np.nan

        hp.mollview(map)
        #plt.show()
        plt.savefig('HCount'+tagname+'.png')
        
        return (90. - np.rad2deg(hp.pix2ang(nside, np.argmax(map))[0])), np.rad2deg(hp.pix2ang(nside, np.argmax(map))[1]), map[np.argmax(map)], cutsources, tmrsremsources

fout = open("NVSSandSUMSSasciiresultsummaryno2MRSxcorrectedcrwford.txt", "w")
fout.write("Angular tolerance within which 2MRS sources are looked for| Sourced removed because of 2MRS overlap|Supergalactic latitude cut(+/-)|Galactic Latitude Cut(+/-)|NVSS-SUMSS patch declination|Flux threshold (mJy)|Total Number of sources left after cuts|dipole DEC(deg)|dipole RA(deg)|dipole value|Galactic b (deg)|Galactic(l)(deg)|Velocity(km/s)\n")

counter =0
for tmrscut in [0.00028, 0.0028, 0.01]:
    for gbcut in [10, 15, 20]:
        for sgcut in [0,5,10]:
            for tpatch in [-40, -30]:
                for tcut in [10,15,20,25, 30, 35, 40, 45, 50]:
                    outfile = 'nvsumss_patch'+str(float(tpatch))+'ascii'
                    outfile = outfile+"_gbcut"+str(gbcut)
                    outfile = outfile+"_fluxcut"+str(tcut)
                    if sgcut:
                        outfile = outfile+"_sgbcut"+str(sgcut)
                    outfile = outfile+'_2mrscut'+str(tmrscut)
                    #decmax, ramax, valmax, nsources, tmrssources = process(outfile)
                    flist = glob("NVSS/"+outfile+"_cwrfrd_result.txt")
                    print len(flist), "files found"
                    lines = open(flist[0]).readlines()
                    nsources = float(lines[1].split(",")[-1])
                    tmrssources = float(lines[1].split(",")[-3]) - float(lines[1].split(",")[-2])
                    ramax, decmax, valmax = float(lines[2].split(",")[0]) , float(lines[2].split(",")[1]), float(lines[2].split(",")[2])
                    
                    print sgcut, gbcut, tcut, nsources, decmax, ramax, valmax
                    gb = SkyCoord(ra = ramax*u.degree, dec=decmax*u.degree).galactic.b.value
                    gl = SkyCoord(ra = ramax*u.degree, dec=decmax*u.degree).galactic.l.value
                    v = float(lines[2].split(",")[3])
                    fout.write(str(tmrscut)+"|"+str(tmrssources)+"|"+str(sgcut)+"|"+str(gbcut)+"|"+str(float(tpatch))+"|"+str(tcut)+"|"+str(nsources)+"|"+str(decmax)+"|"+str(ramax)+"|"+str(valmax)+"|"+str(gb)+"|"+str(gl)+"|"+str(v)+"\n")
                    counter+=1

print "Total ", counter
            


            
