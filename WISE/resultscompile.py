import numpy as np
from glob import glob
import healpy as hp
import matplotlib.pyplot as plt
import sys

nside = int(sys.argv[1])
#glcut = sys.argv[2]

npix = hp.nside2npix(nside)


def process(glcut, jcut):
        
        flist = sorted(glob('*'+str(glcut)+"glcut"+str(jcut)+"jlcut_result.txt"))
        if not len(flist):
            if jcut:
                flist = sorted(glob('*'+str(glcut)+"glcut"+"jcut_result.txt"))
            else:
                flist = sorted(glob('*'+str(glcut)+"glcut"+"_result.txt"))
        totsources=0
        cutsources=0

        totUH = np.zeros(npix)
        totLH = np.zeros(npix)
        for f in flist:
            fline = open(f).readlines()
            totsources = totsources + float(fline[1].split(",")[0])
            cutsources = cutsources + float(fline[1].split(",")[-1])
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
            if ((mapUH+mapLH)[0:npix].sum()/npix - float(fline[1].split(",")[-1])):
                print "wtf 3"
                print f
                print (mapUH + mapLH)[0:npix].sum()/(npix)
                print float(fline[1].split(",")[-1])
                print ((mapUH+mapLH)[0:npix/2].sum()*2./npix - float(fline[1].split(",")[-1]))

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


        hp.mollview(map, rot=[np.pi, np.pi/2.])
        #plt.show()
        plt.savefig('HCount'+str(glcut)+'glcut'+str(jcut)+'jcut.png')
        
        return (90. - np.rad2deg(hp.pix2ang(nside, np.argmax(map))[0])), np.rad2deg(hp.pix2ang(nside, np.argmax(map))[1]), map[np.argmax(map)]

for jcut in (0, 1):
    for glcut in np.arange(10, 21):
        decmax,ramax, valmax = process(glcut, jcut)
        print jcut, glcut, decmax, ramax, valmax
