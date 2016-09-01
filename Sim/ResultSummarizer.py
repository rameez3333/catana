import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from ROOT import TH1D, TLegend



color = {}
color[10000] = 'green'
color[100000] = 'blue'
color[1000000] = 'red'
color[600000] = 'brown'
color[2500000] = 'cyan'
ncats = [100000, 1000000]

ncatc = [100000, 1000000, 2500000]

ncat=10000
flist = glob('*Sim_125.0_-30.0_*_'+str(ncat)+'_summary.txt')
ras = []
decs = []
vels = []
for f in flist:
    lin = open(f).readlines()[4]
    ras.append(float(lin.split("|")[3]))
    decs.append(float(lin.split("|")[2]))
    vels.append(float(lin.split("|")[4]))
ras = np.asarray(ras)
decs = np.asarray(decs)

def makevhist(ncat, c=False):
    if not c:
        flist = glob('*Sim_125.0_-30.0_*_'+str(ncat)+'_summary.txt')
    if c:
        flist = glob('*Simcrwfrd_168.0_-7.0_*_'+str(ncat)+'_summary.txt')
    print len(flist), 'files found'
    ras = []
    decs = []
    vels = []
    for f in flist[0:400]:
        try: 
            lin = open(f).readlines()[4]
            ras.append(float(lin.split("|")[3]))
            decs.append(float(lin.split("|")[2]))
            vels.append(float(lin.split("|")[4])*2./(2.+1.5*(1+0.75)))
        except:
            print 'Weird:', f
    ras = np.asarray(ras)
    decs = np.asarray(decs)
    vels = np.asarray(vels)
    velhist = TH1D("hist"+str(ncat), 'Catalog size '+str(ncat), len(flist)/10,np.min(vels),np.max(vels))
    velhist.FillN(len(vels), vels, np.ones(len(vels)))
    return vels

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
    ax.scatter(np.radians(x),np.radians(Dec), alpha=0.5, s=5)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
    
plot_mwd(ras, decs, title="Random dipole directions on Catalogsize "+str(ncat))
plt.show()
#plt.show()
#plt.savefig('randdirectionsn'+str(ncat)+'.png')

histdict={}
histcdict={}
legend = TLegend(0.1,0.7,0.48,0.9)


for ncat in ncats:
    histdict[ncat] = makevhist(ncat)
for ncat in ncatc:
    histcdict[ncat] = makevhist(ncat, True)
    #legend.AddEntry(histdict[ncat], "Catalog size"+str(ncat))


for ncat in ncats:
    plt.hist(histdict[ncat], 20, normed=1 ,facecolor=color[ncat], label='HemiSphere : catalog size '+str(ncat), alpha=0.5)
for ncat in ncatc:
    plt.hist(histcdict[ncat], 20, normed=1 ,facecolor=color[ncat], label='Crawford : catalog size '+str(ncat), alpha=0.3)

#histdict[ncat].Draw()
#for ncat in ncats:
    #histdict[ncat].Draw("Same")
    
#legend.Draw("Same")

plt.xscale('log')
plt.xlabel('velocity(km/s)')
plt.legend(loc='best', fontsize = 16, ncol=2)
plt.show()
plt.savefig('VelocityDist.png')

#raw_input('Whatever')
