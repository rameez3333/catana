import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from glob import glob

flist = sorted(glob("AllWISE/wise-allwise-cat-part??_gbcut15_gselect_sgbcut5_2mrscut0.00028_2mrszcut0.02_lrem_galaxyselection.txt"))

arr = np.genfromtxt(flist[0], delimiter="|")

for fname in flist[1:]:
    print fname
    arr = np.append(arr, np.genfromtxt(fname, delimiter="|"), axis=1)
    

arr[9][np.isnan(arr[9])] = np.nanmax(arr[9])

flux1 = np.power(10., -1.*arr[3])
flux2 = np.power(10., -1.*arr[9])

print flux1

print flux2

flux=flux1

N=[]

S = np.power(np.linspace(np.log10(np.min(flux)), np.log10(np.max(flux)), 1000),10)


for i in S:
    N.append(len(flux[flux>i]))
    
print S
print N
    
def func(x, k, i):
    return k*np.power(x, -1*i)

popt, pcov = curve_fit(func, S, N)

x = np.linspace(np.min(flux), np.max(flux), 2000)

y = func(x, popt[0], popt[1])

plt.scatter(S,N, s=10.)

plt.plot(x,y, linewidth=2, color='red')

plt.xlabel('Flux Threshold value', fontsize=20.)
plt.ylabel('Sources above threshold', fontsize=20.)

#plt.loglog()
plt.show()
