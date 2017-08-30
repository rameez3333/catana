from numpy import *
from scipy import interpolate
from hmf import filters
from hmf import transfer
from hankel import SymmetricFourierTransform
from scipy.integrate import quad
from scipy.integrate import simps
from scalpy.fluids import LCDM
from hmf import growth_factor 
from scipy.special import erf, legendre
from collections import deque
from scipy.misc import derivative
import os, sys
from scipy.interpolate import InterpolatedUnivariateSpline
import itertools
import matplotlib.pyplot as plt

"""

This script produces interpolation tables for F(r, phi), by multipole expanding \delta as prescribed at 

http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node106.html#e13.33 and producing a numpy array of different F_n s to interpolate from

This is an intermediate commit. Not in final shape. delta_0 is hardcoded.


"""

NLEGENDRE = 25

rf = 20.

ASYMMETRYQ = 1.0

TSz0=transfer.Transfer(sigma_8=0.83,n=0.96,z=0,lnk_min=-35,lnk_max=35,transfer_model="BBKS")
GF=filters.Gaussian(TSz0.k,TSz0.power)
fft=SymmetricFourierTransform(ndim=3,N=int(1e5),h=1e-8)
Power=interpolate.UnivariateSpline(TSz0.k,TSz0.power,k=1,s=0) 


def Power2(k,Rf,AA,BB,sigma0):                                                           
	return Power(k)*exp(-(k*Rf)**2.0/2.0)**2.0*(AA+BB*k**2.0)/sigma0

def sigma0_Rstar_gamma(Rf):
    sigma0=GF.sigma(Rf,order=0)
    sigma1=GF.sigma(Rf,order=1)        
    sigma2=GF.sigma(Rf,order=2)
    gamma=sigma1**2.0/sigma0/sigma2
    Rstar=sqrt(3.0)*sigma1/sigma2
    return sigma0,Rstar,gamma


def delta(rr,Rf,delta0):
    sigma0,Rstar,gamma=sigma0_Rstar_gamma(Rf)
    nu=delta0/sigma0
    Theta=(3.0*(1-gamma**2.0)+(1.216-0.9*gamma**4.0)*exp(-(gamma/2.0)*(gamma*nu/2.0)**2.0))/((3.0*(1.0-gamma**2.0)+0.45+(gamma*nu/2.0)**2.0)**0.5+gamma*nu/2.0)

    AA=(nu-gamma**2.0*nu-gamma*Theta)/(1-gamma**2.0)
    BB=Theta*Rstar**2.0/(3.0*gamma*(1-gamma**2.0))
    Power22=lambda k: Power2(k,Rf,AA,BB,sigma0)
    Eta=fft.transform(Power22,rr,ret_err=False)/2.0/pi**2*sqrt(pi/2.0)/(2*pi)**1.5
    return Eta,sigma0




def PhiInterpConstructor():
    rs= power(10., linspace(1.e-2, log10(TSz0.cosmo.comoving_distance(1100).value), 10.))
    delta_0 = 0.222
    dlss=TSz0.cosmo.comoving_distance(1100)
    d=TSz0.cosmo.comoving_distance(0.52)
    phis = linspace(0, pi, 181.)
    rr=logspace(log10(0.01),log10(dlss.value),250)
    Delta,Sigma0=delta(rr/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)
    Delta=interpolate.UnivariateSpline(rr/TSz0.cosmo.h,Delta,k=1,s=0)
    
    def DeltaAs(rr, phi, q=2.6):
        return Delta(sqrt(power(rr*cos(phi)/ASYMMETRYQ,2.) + power(rr*sin(phi),2.)))
    
    args = asarray(list(itertools.product(rs, phis)))
    
    def DeltaIntF(phi, rr, n):
        print rr, phi, 'coords'
        return DeltaAs(rr, phi)*legendre(n)(cos(phi))*sin(phi)
    
    pPhiinte={}
    partPhi1={}
    partPhi2={}
    fullPhi={}
    PhineR = {}
    PhineRSpline = {}
    for i in range(NLEGENDRE):
        partPhi1[i] = lambda r : 0.5*r**(i+2)*quad(DeltaIntF, 0, pi, args=(r, i))[0]
        partPhi2[i] = lambda r : 0.5*r**(1-i)*quad(DeltaIntF, 0, pi, args=(r, i))[0]
        fullPhi[i] = lambda rr : rr**(-1.-i)*quad(partPhi1[i], 0, rr)[0] + rr**i*quad(partPhi2[i], rr, TSz0.cosmo.comoving_distance(1100).value)[0]
    
    for i in range(NLEGENDRE):
        print 'Step :', i
        PhineR[i] = asarray(map(fullPhi[i], rs))
        PhineRSpline[i] = InterpolatedUnivariateSpline(rs, PhineR[i], ext=3.)
    
    def Phirphi(coord):
        r, phi = coord[0], coord[1]
        print r, phi, 'coords'
        retval = 0
        for i in range(NLEGENDRE):
            retval = retval + PhineRSpline[i](r)*legendre(i)(cos(phi))
        return retval
    
    resarr = map(Phirphi, args)
    
    print resarr
    

def PhiConstructor():
    outl = []
    rs= power(10., linspace(-2., log10(TSz0.cosmo.comoving_distance(1100).value), 100.))
    outl.append(rs)
    print rs
    delta_0 = 0.222
    dlss=TSz0.cosmo.comoving_distance(1100)
    d=TSz0.cosmo.comoving_distance(0.52)
    phis = linspace(0, pi, 181.)
    rr=logspace(log10(0.01),log10(dlss.value),250)
    Delta,Sigma0=delta(rr/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)
    Delta=interpolate.UnivariateSpline(rr/TSz0.cosmo.h,Delta,k=1,s=0)
    
    def DeltaAs(rr, phi, q=2.6):
        return Delta(sqrt(power(rr*cos(phi)/ASYMMETRYQ,2.) + power(rr*sin(phi),2.)))

    deltan={}
    #deltanintegrater={}
    deltanlambdas={}
    deltaninterpolator = {}
    def DeltaNInteg(phi, r, n):
        return DeltaAs(r, phi)*legendre(n)(cos(phi))*sin(phi)

    plt.figure(0)
    for i in range(NLEGENDRE):
        print 'Step :', i
        deltanlambdas[i] = lambda r : (i+0.5)*quad(DeltaNInteg, 0, pi, args=(r,i))[0]
        deltan[i] = asarray(map(deltanlambdas[i], rs))
        deltaninterpolator[i] = InterpolatedUnivariateSpline(rs, deltan[i], ext=3.)
        rtry = power(10., linspace(1.e-2, log10(TSz0.cosmo.comoving_distance(1100).value), 10000.))
        plt.plot(rs, deltan[i],  ls='None', marker='*')
        plt.plot(rtry, deltaninterpolator[i](rtry), label=str(i))
        outl.append(deltan[i])
    
    plt.xscale('log')
    plt.legend(loc='best')
    plt.savefig('InterpEval.png')
    plt.show()
    return vstack(outl)






#PhiInterpConstructor()
arr = PhiConstructor()
savetxt('Deltac0.222_Rf20.txt', arr, delimiter='|')

#arr = genfromtxt('Deltac0.222_Rf20.txt', delimiter='|')

def DeltaApproximator(arr):
    
    def DeltaAs(rr, phi):
        interps={}
        retval = rr*0.
        for i in range(NLEGENDRE):
            interps[i] = InterpolatedUnivariateSpline(arr[0], arr[i+1], ext=3.)
            retval = retval + interps[i](rr)*legendre(i)(cos(phi))
        return retval
    return DeltaAs


def FACalculator(arr):
        interps={}
        integs={}
        phin={}
        phinterps={}
        outl=[]
        #retval = rr*0.
        
        rs= power(10., linspace(-2., log10(TSz0.cosmo.comoving_distance(1100).value), 1000.))
        outl.append(rs)
        for i in range(NLEGENDRE):
            interps[i] = InterpolatedUnivariateSpline(arr[0], arr[i+1], ext=3.)
        
        def integrand1(r, i):
            return r**(i+2)*interps[i](r)
        
        def integrand2(r, i):
            return r**(1-i)*interps[i](r)
        
        for i in range(NLEGENDRE):
            print 'Step ', i
            integs[i] = lambda r : quad(integrand1, 0, r, args=(i))[0]/((i+0.5)*r**(i+1)) + quad(integrand2, r, TSz0.cosmo.comoving_distance(1100).value, args=(i))[0]*r**(i)/(i+0.5)
            phin[i] = asarray(map(integs[i], rs))
            outl.append(phin[i])
        
        return vstack(outl)
        

phiarr  = FACalculator(arr)

savetxt('Phi_Rf'+str(rf)+'_q'+str(ASYMMETRYQ)+'_N'+str(NLEGENDRE)+'.txt', phiarr, delimiter='|')

def PhiApproximator(arr):
    
    def DeltaAs(rr, phi):
        interps={}
        retval = rr*0.
        for i in range(NLEGENDRE):
            interps[i] = InterpolatedUnivariateSpline(arr[0], arr[i+1], ext=3.)
            retval = retval + interps[i](rr)*legendre(i)(cos(phi))
        return retval
    return DeltaAs


PhiAs = PhiApproximator(phiarr)


rs= power(10., linspace(-2., log10(TSz0.cosmo.comoving_distance(1100).value), 9000.))

plt.figure(1)
plt.plot(rs, PhiAs(rs, 0.), label='phi = 0')
plt.plot(rs, PhiAs(rs, pi/3), label='phi = pi/3.')
plt.plot(rs, PhiAs(rs, pi/2), label='phi = pi/2.')
plt.plot(rs, PhiAs(rs, 5.*pi/6.), label='phi = 5.pi/6.')
plt.xscale('log')
plt.legend(loc='best')
plt.show()

#DeltaAsDipole = DeltaApproximator(arr)
    
    
#delta_0 = 0.222
#dlss=TSz0.cosmo.comoving_distance(1100)
#d=TSz0.cosmo.comoving_distance(0.52)
#phis = linspace(0, pi, 181.)
#rr=logspace(log10(0.01),log10(dlss.value),250)
#Delta,Sigma0=delta(rr/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)
#Delta=interpolate.UnivariateSpline(rr/TSz0.cosmo.h,Delta,k=1,s=0)

#def DeltaAs(rr, phi, q=2.6):
    #return Delta(sqrt(power(rr*cos(phi)/ASYMMETRYQ,2.) + power(rr*sin(phi),2.)))
    
    
#rs = power(10, random.uniform(1.e-2, log10(TSz0.cosmo.comoving_distance(1100).value), 10000))

#phis = random.uniform(0, pi, 10000)


#provals = DeltaAs(rs, phis)
#approxvals = DeltaAsDipole(rs, phis)

#err = approxvals/provals - 1.0

#plt.figure(1)
#plt.hist(err, bins=linspace(min(err), max(err), 3000))
#plt.show()
