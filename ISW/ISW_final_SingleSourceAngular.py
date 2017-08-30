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
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
import itertools

"""
This scipt adapts Xiaoyuans code to calculate the ISW temperature profile, but for a single void, elongated along the line of site. It relies on files produced by PotentialInterpolatormaker.py to interpolate for F(r, phi), since calculating it on the fly from the multipole expansion inside the integral seems to take impossibly long



"""



ASYMMETRYQ = 2.6
NLEGENDRE = 25


TSz0=transfer.Transfer(sigma_8=0.83,n=0.96,z=0,lnk_min=-35,lnk_max=35,transfer_model="BBKS")

def Power1(k,Rf):                                                           
	return Power(k)*exp(-(k*Rf)**2.0/2.0)**2.0 
def eta0(k,Rf):                                                             
	k=exp(k)                                    
	return k**3.0/(2*pi**2.0)*Power1(k,Rf)

def delta_random(rr,Rf,delta0):
	Power11=lambda k:Power1(k,Rf)
	Eta=fft.transform(Power11,rr,ret_err=False)/2.0/pi**2*sqrt(pi/2.0)/(2*pi)**1.5
	R=[0]+list(r)
	R=array(R)
	Eta=[quad(eta0,log(TSz0.k[0]),log(TSz0.k[-1]),args=Rf,limit=10000000)[0]]+list(Eta)
	Eta=array(Eta)
	phi=Eta/Eta[0]
	return R,phi*delta0

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
#	sigma0=GF.sigma(Rf,order=0)
#	sigma1=GF.sigma(Rf,order=1)        
#	sigma2=GF.sigma(Rf,order=2)
#	gamma=sigma1**2.0/sigma0/sigma2
#	Rstar=sqrt(3.0)*sigma1/sigma2
	sigma0,Rstar,gamma=sigma0_Rstar_gamma(Rf)
	nu=delta0/sigma0
	Theta=(3.0*(1-gamma**2.0)+(1.216-0.9*gamma**4.0)*exp(-(gamma/2.0)*(gamma*nu/2.0)**2.0))/((3.0*(1.0-gamma**2.0)+0.45+(gamma*nu/2.0)**2.0)**0.5+gamma*nu/2.0)

	AA=(nu-gamma**2.0*nu-gamma*Theta)/(1-gamma**2.0)
	BB=Theta*Rstar**2.0/(3.0*gamma*(1-gamma**2.0))
	Power22=lambda k: Power2(k,Rf,AA,BB,sigma0)
	Eta=fft.transform(Power22,rr,ret_err=False)/2.0/pi**2*sqrt(pi/2.0)/(2*pi)**1.5
	return Eta,sigma0


def w_critical(delta_0):
	Delta,Sigma0=delta(r/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)
	delta_sl=Delta[sum(diff(Delta)/diff(r)<0)]
	return (-b*delta_sl+1.0)/(-b*delta_0+1.0)


def dNddelta(delta0,Rf):
	sigma0,Rstar,gamma=sigma0_Rstar_gamma(Rf)
	nu=delta0/sigma0
#	print sigma0
	A=2.5/(9.0-5.0*gamma**2.0)
	B=432.0/((10.0*pi)**0.5*(9.0-5.0*gamma**2.0)**2.5)
	C1=1.84+1.13*(1.0-gamma**2.0)**5.72
	C2=8.91+1.27*exp(6.51*gamma**2.0)
	C3=2.58*exp(1.05*gamma**2.0)
	w=gamma*nu
	G=(w**3.0-3.0*gamma**2.0*w+(B*w**2.0+C1)*exp(-A*w**2.0))/(1.0+C2*exp(-C3*w))
	return 1/(2*pi)**2.0/Rstar**3.0*exp(-nu**2.0/2.0)*G/sigma0


def F(z,theta,rf,delta_0):
	d=TSz0.cosmo.comoving_distance(0.52)
	r=TSz0.cosmo.comoving_distance(z)
	rr=sqrt(d**2.0+r**2.0-2*d*r*cos(theta)).value
	F1=lambda rp: rp**2.0/rr*delta(rp/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)[0]
	F2=lambda rp: rp*delta(rp/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)[0]
	return quad(F1,0,rr)[0]+quad(F2,rr,TSz0.cosmo.comoving_distance(1100).value)[0] #Mpc^2

def F_test(z,theta,delta,d,dlss):
	r=TSz0.cosmo.comoving_distance(z)
	rr=sqrt(d**2.0+r**2.0-2*d*r*cos(theta)).value
	rp1=logspace(log10(0.1),log10(rr),100)                       ###careful, here rp with unit Mpc, but the corresponding delta is for Mpc/h????
	rp2=logspace(log10(rr),log10(dlss.value),100)
	delta1=delta(rp1)#/TSz0.cosmo.h)
	F1=simps(rp1**2.0/rr*(delta1*(abs(delta1)>1e-20)),rp1)
	delta2=delta(rp2)#/TSz0.cosmo.h)
	F2=simps(rp2*(delta2*(abs(delta2)>1e-20)),rp2)
#	F1=simps(rp1**2.0/rr*delta(rp1/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)[0],rp1)
#	F2=simps(rp2*delta(rp2/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)[0],rp2)
#	print rr,F1,F2

	return F1+F2 #Mpc^2	

def PhiApproximator(arr):
    
    def DeltaAs(rr, phi):
        interps={}
        retval = rr*0.
        for i in range(NLEGENDRE):
            interps[i] = InterpolatedUnivariateSpline(arr[0], arr[i+1], ext=3.)
            retval = retval + interps[i](rr)*legendre(i)(cos(phi))
        return retval
    
    rs= power(10., linspace(-2., log10(TSz0.cosmo.comoving_distance(1100).value), 100.))
    phis = linspace(0, pi, 181.)
    #cross = asarray(list(itertools.product(rs, phis))).transpose()
    pphis, rrs = meshgrid(phis, rs)
    Pot = DeltaAs(rrs, pphis)
    DeltaASap = RectBivariateSpline(rs, phis, Pot)
    return DeltaASap


arr1 = genfromtxt('Phi_Rf20.0_q1.0_N25.txt', delimiter='|')
PhiAS1 = PhiApproximator(arr1)

def FInterp_test1(z,theta,delta,d,dlss):
	r=TSz0.cosmo.comoving_distance(z)
	rr=sqrt(d**2.0+r**2.0-2*d*r*cos(theta)).value
	rp1=logspace(log10(0.1),log10(rr),100)                       ###careful, here rp with unit Mpc, but the corresponding delta is for Mpc/h????
	#phi1 = arcsin(d*sin(theta)/rp1)-theta
	rp2=logspace(log10(rr),log10(dlss.value),100)
	#phi2 = arcsin(d*sin(theta)/rp2)-theta
	#delta1=delta(rp1, phi1)#/TSz0.cosmo.h)
	#delta2=delta(rp2, phi2)#/TSz0.cosmo.h)
	
	phi = arcsin(d.value*sin(theta)/rr)-theta

	return PhiAS1(rr, phi)[0][0]

arr26 = genfromtxt('Phi_Rf20.0_q2.6_N25.txt', delimiter='|')
PhiAS26 = PhiApproximator(arr26)

def FInterp_test26(z,theta,delta,d,dlss):
	r=TSz0.cosmo.comoving_distance(z)
	rr=sqrt(d**2.0+r**2.0-2*d*r*cos(theta)).value
	rp1=logspace(log10(0.1),log10(rr),100)                       ###careful, here rp with unit Mpc, but the corresponding delta is for Mpc/h????
	#phi1 = arcsin(d*sin(theta)/rp1)-theta
	rp2=logspace(log10(rr),log10(dlss.value),100)
	#phi2 = arcsin(d*sin(theta)/rp2)-theta
	#delta1=delta(rp1, phi1)#/TSz0.cosmo.h)
	#delta2=delta(rp2, phi2)#/TSz0.cosmo.h)
	
	phi = arcsin(d.value*sin(theta)/rr)-theta

	return PhiAS26(rr, phi)[0][0]






def F2_test(z,theta,delta,d,dlss):
	r=TSz0.cosmo.comoving_distance(z)
	rr=sqrt(d**2.0+r**2.0-2*d*r*cos(theta)).value
	rp1=logspace(log10(0.1),log10(rr),100)                       ###careful, here rp with unit Mpc, but the corresponding delta is for Mpc/h????
	#phi1 = arcsin(d*sin(theta)/rp1)-theta
	rp2=logspace(log10(rr),log10(dlss.value),100)
	#phi2 = arcsin(d*sin(theta)/rp2)-theta
	#delta1=delta(rp1, phi1)#/TSz0.cosmo.h)
	#delta2=delta(rp2, phi2)#/TSz0.cosmo.h)
	
	phi = arcsin(d.value*sin(theta)/rr)-theta
	
	partPhi1={}
	partPhi2={}
	fullPhi1={}
	fullPhi2={}
	pPhiinte={}
	
	
	def DeltaIntF(phi, rr, n):
	    return delta(rr, phi)*legendre(n)(cos(phi))*sin(phi)
	
	for i in range(NLEGENDRE):
	    print 'Theta : ',rad2deg(theta), 'Step :', i
	    pPhiinte[i] = lambda r : 0.5*quad(DeltaIntF, 0, pi, args=(r, i))[0]
	    partPhi1[i] = asarray(map(pPhiinte[i], rp1))
	    partPhi2[i] = asarray(map(pPhiinte[i], rp2))
	    fullPhi1[i] = simps(rp1**(2.+i)/rr**(i+1)*(partPhi1[i]*(abs(partPhi1[i]))>1e-20),rp1)
	    fullPhi2[i] = simps(rp2**(1.-i)*rr**(i)*(partPhi2[i]*(abs(partPhi2[i]))>1e-20),rp2)
    
	retval = 0.
	for i in range(NLEGENDRE):
	    retval = retval + (fullPhi1[i] + fullPhi2[i])*legendre(i)(cos(phi))
	return retval #Mpc^2	

def dPhidt(z,theta,delta,d,dlss, F_test = F_test):
	z=exp(z)-1
#	G=(TSz0.cosmo.H(z)*(1.0-growth_rate(z))*x1.growth_factor(z)).value #km / (Mpc s)
	G=GG(z)
#	print z,G
#	return G*F_test(z,theta,delta,d,dlss)*TSz0.cosmo.inv_efunc(z)*(1.0+z)
	return 1.5*TSz0.cosmo.Om0*TSz0.cosmo.H0.value**2.0*G*F_test(z,theta,delta,d,dlss)*TSz0.cosmo.inv_efunc(z)*TSz0.cosmo.hubble_distance.value*(1.0+z) 
#	return 1.5*TSz0.cosmo.Om0*TSz0.cosmo.H0.value**2.0*G*F_test(z,theta,rf,delta_0)*TSz0.cosmo.hubble_distance.value*TSz0.cosmo.inv_efunc(z)#*(1.0+z) #(km/s)^3


def W_theta(theta):
	return theta<4.0 and 1 or -1

def growth_rate(z):
	return -derivative(x1.growth_factor,z,dx=1e-9)/x1.growth_factor(z)*(1+z)

def T(delta_0,rf):
	dlss=TSz0.cosmo.comoving_distance(1100)
	d=TSz0.cosmo.comoving_distance(0.52)
	rr=logspace(log10(0.01),log10(dlss.value),250)
	Delta,Sigma0=delta(rr/TSz0.cosmo.h,rf/TSz0.cosmo.h,delta_0)
	Delta=interpolate.UnivariateSpline(rr/TSz0.cosmo.h,Delta,k=1,s=0)
	T_theta=lambda theta:  W_theta(theta/pi*180.0)*2*pi*(theta)*quad(dPhidt,0,log(1100),args=(theta,Delta,d,dlss),epsrel=0.01)[0]
	A1=array(quad(T_theta,0,4.0/180*pi*0.999,epsrel=0.01))
	A2=array(quad(T_theta,4.0/180*pi*1.001,sqrt(2.0)*4.0/180*pi,epsrel=0.01))
	return (A1+A2)/(pi*(4.0/180.0*pi)**2.0)/(299792458.0/1e3)**3.0/1e-6*2.725*2.25*2.0
#	return quad(dPhidt,0,log(1100),args=(theta/180*pi,Delta,d,dlss),epsrel=0.01)[0]

def dNdx(x,delta0,Rf):
	x=exp(x)-1
	sigma0,Rstar,gamma=sigma0_Rstar_gamma(Rf)
	nu=delta0/sigma0
	A=2.5/(9.0-5.0*gamma**2.0)
	B=432.0/((10.0*pi)**0.5*(9.0-5.0*gamma**2.0)**2.5)
	C1=1.84+1.13*(1.0-gamma**2.0)**5.72
	C2=8.91+1.27*exp(6.51*gamma**2.0)
	C3=2.58*exp(1.05*gamma**2.0)
	w=gamma*nu
	G=(w**3.0-3.0*gamma**2.0*w+(B*w**2.0+C1)*exp(-A*w**2.0))/(1.0+C2*exp(-C3*w))
	F=exp(-(x-w)**2.0/2.0/(1-gamma**2.0))/(2*pi*(1-gamma**2.0))**0.5
	f=(x**3.0-3*x)*(erf((2.5)**0.5*x)+erf((2.5)**0.5*x/2.0))/2.0+(0.4/pi)**0.5*((31.0*x**2.0/4.0+1.4)*exp(-0.625*x**2.0)+(0.5*x**2.0-1.4)*exp(-2.5*x**2.0))
	return F*f/G*(1.0+x)


def dNdx_test(x,delta0,Rf):
	sigma0,Rstar,gamma=sigma0_Rstar_gamma(Rf)
	nu=delta0/sigma0
	A=2.5/(9.0-5.0*gamma**2.0)
	B=432.0/((10.0*pi)**0.5*(9.0-5.0*gamma**2.0)**2.5)
	C1=1.84+1.13*(1.0-gamma**2.0)**5.72
	C2=8.91+1.27*exp(6.51*gamma**2.0)
	C3=2.58*exp(1.05*gamma**2.0)
	w=gamma*nu
	G=(w**3.0-3.0*gamma**2.0*w+(B*w**2.0+C1)*exp(-A*w**2.0))/(1.0+C2*exp(-C3*w))
	F=exp(-(x-w)**2.0/2.0/(1-gamma**2.0))/(2*pi*(1-gamma**2.0))**0.5
	f=(x**3.0-3*x)*(erf((2.5)**0.5*x)+erf((2.5)**0.5*x/2.0))/2.0+(0.4/pi)**0.5*((31.0*x**2.0/4.0+1.4)*exp(-0.625*x**2.0)+(0.5*x**2.0-1.4)*exp(-2.5*x**2.0))
	return F*f/G

def delta1(rr,Rf):
	Power11=lambda k:Power1(k,Rf)
	Eta=fft.transform(Power11,rr,ret_err=False)/2.0/pi**2*sqrt(pi/2.0)/(2*pi)**1.5
	return Eta

def delta2(rr,Rf):
	sigma0,Rstar,gamma=sigma0_Rstar_gamma(Rf)
	AA=0
	BB=1
	Power22=lambda k: Power2(k,Rf,AA,BB,sigma0)
	Eta=fft.transform(Power22,rr,ret_err=False)/2.0/pi**2*sqrt(pi/2.0)/(2*pi)**1.5
	return Eta

def R_min(deltac,rf):

#	r=arange(1,100,0.5)*1.0
	r=arange(1,400,0.5)*1.0
	Ar=delta1(r/TSz0.cosmo.h,rf/TSz0.cosmo.h)
	Br=delta2(r/TSz0.cosmo.h,rf/TSz0.cosmo.h)

	sigma0,Rstar,gamma=sigma0_Rstar_gamma(rf/TSz0.cosmo.h)

	delta_min=deltac
	delta_max=1   
	delta=arange(deltac,1.0,1e-5)
	a=dNddelta(delta,rf/TSz0.cosmo.h) 
	dNddelta_min=a.min()
	dNddelta_max=a.max()

	n=1e5
	Delta=random.random(int(n))    
	DNddelta=random.random(int(n))
	delta=delta_min+(delta_max-delta_min)*Delta 
	dnddelta=dNddelta_min+(dNddelta_max-dNddelta_min)*DNddelta  
	delta=delta*(dnddelta<dNddelta(delta,rf/TSz0.cosmo.h))
	delta=delta[nonzero(delta)]
	n_sample=2000#delta.shape[0]/2
	ind=random.random_integers(0,delta.shape[0]-1,size=n_sample)
	RV=[]
	for i in ind:
		delta0=delta[i]

		x_min=0
		x_max=2.5   
		x_sample=arange(0,2.5,1e-5)
		a=dNdx(x_sample,delta0,rf/TSz0.cosmo.h) 
		dNdx_min=a.min()
		dNdx_max=a.max()


		n=1e5
		X=random.random(int(n))    
		DNDX=random.random(int(n))
		x_range=x_min+(x_max-x_min)*X 
		dNdx_range=dNdx_min+(dNdx_max-dNdx_min)*DNDX  
		x_range=x_range*(dNdx_range<dNdx(x_range,delta0,rf/TSz0.cosmo.h))
#x_range=x_range[nonzero(x_range)]
		x_range=exp(x_range[nonzero(x_range)])-1.0
		 
		A=(delta0-x_range*sigma0*gamma)/(1-gamma**2.0)/sigma0**2.0
		B=(gamma*delta0*Rstar**2.0-x_range*sigma0*Rstar**2.0)/(3*gamma*(1.0-gamma**2.0))/sigma0

		phi=Ar[None,:]*A[:,None]-Br[None,:]*B[:,None]
		Rv=r[argmin(phi,axis=1)] 
		RV.extend(list(Rv))
#		RV.extend(Rv)
	RV=array(RV)
	rv=percentile(RV,q=5.0)
	return rv,RV

def ISW(deltac,rf):
    
	Ta=lambda delta0:T(delta0,rf)[0]*dNddelta(delta0,rf/TSz0.cosmo.h)
	TISW=quad(Ta,deltac,1.0,epsabs=1,epsrel=0.01,limit=30,maxp1=30,limlst=30)[0]/quad(dNddelta,deltac,1.0,args=rf/TSz0.cosmo.h)[0]
	T2=lambda delta0:T(delta0,rf)[0]**2.0*dNddelta(delta0,rf/TSz0.cosmo.h)
	T2ISW=quad(T2,deltac,1.0,epsabs=1,epsrel=0.01,limit=30,maxp1=30,limlst=30)[0]/quad(dNddelta,deltac,1.0,args=rf/TSz0.cosmo.h)[0]
	DeltaT=sqrt(T2ISW-TISW**2.0)
	N1=quad(dNddelta,deltac,1.0,args=rf/TSz0.cosmo.h,limit=10000000)*(TSz0.cosmo.comoving_volume(0.75)-TSz0.cosmo.comoving_volume(0.4))

	N2=quad(dNddelta,deltac,1.0,args=rf/TSz0.cosmo.h,limit=10000000)[0]*5e9/TSz0.cosmo.h**3.0
	rv,RV=R_min(deltac,rf)
	return deltac,rf,TISW,DeltaT,rv,N2



Power=interpolate.UnivariateSpline(TSz0.k,TSz0.power,k=1,s=0) 
fft=SymmetricFourierTransform(ndim=3,N=int(1e5),h=1e-8)
GF=filters.Gaussian(TSz0.k,TSz0.power)

x1=growth_factor.GrowthFactor(TSz0.cosmo)  
x=LCDM(Om0=0.29,h=0.69)
growth_rate=x.growth_rate_z
	 

z=logspace(0,log10(1101),1000)-1.0
G=z*1
for i in range(len(z)):
	G[i]=(TSz0.cosmo.H(z[i])*(1.0-growth_rate(z[i]))*x1.growth_factor(z[i])).value
	
GG=interpolate.UnivariateSpline(z,G,k=1,s=0)

rf=20
deltac=0.222
#TISW,DeltaT,rv,N2=ISW(deltac,rf)
#print TISW,DeltaT,rv,N2
#result=[]
#for rf in arange(15,60,1):
    
    
    #res = ISW(deltac,rf)
    #print rf, res
    #result.append(res)

#save('ISW.dat',result)

dlss=TSz0.cosmo.comoving_distance(1100)
d=TSz0.cosmo.comoving_distance(0.52)
rr=logspace(log10(0.01),log10(dlss.value),250)
Delta,Sigma0=delta(rr/TSz0.cosmo.h,rf/TSz0.cosmo.h,deltac)
Delta=interpolate.UnivariateSpline(rr/TSz0.cosmo.h,Delta,k=1,s=0)
T_theta=lambda theta: quad(dPhidt,0,log(1100),args=(theta,Delta,d,dlss),epsrel=0.01)[0]/(299792458.0/1e3)**3.0/1e-6*2.725*2.25*2.0


def DeltaAs(rr, phi, q=2.6):
    return Delta(sqrt(power(rr*cos(phi)/ASYMMETRYQ,2.) + power(rr*sin(phi),2.)))

T_thetaAs=lambda theta:  quad(dPhidt,0,log(1100),args=(theta,DeltaAs,d,dlss,F2_test),epsrel=0.01)[0]/(299792458.0/1e3)**3.0/1e-6*2.725*2.25*2.0


thetas = linspace(0, 20, 100)


Toftheta = asarray(map(T_theta, deg2rad(thetas)))
print 'Doing the Asymmetric bit now'




T_thetaAs1=lambda theta:  quad(dPhidt,0,log(1100),args=(theta,DeltaAs,d,dlss,FInterp_test1),epsrel=0.01)[0]/(299792458.0/1e3)**3.0/1e-6*2.725*2.25*2.0

T_thetaAs26=lambda theta:  quad(dPhidt,0,log(1100),args=(theta,DeltaAs,d,dlss,FInterp_test26),epsrel=0.01)[0]/(299792458.0/1e3)**3.0/1e-6*2.725*2.25*2.0

#The 0.5 below is to account for a factor of 2 in going from Eq 2.7 in Nadathur et al 2012 to the multipole expansions

TofthetaAs26 = 0.5*asarray(map(T_thetaAs26, deg2rad(thetas)))
TofthetaAs1 = 0.5*asarray(map(T_thetaAs1, deg2rad(thetas)))

#Rescaling the temperature profile of the elongated void so that far away from the center, it is the same as that of the symmetric void. This is necessary because the asymmetric density profile (see def DeltaAs above) has been modified in such a way that the considering one void in isolation, the universe now is no longer isotropic. 

TofthetaAs26 = TofthetaAs26*TofthetaAs1[-1]/TofthetaAs26[-1]

print TofthetaAs1/Toftheta


plt.plot(thetas, Toftheta-Toftheta[-1], label = 'Eq 2.7 of Nadathur et al', ls='--', lw=3)
plt.plot(thetas, TofthetaAs1-TofthetaAs1[-1], label = 'Q=1')
plt.plot(thetas, TofthetaAs26-TofthetaAs26[-1], label = 'Q=2.6')
plt.legend(loc='best')
plt.ylabel(r'$\Delta T[\mu K]$', fontsize=15)
plt.xlabel(r'$\theta [\degree]$', fontsize=15)
plt.title(r'$\mathrm{Single Void} \delta_0 = 0.222, R_f = 20 h^{-1}Mpc, z=0.52$', fontsize=15)
#plt.ylabel()
#plt.yscale('log')
plt.savefig('DeltaTvsThetaSingleVoid.png')
plt.show()
