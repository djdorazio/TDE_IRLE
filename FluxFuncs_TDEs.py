import numpy as np
import math as ma
#import numexpr as ne
import scipy as sc

#from scipy.optimize import brentq #fmin

import scipy.integrate as intg
import scipy.signal as sgn
#from apw import Tfunc

#Goal is to output the IR flux from irradiated dust geometry surrounding TDE.


# ###FOR TRAP INT
# Ntrap_ph = 10.
# Ntrap_th = 10.
# Ntrap_nu = 10.

#### INTEGRATION ERROR TOLS
myrel = 1.e-8
myabs = 1.e-8#1.e-60
reclim = 1#100
limlst = 1
maxp1 = 1
fo = 1

##GLOBAL PHYSICS CONSTANTS (cgs):
c = 2.9979*10**(10)
G = 6.673*10**(-8)
Msun = 1.998*10**(33)

kb = 1.3807*10**(-16)
h = 6.62607*10**(-27) 
sigSB = 5.670*10**(-5)

pc2cm = 3.08567758*10**(18)
numicron = c/(10**(-4))

yr2sec = 3600.*24.*365.25
days2yrs = 1./365.



## Back body specific intensity
def Bv(nu, T):
	return 2.*h*nu*nu*nu/(c*c)*1./(np.exp( h*nu/(kb*T) ) - 1.)

## Dust absorption/emission efficiency
def Qv(nu, nu0, nn):
	qv = (nu/nu0)**(nn)
	qv = np.minimum(qv, 1.0*nu/nu)
	return qv


def QvBv(nu, T, nu0, nn):
	qv = (nu/nu0)**(nn)
	qv = np.minimum(qv, 1.0*nu/nu)
	return 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*T    )  ) - 1.) * qv
	       


####################################################
### Dust number densityprofile - IGNORE BELOW FOR NOW
####################################################
def nDust(x,y,z, n0, Rd, p, thetT, JJ):
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	r  = (x*x + y*y + z*z)**(0.5)

	throt = np.arctan2( (xrot*xrot + y*y)**(0.5), zrot)   ##arctan of arg1/arg2 arg1 always positive so btwn 0, pi
	throt = np.array(throt)

	nprof = 0.0			
	if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
		nprof = n0*(Rd/r)**(p)	

	return nprof





def nprof(r, n0, Rd, p):
	return n0*(Rd/r)**(p)

def nDust_pcwse(x,y,z, n0, Rd, p, thetT, JJ):
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	r  = (x*x + y*y + z*z)**(0.5)

	throt = np.arctan2( (xrot*xrot + y*y)**(0.5), zrot)   ##arctan of arg1/arg2 arg1 always positive so btwn 0, pi


	#nprof = n0*(Rd/r)**(p)
	#BoxCar = Heaviside[throt - thetT] - Heaviside[throt - (ma.pi - thetT)]
	#return nprof * Heaviside[r-Rd] * BoxCar

	#nd = np.piecewise(r, [r < Rd, r >= Rd], [lambda r:0.0, lambda r:n0*(Rd/r)**(p)])
	nd = np.piecewise(throt, [throt>=(np.pi - thetT), throt<=(np.pi - thetT)], [lambda throt:0.0, lambda throt:n0])
	nd = np.piecewise(throt, [throt<thetT, throt>=thetT], [lambda throt:0.0, lambda throt:nd])

	#or throt>=(np.pi - thetT)
	#and throt<=(np.pi - thetT)
	return nd
####################################################
### Dust number densityprofile - IGNORE ABOVE FOR NOW
####################################################







## equation to tabulate RHS and T
def T_RHS(Td, nu0, nn):
	# 4 for difference in cross sectional area and surface area, pi for isotropic flux from Grain
	RHS = 4.*ma.pi* intg.quad(QvBv  ,0., 3.0*numicron, args=(Td, nu0, nn), epsabs=0.0 )[0] 
	return RHS












### optical depth - IGNORE FOR NOW
def tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn):
	xe     = (Rout*Rout - (z*z + y*y))**(0.5)
	return ma.pi*aeff*aeff * Qv(nu, nu0, nn) * intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	



















####################################################
### optical/UV source 
####################################################
def Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac):
	F0 = Lavg/(4.*ma.pi*r*r)
	FQ = F0/100.#FQfac
	gam = 2.0


	if ((t-t0)<=0.0):
		return FQ
	else:
		Fsrc = F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam)
		if (Fsrc<=0.0):
			Fsrc=FQ
		return Fsrc

#Fsrc_Anl = np.vectorize(Fsrc_Anl)


def Fsrc_Anl_Fit(t, r, Lavg, tfb, t0, gam, FQfac):
	F0 = Lavg/(4.*ma.pi*r*r)
	FQ = F0/FQfac
	#nuVbnd = c/(5.45*10**(-5))
	#FVbndRel = 3.636*10**(-20)*nuVbnd 
	#BC=10.0
	
	if (type(t) is float or type(t) is np.float64):
		if (t-t0<=0.0):
			return FQ 
		else:
			Fsrc = F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam)
			if (Fsrc<=0.0):
				Fsrc=FQ
			return Fsrc 

	else:
		res = []
		for i in range (len(t)):
			if (t[i]-t0<=0.0):
				Fsrc = FQ
			else:
				Fsrc = F0*(1.+np.log( (t[i]-t0)/tfb) ) * ((t[i]-t0)/tfb)**(-gam)
				if (Fsrc<=0.0):
					Fsrc=FQ

			res.append(Fsrc)
		return np.array(res)



####################################################
### Compute Dust temperature from Therm Eql
####################################################
def TDust(t,r,thet, ph, args, RHStable, Ttable):
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR = args
	#Rd = 0.0
	#Loft = Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac) * 4.*ma.pi*r*r
	#r = etaR * 0.5 * pc2cm * (Loft/(10.**(46)))**(0.5) ##sublimatin front

	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	
	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	Td = 1.0*r/r  ##T=1 is very small

	if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###
		Fsrc = Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac)
		


		###-----------------###
		### Compute taudust ###
		###-----------------###
		Qbar=1. 
		#if r=Rd, tau is just 0
		tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
		### if flux is greater than RHS max at which T > Tsub~1800K, then dust sublimates
		LHS = Fsrc #* np.exp(-tauDust)
		# if (type(Fsrc) is np.ndarray):
		# 	if (len(Fsrc) > len(RHStable)):
		# 		LHS = sgn.resample(LHS, len(RHStable))
		# 	elif (len(Fsrc) < len(RHStable)):
		# 		RHStable = sgn.resample(RHStable, len(Fsrc))
		RHS_mx = RHStable[len(RHStable)-1]
		RHS_mn = RHStable[0]

		# if (type(t) is float or type(t) is np.float64):
		# 	if (LHS > RHS_mx or LHS <= RHS_mn):
		# 		Td = 1.
		# 	else:
		# 		istar = np.where( LHS <= RHStable )[0].min()
		# 		Td = Ttable[istar]
		# 	return Td
		# else:
		# 	res = []
		# 	for i in range(len(LHS)):
		# 		if (LHS[i] > RHS_mx or LHS[i] <= RHS_mn):
		# 			Td = 1.
		# 		else:
		# 			istar = np.where( LHS[i] <= RHStable )[0].min()
		# 			Td = Ttable[istar]

		# 		res.append(Td)

		# 	return res

		if (LHS > RHS_mx or LHS <= RHS_mn):
			Td = 1.
		else:
			istar = np.where( LHS <= RHStable )[0].min()
			Td = Ttable[istar]

	return Td

#TDust = np.vectorize(TDust, excluded=[4,5,6])

####################################################
### T is analytic when Q_nu=1
####################################################
def TDust_Anl(t,r,thet, ph, args):
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR = args

	#Loft = Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac) * 4.*ma.pi*r*r
	#r = etaR * 0.5 * pc2cm * (Loft/(10.**(46)))**(0.5) ##sublimatin front

	
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	Td = 1.*r/r  ##T=1 is very small

	if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###

		Fsrc = Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac)
		#Fsrc = Fsrc_Anl_subl(t, r, Lavg, tfb, t0, FQfac)


		###-----------------###
		### Compute taudust ###
		###-----------------###
		Qbar=1. 
		tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
		### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
		LHS = Fsrc #* np.exp(-tauDust)
		Td = (LHS/(4.*sigSB))**(1./4.)

	return Td





####################################################
### This is the spherical shell case resulting  
### from the assumption that dust is Opt Thick 
### to optical/UV and Opt Thin to IR
####################################################

def Fnuint_UVthick_IRThin_Iso(ph, thet, nu, t, Dist, args, RHStable, Ttable):
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR = args

	Loft = Fsrc_Anl(t, Rd, Lavg, tfb, t0, FQfac) * 4.*ma.pi*Rd*Rd
	Rd = etaR * 0.5 * pc2cm * (Loft/(10.**(46)))**(0.5)
###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted from dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))



	# Tdust 
	Tdust = TDust(tem, Rd, thet, ph,  args, RHStable, Ttable)
	#Tdust = TDust_Anl(tem, Rd, thet, ph, args)


	if (Tdust==1.0):
		fint=1.e-16
	else:
		# surface density in optically thick limit
		Surf_nd = 1./(ma.pi*aeff*aeff)
		# Rd is the inner edge of the shell
		fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
		fint = fint* Rd*Rd* np.sin(thet) * Surf_nd #mult by spherical Jacobian and surface denisty


	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist * fint


#Fnuint_UVthick_IRThin_Iso = np.vectorize(Fnuint_UVthick_IRThin_Iso, excluded=[5,6,7])


#################################################################
### For dust which is optically thin to opt/UV ##IGNORE FOR NOW
#################################################################
## Allow the optical/UV to penetrate the dust, so optically think to opical/UV out to some distance
## must still assume opt thin to IR
def Fnuint_OptThin_IRThin_Iso(ph, thet, r, nu, t, Dist, args, RHStable, Ttable): #, tauGrid):	
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR = args
	## where UV penetrates out to (tau_UV=1)
	Rout = 2. * Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) #+ 0.01*Rd



###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
	## retarded time - time light emitted form dust
	tem = t - r/c*(1. - np.sin(thet)*np.cos(ph))



	# Tdust for doppler source
	Tdust = TDust(tem, Rd, thet, ph,  args, RHStable, Ttable)
	#Tdust = TDust_Anl(tem, Rd, thet, ph, args)

	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)
	fint = fint* r*r* np.sin(thet) * nDust(x,y,z, n0, Rd, p, thetT, JJ)


	return ma.pi* aeff*aeff/Dist/Dist *fint





























































####################################################
###   INTEGRATION   
####################################################





####################################################
###   GAUSSIAN QUADRATURE
####################################################
##################
### F_ISO
##################
#########
###Torus Shell - OPT ThICK to optical/UV, OPT THIN to IR
#########
def FThnu_UVthick_IRThin_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_UVthick_IRThin_Iso, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]

def Fnu_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThnu_UVthick_IRThin_QuadInt, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

def F_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	if (type(t) is float or type(t) is np.float64):
		return intg.quad(Fnu_UVthick_IRThin_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]
	else:
		res = []
		for i in range(len(t)):
			res.append(intg.quad(Fnu_UVthick_IRThin_QuadInt, numin, numax, args=(t[i], Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0])
		return np.array(res)


# ####USE TRAP RULE
# def FThnu_UVthick_IRThin_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable):
# 	phis = np.linspace(0.0,2.*ma.pi, Ntrap_ph)
# 	return intg.trapz( Fnuint_UVthick_IRThin_Iso(phis, thet, nu, t, Dist, Aargs, RHStable, Ttable) )
# 	# #phis = np.arange(0.0, 2.*ma.pi, 0.1)
# 	# arg = []#np.zeros(len(phis))
# 	# for i in range(len(phis)):
# 	# 	#arg.append(Fnuint_UVthick_IRThin_Iso(phis[i], 0.0, numicron, tt, Dst, argW1, RHS_table, T_table))
# 	# 	arg.append(Fnuint_UVthick_IRThin_Iso(phis[i], thet, nu, t, Dist, Aargs, RHStable, Ttable))
# 	# 	#ans = intg.trapz(np.array(arg))

# 	# return intg.trapz( np.array(arg) )
# FThnu_UVthick_IRThin_QuadInt = np.vectorize(FThnu_UVthick_IRThin_QuadInt, excluded=[4,5,6])


# def Fnu_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
# 	ths = np.linspace(0.0, ma.pi, Ntrap_th)
# 	return intg.trapz( FThnu_UVthick_IRThin_QuadInt(ths, nu, t, Dist, Aargs, RHStable, Ttable) )
# 	# return intg.trapz( np.array(arg) ) )
# 	# arg = []#np.zeros(len(ths))
# 	# for i in range(len(ths)):
# 	# 	arg.append(FThnu_UVthick_IRThin_QuadInt(ths[i], nu, t, Dist, Aargs, RHStable, Ttable))
# 	# return intg.trapz( np.array(arg) )

# Fnu_UVthick_IRThin_QuadInt = np.vectorize(Fnu_UVthick_IRThin_QuadInt, excluded=[3,4,5])


# def F_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
# 	nus = np.linspace(numin, numax, Ntrap_nu)
# 	#if (type(t) is float or type(t) is np.float64):
# 	return intg.trapz( Fnu_UVthick_IRThin_QuadInt(nus, t, Dist, Aargs, RHStable, Ttable))
# 	# else:
# 	# 	res = []
# 	# 	for i in range(len(t)):
# 	# 		res.append(intg.trapz(Fnu_UVthick_IRThin_QuadInt(nus, t[i], Dist, Aargs, RHStable, Ttable)))
# 	# 	return np.array(res)

# F_ShTorOptThin_Iso_QuadInt = np.vectorize(F_ShTorOptThin_Iso_QuadInt, excluded=[4,5,6])




#########
###THICK
#########
def FThRnu_OptThin_IRThin_QuadInt(thet, r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_OptThin_IRThin_Iso, 0.,2.*ma.pi, args=(thet, r, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def FRnu_OptThin_IRThin_QuadInt(r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThRnu_OptThin_IRThin_QuadInt, 0., ma.pi, args=(r, nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def Fnu_OptThin_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	#Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	n0 = Aargs[4]
	Rd = Aargs[5]
	p  = Aargs[6]
	aeff = Aargs[9]
	Rtau1 = Rd*(  1. - 3.*(p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) #+ 0.01*Rd
	return intg.quad(FRnu_OptThin_IRThin_QuadInt, Rd, Rtau1, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def F_TOptThin_IRThin_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_OptThin_IRThin_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]






















