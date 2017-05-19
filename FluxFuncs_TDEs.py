import numpy as np
import math as ma
#import numexpr as ne
import scipy as sc

#from scipy.optimize import brentq #fmin

import scipy.integrate as intg
import scipy.signal as sgn
#from apw import Tfunc


####GET args to pass with out passing arguments to functions
#from MCMC_IR import Args2Pass

#Goal is to output the IR flux from irradiated dust geometry surrounding TDE.


# ###FOR TRAP INT
Ntrap_ph = 181#181
Ntrap_th = 81#81
Ntrap_nu = 31#31
N_RHS = 101

#### INTEGRATION ERROR TOLS
myrel = 1.e-8
myabs = 1.e-8#1.e-60
reclim = 10#100
limlst = 10
maxp1 = 10
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
day2sec = 3600.*24.

zTDE = 0.11 ##HARDCODE FOR NOW



## Back body specific intensity
def Bv(nu, T):
	#return 2.*h*nu*nu*nu/(c*c)*1./(np.exp( h*nu/(kb*T) ) - 1.)
	return 2.*h*nu*nu*nu/(c*c)*1./(np.expm1( h*nu/(kb*T) ) )

## Dust absorption/emission efficiency
def Qv(nu, nu0, nn):
	qv = (nu/nu0)**(nn)
	qv = np.minimum(qv, 1.0*nu/nu)
	return qv


def QvBv(nu, T, nu0, nn):
	qv = (nu/nu0)**(nn)
	qv = np.minimum(qv, 1.0*nu/nu)
	#return 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*T)  ) - 1.) * qv
	return 2.*h*nu*nu*nu/(c*c)*1./(np.expm1(  h*nu/(kb*T)  ) ) * qv
	       


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
def T_RHS_Quad(Td, nu0, nn):
	# 4 for difference in cross sectional area and surface area, pi for isotropic flux from Grain
	RHS = 4.*ma.pi* intg.quad(QvBv  ,0., 3.0*numicron, args=(Td, nu0, nn), epsabs=0.0 )[0] 
	return np.array(RHS)

def T_RHS(nus, Td, nu0, nn):
	# 4 for difference in cross sectional area and surface area, pi for isotropic flux from Grain
	return 4.*ma.pi* np.trapz(np.exp(nus)*QvBv(np.exp(nus), Td, nu0, nn),  nus, axis=0 ) 


# def Temp_min(Td, t, r, Lavg, tfb, t0, gam, FQfac, nu0, nn):
# 	return Fsrc_Anl(t, r, Lavg, tfb, t0, gam, FQfac)**2-T_RHS(Td, nu0, nn)**2


# def Temp_Num(t, r, Lavg, tfb, t0, gam, FQfac, nu0, nn):
# 	return sc.brentq(Temp_min, args=(t, r, Lavg, tfb, t0, gam, FQfac, nu0, nn), (Td,100.,1800.))








### optical depth - IGNORE FOR NOW
def tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn):
	xe     = (Rout*Rout - (z*z + y*y))**(0.5)
	return ma.pi*aeff*aeff * Qv(nu, nu0, nn) * intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	
















# def FfQ(t, FQ):
# 	return FQ
# def Ftgt0(t, tfb, t0, gam, F0, FQ):
# 	return np.maximum(F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam),FQ)
####################################################
### optical/UV source 
####################################################
def Fsrc_Anl(t, r, Lavg, tfb, t0, gam, FQfac):
	F0 = Lavg/(4.*ma.pi*r*r)
	FQ = 5.0*10.**40/(4.*ma.pi*r*r) * 1./0.042 ##Xray Lum from Kaite and BC from https://ned.ipac.caltech.edu/level5/March04/Risaliti/Risaliti2_5.html
	#gam = 2.3

	# Fsrc = np.piecewise(  t, [t<=t0, t>t0], [lambda t:FQ, lambda t:np.maximum(F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam),FQ)]  )
	# # Fsrc = np.piecewise(  t, [t<=t0, t>t0], [FQ, Ftgt0(t, tfb, t0, gam, F0, FQ)]  )
	# return Fsrc 

	#return F0*np.sin(2.*ma.pi*(t-t0)/(5.*tfb) )  #exp( -((t-t0)/tfb)**2 ) *  ( (t-t0)/tfb )**(-gam)
	

	if (type(t) is float or type(t) is np.float64):
		#if ((t-t0)<=0.0):
		# tLEQt0    = 0.5*(np.sign(t0 - t) + 1.0 )  
		# tGEQt0    = 0.5*(np.sign(t - t0) + 1.0 ) 
		# FTDE = np.maximum( F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)* tGEQt0/tfb)**(-gam), FQ )
		
		# Fsrc = tLEQt0 * FQ + tGEQt0 * FTDE
		# return np.nan_to_num(Fsrc)
		
		# else:
		# 	Fsrc = np.maximum( F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam), FQ )
		# 	# Fsrc =  F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam)
		# 	# if (Fsrc<=0.0):
		# 	# 	Fsrc=FQ
		# 	return Fsrc 
		return np.maximum(  np.nan_to_num(   F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam)   ), FQ ) 
		#return np.nan_to_num(np.maximum( F0*( 1.+np.log( (t-t0)/tfb ) ) * ((t-t0)/tfb)**(-gam), FQ ))

	else:
		NL = len(t)
		#res = np.empty(NL)
		res = []
		for i in range (NL):
		# 	tLEQt0    = 0.5*(np.sign(t0 - t[i]) + 1.0 ) 
		# 	tGEQt0    = 0.5*(np.sign(t[i] - t0) + 1.0 )  
		# 	FTDE = np.maximum( F0*(1.+np.log( (t[i]-t0)/tfb) ) * ((t[i]-t0)* tGEQt0/tfb)**(-gam), FQ )
		# 	res.append(tLEQt0 * FQ + tGEQt0 * FTDE)
		# return np.array(np.nan_to_num(res))


			# res = np.array(res)

			# if (t[i]-t0<=0.0):
			# 	Fsrc = FQ
			# else:
			# Fsrc = np.maximum(F0*(1.+np.log( (t[i]-t0)/tfb) ) * ((t[i]-t0)/tfb)**(-gam), FQ)
				# Fsrc = F0*(1.+np.log( (t[i]-t0)/tfb) ) * ((t[i]-t0)/tfb)**(-gam)
				# if (Fsrc<=0.0):
				# 	Fsrc=FQ

		# 	res.append(Fsrc)
		
			#res[i] = np.nan_to_num(np.maximum( F0*(1.+np.log( (t[i]-t0)/tfb) ) * ((t[i]-t0)/tfb)**(-gam), FQ ))  

			res.append(   np.maximum(  np.nan_to_num(   F0*(1.+np.log( (t[i]-t0)/tfb) ) * ((t[i]-t0)/tfb)**(-gam)   ), FQ )  ) 

			#res.append(  np.nan_to_num(np.maximum( F0*(1.+np.log( (t[i]-t0)/tfb) ) * ((t[i]-t0)/tfb)**(-gam), FQ ))  )
		return np.array(res)
		



# 	if ((t-t0)<=0.0):
# 		return FQ
# 	else:
# 		Fsrc = F0*(1.+np.log( (t-t0)/tfb) ) * ((t-t0)/tfb)**(-gam)
# 		if (Fsrc<=0.0):
# 			Fsrc=FQ
# 		return Fsrc

#Fsrc_Anl = np.vectorize(Fsrc_Anl)


def Fsrc_Anl_Fit(t, r, Lavg, tfb, t0, gam, FQfac):
	F0 = Lavg/(4.*ma.pi*r*r)
	FQ = 5.*10.**40/(4.*ma.pi*r*r) * 1./0.042 ##Xray Lum from Kaite and BC from https://ned.ipac.caltech.edu/level5/March04/Risaliti/Risaliti2_5.html
	#gam = 2.3
	#FQ = F0/100.0
	#gam = 2.3
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
#def TDust(t,r,thet, ph, args, RHStable, Ttable, RHS_mx, RHS_mn):
def TDust(t,r,thet, ph, args, RHStable, Td_interp, RHS_mx, RHS_mn):
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args

	sinth = np.sin(thet)
	cosJ = np.cos(JJ)
	sinJ = np.sin(JJ)
	x = r*sinth*np.cos(ph)
	y = r*sinth*np.sin(ph)
	z = r*np.cos(thet)
	
	xrot = x*cosJ + z*sinJ
	zrot = z*cosJ - x*sinJ
	
	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	#Td = 1.0*r/r  ##T=1 is very small

	#if (throt>=thetT and throt<=(ma.pi - thetT) and r>=Rd ):
	#Thf_p = float( int(0.5*(np.sign(throt - thetT) + 1.0 ) ) ) 
	#Thf_m = float( int(0.5*(np.sign((ma.pi - thetT) - throt) + 1.0 ) ) )
	# Rf    = float( int(0.5*(np.sign(r - Rd) + 1.0 ) ) )

	Thf_p =  0.5*(np.sign(throt - thetT) + 1.0 ) 
	Thf_m =  0.5*(np.sign((ma.pi - thetT) - throt) + 1.0 ) 
	xbak  =  0.5*(np.sign(-xrot) + 1.0 ) 

		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###
	Fsrc = Fsrc_Anl(t, r, Lavg, tfb, t0, gam, FQfac)

	#0.5*(np.sign(throt - thetT) + 1.0 ) 
	#0.5*(np.sign((np.pi - thetT) - throt) + 1.0 )


	#upp = 0.5*(np.sign(Fsrc - RHS_mx) + 1.0 ) 
	#dwn = 0.5*(np.sign(RHS_mn-Fsrc) + 1.0 ) 

	#Td = 1. * upp

##FOR QUAD
	# if (Fsrc > RHS_mx or Fsrc <= RHS_mn):
	# 	Td = 1.
	# else:
	# #for i in range(len(RHStable)-1):
	# 	istar = np.where( Fsrc <= RHStable )[0].min()  ### Fsrc does not depend on theta phi or nu!
	# 	Td = Ttable[istar]
	# return Td

##FOR TRAP
	Td = np.ones(np.shape(Fsrc))
	
	#Fflat = Fsrc.ravel()
	#then np.where() on 1D array
	#then unravel with correct shape
	##
	# Make 1D interp -> Tdust(RHS) -> Tdust (Fsrc gives answer)
	if ( JJ > (np.pi/2. - thetT) or JJ < -(np.pi/2. - thetT)):
		Fsrc_Td = np.minimum(np.maximum(Fsrc, RHS_mn), RHS_mx)
		Td = Td_interp(Fsrc)  #Td_Interp(Fsrc)
		# for i in range( np.shape(Fsrc)[0]-1 ):
		# 	for j in range( np.shape(Fsrc)[1]-1 ):
		# 		for k in range( np.shape(Fsrc)[2]-1 ):					
		# 			if (xrot[i][j][k] > 0.0 or throt[i][j][k] < thetT or throt[i][j][k]>(np.pi - thetT) or Fsrc[i][j][k] >= RHS_mx or Fsrc[i][j][k] <= RHS_mn ):  #)         or       (throt[i][j][k] > -thetT and throt[i][j][k] < 0.0) or (throt[i][j][k] < -(ma.pi - thetT) and throt[i][j][k] < 0.0)):
		# 				Td[i][j][k] = 1.0


	# if ( JJ > (np.pi/2. - thetT) or JJ < -(np.pi/2. - thetT) ):
	# 	for i in range( np.shape(Fsrc)[0]-1 ):
	# 		for j in range( np.shape(Fsrc)[1]-1 ):
	# 			for k in range( np.shape(Fsrc)[2]-1 ):
	# 				istar = np.where( Fsrc[i][j][k] <= RHStable )[0].min()  #DO this faster?
	# 				Td[i][j][k] = Ttable[istar]
					
	# 				if (xrot[i][j][k] > 0.0 or throt[i][j][k] < thetT or throt[i][j][k]>(np.pi - thetT) or Fsrc[i][j][k] >= RHS_mx or Fsrc[i][j][k] <= RHS_mn ):  #)         or       (throt[i][j][k] > -thetT and throt[i][j][k] < 0.0) or (throt[i][j][k] < -(ma.pi - thetT) and throt[i][j][k] < 0.0)):

	# 				#if ( throt[i][j][k] < thetT or throt[i][j][k]>(np.pi - thetT) or Fsrc[i][j][k] >= RHS_mx or Fsrc[i][j][k] <= RHS_mn ):#)         or       (throt[i][j][k] > -thetT and throt[i][j][k] < 0.0) or (throt[i][j][k] < -(ma.pi - thetT) and throt[i][j][k] < 0.0)):
	# 					Td[i][j][k] = 1.0

					# if (throt[i][j][k] < 0.0 or throt[i][j][k] > ma.pi):
					# 	print "thetrot out of range"


	#return Td

	return np.minimum(np.maximum(Td * Thf_p * Thf_m * xbak, 1.0), 1800.0)
	#return np.maximum(Td * Thf_p * Thf_m, 1.0)









def TDust_old(t,r,thet, ph, args, RHStable, Ttable):
#TDust(tem, Rd, thet, ph, args, RHStable, Ttable)
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args
	#Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args
	#Rd = 0.0
	#Loft = Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac) * 4.*ma.pi*r*r
	#r = etaR * 0.5 * pc2cm * (Loft/(10.**(46)))**(0.5) ##sublimatin front


	sinth = np.sin(thet)
	cosJ = np.cos(JJ)
	sinJ = np.sin(JJ)
	x = r*sinth*np.cos(ph)
	y = r*sinth*np.sin(ph)
	z = r*np.cos(thet)
	
	xrot = x*cosJ + z*sinJ
	zrot = z*cosJ - x*sinJ
	
	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	#Td = 1.0*r/r  ##T=1 is very small

	if (throt>=thetT and throt<=(ma.pi - thetT) and r>=Rd ):
		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###
		Fsrc = Fsrc_Anl(t, r, Lavg, tfb, t0, gam, FQfac)
		


		###-----------------###
		### Compute taudust ###
		###-----------------###
		#Qbar=1. 
		#if r=Rd, tau is just 0
		#tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
		
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
			return 1.
		else:
			istar = np.where( LHS <= RHStable )[0].min()
			return Ttable[istar]

	return r/r

#TDust_old = np.vectorize(TDust_old, excluded=[4,5,6])

####################################################
### T is analytic when Q_nu=1
####################################################
def TDust_Anl(t,r,thet, ph, args):
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args

	#Loft = Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac) * 4.*ma.pi*r*r
	#r = etaR * 0.5 * pc2cm * (Loft/(10.**(46)))**(0.5) ##sublimatin front

	
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	#Td = 1.*r/r  ##T=1 is very small


	Thf_p = 0.5*(np.sign(throt - thetT) + 1.0 ) 
	Thf_m = 0.5*(np.sign((np.pi - thetT) - throt) + 1.0 )
	#Rf    = 0.5*(np.sign(r+100.0 - Rd) + 1.0 ) 
	#if (throt>=thetT and throt<=(ma.pi - thetT) and r>=Rd ):
	###-----------------###
	### COMPUTE Fsrc    ###
	###-----------------###

	Fsrc = Fsrc_Anl(t, r, Lavg, tfb, t0, gam, FQfac)
	#Fsrc = Fsrc_Anl_subl(t, r, Lavg, tfb, t0, FQfac)


	###-----------------###
	### Compute taudust ###
	###-----------------###
	#Qbar=1. 
	#tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
	### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
	

	# nuLEQnu0 = h/kb * ( Fsrc*c*c*nu0**(gam)/(8.*ma.pi * h * special.Gamma[4.+nn] * special.PolyLog[4.+gam, 1]) )**(1./4.+nn)
	# nuGEQnu0 = (Fsrc/(4.*sigSB))**(1./4.)

	# Td = nuGEQnu0*0.5*(np.sign(nu-n0)+1.0) + nuLEQnu0* 0.5*(np.sign(nu0-nu)+1.0)

	#(h (8 \[Pi])^(1/(-4 - gam)) ((c^2 Fsrc (1/nu0)^-gam)/(h Gamma[4 + gam] PolyLog[4 + gam, 1]))^(-(1/(-4 - gam))))/K

	#LHS = Fsrc #* np.exp(-tauDust)
	Td = (Fsrc/(4.*sigSB))**(1./4.)

	return np.minimum(np.maximum(Td * Thf_p * Thf_m, 1.0), 1800.0) #

#TDust_Anl = np.vectorize(TDust_Anl, excluded=[4])


####################################################
### This is the spherical shell case resulting  
### from the assumption that dust is Opt Thick 
### to optical/UV and Opt Thin to IR
####################################################

def Fnuint_UVthick_IRThin_Iso(ph, thet, nu, t, Dist, args, RHStable, Ttable, RHS_mx, RHS_mn):
	### allow dust to track luminosity increse and decrease
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args
	Loft = Fsrc_Anl(t, Rd, Lavg, tfb, t0, gam, FQfac) * 4.*ma.pi*Rd*Rd  ##Rd does not matter here! because Fsrc has (/4.*ma.pi*Rd*Rd in it)
	#Rd = etaR * 0.5 * pc2cm *  (Loft/(10.**(46)))**(0.5)

	# ### allow dust to only to track luminosity increase
	tpeak = t0 + np.exp(1./gam - 1.) * tfb
	Lpeak = Fsrc_Anl(tpeak, Rd, Lavg, tfb, t0, gam, FQfac) * 4.*ma.pi*Rd*Rd

	tLEQpk = 0.5 * (np.sign(tpeak-t) + 1.0)
	tGEQpk = 0.5 * (np.sign(t-tpeak) + 1.0)
	Rd = etaR * 0.5 * pc2cm * (   (Loft/(10.**(46)))**(0.5) * tLEQpk  +  (Lpeak/(10.**(46)))**(0.5) * tGEQpk )

	

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted from dust
	sinth = np.sin(thet)
	tem = t - Rd/c*(1. - sinth*np.cos(ph)) ## t is at z=0 / rest frame but still observed light travel time



	# Tdust - even though working in rest frame time, need to redshift the emitted dust BBs
	Tdust = TDust(tem, Rd, thet, ph,  args, RHStable, Ttable, RHS_mx, RHS_mn) / (1.+zTDE)
	#Tdust = TDust_Anl(tem, Rd, thet, ph, args) / (1.+zTDE)


	# if (Tdust==1.0/(1.+zTDE)):
	# 	fint=1.e-16
	# else:
		# surface density in optically thick limit
	Surf_nd = 1./(ma.pi*aeff*aeff)
	# Rd is the inner edge of the shell
	#fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.expm1(  h*nu/(kb*Tdust)  ) )	
	fint = fint* Rd*Rd* sinth * Surf_nd #mult by spherical Jacobian and surface denisty


	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist * fint


#Fnuint_UVthick_IRThin_Iso = np.vectorize(Fnuint_UVthick_IRThin_Iso, excluded=[5,6,7])

def Fnuint_fxdR_UVthick_IRThin_Iso(ph, thet, nu, t, Dist, args, RHStable, Ttable, RHS_mx, RHS_mn):
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args
  #[Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]

	#Loft = Fsrc_Anl(t, Rd, Lavg, tfb, t0, FQfac) * 4.*ma.pi*Rd*Rd
	#Rd = etaR * 0.5 * pc2cm * (Loft/(10.**(46)))**(0.5)
###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted from dust
	sinth = np.sin(thet)
	tem = t - Rd/c*(1. - sinth*np.cos(ph))



	# Tdust 
	Tdust = TDust(tem, Rd, thet, ph, args, RHStable, Ttable, RHS_mx, RHS_mn) / (1.+zTDE)
	#Tdust = TDust_Anl(tem, Rd, thet, ph, args)/ (1.+zTDE)


	# if (Tdust==1.0/(1.+zTDE)):
	# 	fint=1.e-16
	# else:




		# surface density in optically thick limit
	Surf_nd = 1./(ma.pi*aeff*aeff)
	# Rd is the inner edge of the shell
	#fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.expm1(  h*nu/(kb*Tdust)  ) )	
	fint = fint* Rd*Rd* sinth * Surf_nd #mult by spherical Jacobian and surface denisty


	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist * fint


#Fnuint_fxdR_UVthick_IRThin_Iso = np.vectorize(Fnuint_fxdR_UVthick_IRThin_Iso, excluded=[5,6,7])


#################################################################
### For dust which is optically thin to opt/UV ##IGNORE FOR NOW
#################################################################
## Allow the optical/UV to penetrate the dust, so optically think to opical/UV out to some distance
## must still assume opt thin to IR
def Fnuint_OptThin_IRThin_Iso(ph, thet, r, nu, t, Dist, args, RHStable, Ttable): #, tauGrid):	
	Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args
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
	Tdust = TDust(tem, Rd, thet, ph, args, RHStable, Ttable)/(1.+zTDE)
	#Tdust = TDust_Anl(tem, Rd, thet, ph, args)

	# fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.expm1(  h*nu/(kb*Tdust)  ) )
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
# #########
###QUAD INT
def FThnu_UVthick_IRThin_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
	return intg.quad(Fnuint_UVthick_IRThin_Iso, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]

def Fnu_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
	return intg.quad(FThnu_UVthick_IRThin_QuadInt, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

def F_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
	if (type(t) is float or type(t) is np.float64):
		return intg.quad(Fnu_UVthick_IRThin_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]
	else:
		res = []
		for i in range(len(t)):
			res.append(intg.quad(Fnu_UVthick_IRThin_QuadInt, numin, numax, args=(t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0])
		return np.array(res)

######################
### WORKING W1 and W2 QUAD
######################
# def F_W1_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus):
# def F_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
# 	if (type(t) is float or type(t) is np.float64):
# 		return intg.quad(Fnu_UVthick_IRThin_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]
# 	else:
# 		res = []
# 		for i in range(len(t)):
# 			res.append(intg.quad(Fnu_UVthick_IRThin_QuadInt, numin, numax, args=(t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0])
# 		return np.array(res)
######################
######################


def FThnu_UVthick_IRThin_TrapInt(thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis):
	#phis = np.linspace(0.0,2.*ma.pi, Ntrap_ph)
	Ntrap_ph = len(phis)
	return 2.*ma.pi/(2.*Ntrap_ph) * (2.0 * np.sum([Fnuint_UVthick_IRThin_Iso(ph,thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) for ph in phis]) - Fnuint_UVthick_IRThin_Iso(phis[0],thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) - Fnuint_UVthick_IRThin_Iso(phis[Ntrap_ph-1],thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) )


	
def Fnu_UVthick_IRThin_TrapInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths):
	Ntrap_th = len(ths)
	#ths = np.linspace(0.0, ma.pi, Ntrap_th)
	return ma.pi/(2.*Ntrap_th) * (2.0 * np.sum([FThnu_UVthick_IRThin_QuadInt(th, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis) for th in ths]) - FThnu_UVthick_IRThin_QuadInt(ths[0], nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis) - FThnu_UVthick_IRThin_QuadInt(ths[Ntrap_th-1], nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis) )

##longnu int
# def F_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
# 	nmn = np.log(numin)
# 	nmx = np.log(numax)
# 	lnnus = np.linspace(nmn, nmx, Ntrap_nu)
# 	return (nmx-nmn)/(2.*Ntrap_nu) * (2.0 * np.sum([2.*np.exp(lnnu)*Fnu_UVthick_IRThin_QuadInt(np.exp(lnnu), t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) for lnnu in lnnus]) - np.exp(lnnus[0])*Fnu_UVthick_IRThin_QuadInt(np.exp(lnnus[0]), t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) - np.exp(lnnus[Ntrap_nu-1])*Fnu_UVthick_IRThin_QuadInt(np.exp(lnnus[Ntrap_nu-1]), t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) )

##reg nu int
def F_ShTorOptThin_Iso_TrapInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths, nus):
	#nus = np.linspace(numin, numax, Ntrap_nu)
	Ntrap_nu = len(nus)
	if (type(t) is float or type(t) is np.float64):
		return (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu) * (2.0 * np.sum([Fnu_UVthick_IRThin_TrapInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - Fnu_UVthick_IRThin_TrapInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - Fnu_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )
	else:
		Nl = len(t) 
		res = np.empty(Nl)
		for i in range(Nl):
			res[i] =  (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu) * (2.0 * np.sum([Fnu_UVthick_IRThin_TrapInt(nu, t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - Fnu_UVthick_IRThin_TrapInt(nus[0], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - Fnu_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )  

		return res


##W1 RSR
# def F_W1_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus):
# 	Ntrap_nu = len(nus)
# 	#nus = np.linspace(numin, numax, Ntrap_nu)
# 	#F0nrm = (numax-numin)/(2.*Ntrap_nu) * (2.0 * np.sum([W1RSR_intrp(nus) for nu in nus]) - W1RSR_intrp(nus[0]) -   W1RSR_intrp(nus[Ntrap_nu-1]) )
# 	if (type(t) is float or type(t) is np.float64):
# 		return (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu) * (2.0 * np.sum([W1RSR_intrp(nu)*Fnu_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W1RSR_intrp(nus[0])*Fnu_UVthick_IRThin_QuadInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W1RSR_intrp(nus[Ntrap_nu-1])*Fnu_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )
# 	else:
# 		Nl = len(t) 
# 		res = np.empty(Nl)
# 		for i in range(Nl):
# 			res[i] = (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu) * (2.0 * np.sum([W1RSR_intrp(nu)*Fnu_UVthick_IRThin_QuadInt(nu, t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W1RSR_intrp(nus[0])*Fnu_UVthick_IRThin_QuadInt(nus[0], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W1RSR_intrp(nus[Ntrap_nu-1])*Fnu_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )  

# 		return res

# ##W1 RSR
# def F_W2_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nus):
# 	#nus = np.linspace(numin, numax, Ntrap_nu)
# 	#F0nrm = (numax-numin)/(2.*Ntrap_nu) * (2.0 * np.sum([W1RSR_intrp(nus) for nu in nus]) - W1RSR_intrp(nus[0]) -   W1RSR_intrp(nus[Ntrap_nu-1]) )
# 	Ntrap_nu = len(nus)
# 	if (type(t) is float or type(t) is np.float64):
# 		return  (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu) * (2.0 * np.sum([W2RSR_intrp(nu)*Fnu_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W2RSR_intrp(nus[0])*Fnu_UVthick_IRThin_QuadInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W2RSR_intrp(nus[Ntrap_nu-1])*Fnu_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )
# 	else:
# 		Nl = len(t) 
# 		res = np.empty(Nl)
# 		for i in range(Nl):
# 			res[i]  = (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu) * (2.0 * np.sum([W2RSR_intrp(nu)*Fnu_UVthick_IRThin_QuadInt(nu, t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W2RSR_intrp(nus[0])*Fnu_UVthick_IRThin_QuadInt(nus[0], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W2RSR_intrp(nus[Ntrap_nu-1])*Fnu_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )  

# 		return res


######################
### WORKING W1 and W2 TRAPZ
######################
def F_W1_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus):
	if (type(t) is float or type(t) is np.float64):
		return np.trapz(np.trapz( np.trapz(W1RSR_intrp(nus[2])*Fnuint_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=1), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
	else:
		res = np.empty(len(t))
		for i in range(len(t)):
			res[i] = np.trapz(np.trapz( np.trapz(W1RSR_intrp(nus[2])*Fnuint_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=1), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
		return res


def F_W2_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nus):
	if (type(t) is float or type(t) is np.float64):
		return np.trapz(np.trapz( np.trapz(W2RSR_intrp(nus[2])*Fnuint_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=1), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
	else:
		res = np.empty(len(t))
		for i in range(len(t)):
			res[i] = np.trapz(np.trapz( np.trapz(W2RSR_intrp(nus[2])*Fnuint_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=1), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
		return res
######################
######################



# 	# if (type(t) is float or type(t) is np.float64):
# 	# 	Trap_sub = Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) + Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) 
# 	# 	#Trap_int = -Trap_sub + ma.fsum(2.*Fnu_fxdR_UVthick_IRThin_QuadInt(nus, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) )
# 	# 	Trap_int = 0.0
# 	# 	for i in range(Ntrap_nu):
# 	# 		Trap_int += 2.*Fnu_UVthick_IRThin_QuadInt(nus[i], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) 
# 	# 	Trap_int -= Trap_sub
# 	# 	return (numax-numin)/(2.*Ntrap_nu) * (Trap_int)
# 	# else:
# 	# 	res = []
# 	# 	for i in range(len(t)):
# 	# 		Trap_sub = Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) + Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) 
# 	# 		Trap_int = 0.0
# 	# 		#Trap_int = -Trap_sub + ma.fsum(2.*Fnu_fxdR_UVthick_IRThin_QuadInt(nus, t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) )
# 	# 		for j in range(Ntrap_nu):
# 	# 			Trap_int += 2.*Fnu_UVthick_IRThin_QuadInt(nus[j], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) 
# 	# 		Trap_int -= Trap_sub
# 	# 		res.append((numax-numin)/(2.*Ntrap_nu) * (Trap_int))
# 	# 	return np.array(res)




##FOR FIXED Rsub
###QUADINT
# def FThnu_fxdR_UVthick_IRThin_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
# 	return intg.quad(Fnuint_fxdR_UVthick_IRThin_Iso, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]

# def Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
# 	return intg.quad(FThnu_fxdR_UVthick_IRThin_QuadInt, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

# def F_fxdR_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
# 	if (type(t) is float or type(t) is np.float64):
# 		return intg.quad(Fnu_fxdR_UVthick_IRThin_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]
# 	else:
# 		res = []
# 		for i in range(len(t)):
# 			res.append(intg.quad(Fnu_fxdR_UVthick_IRThin_QuadInt, numin, numax, args=(t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0])
# 		return np.array(res)


###TRAPINT
def FThnu_fxdR_UVthick_IRThin_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis):
	#phis = np.linspace(0.0,2.*ma.pi, Ntrap_ph)
	Ntrap_ph = len(phis)
	return 2.*ma.pi/(2.*Ntrap_ph) * (2.0 * np.sum([Fnuint_fxdR_UVthick_IRThin_Iso(ph, thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) for ph in phis]) - Fnuint_fxdR_UVthick_IRThin_Iso(phis[0],thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) - Fnuint_fxdR_UVthick_IRThin_Iso(phis[Ntrap_ph-1],thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) )


	
def Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths):
	#ths = np.linspace(0.0, ma.pi, Ntrap_th)
	Ntrap_th = len(ths)
	return ma.pi/(2.*Ntrap_th) * (2.0 * np.sum([FThnu_fxdR_UVthick_IRThin_QuadInt(th, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis) for th in ths]) - FThnu_fxdR_UVthick_IRThin_QuadInt(ths[0], nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis) - FThnu_fxdR_UVthick_IRThin_QuadInt(ths[Ntrap_th-1], nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis) )


##longnu int
# def F_fxdR_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn):
# 	nmn = np.log(numin)
# 	nmx = np.log(numax)
# 	lnnus = np.linspace(nmn, nmx, Ntrap_nu)
# 	return (nmx-nmn)/(2.*Ntrap_nu) * (2.0 * np.sum([2.*np.exp(lnnu)*Fnu_fxdR_UVthick_IRThin_QuadInt(np.exp(lnnu), t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) for lnnu in lnnus]) - np.exp(lnnus[0])*Fnu_fxdR_UVthick_IRThin_QuadInt(np.exp(lnnus[0]), t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) - np.exp(lnnus[Ntrap_nu-1])*Fnu_fxdR_UVthick_IRThin_QuadInt(np.exp(lnnus[Ntrap_nu-1]), t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) )

#reg nu int
def F_fxdR_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths, nus):
	#nus = np.linspace(numin, numax, Ntrap_nu)
	#return (numax-numin)/(2.*Ntrap_nu) * (2.0 * np.sum([Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) for nu in nus]) - Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) - Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn) )
	#nus = np.linspace(numin, numax, Ntrap_nu)
	Ntrap_nu = len(nus)
	if (type(t) is float or type(t) is np.float64):
		return (nus[0]-nus[Ntrap_nu-1])/(2.*Ntrap_nu) * (2.0 * np.sum([Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )
	else:
		#res = []
		Nl = len(t) 
		res = np.empty(Nl)
		for i in range(Nl):
			res[i] = (nus[0]-nus[Ntrap_nu-1])/(2.*Ntrap_nu)  * (2.0 * np.sum([Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )  
		return res


# def F_W1_fxdR_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus):
# 	Ntrap_nu = len(nus)
# 	if (type(t) is float or type(t) is np.float64):
# 		return  (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu)  * (2.0 * np.sum([W1RSR_intrp(nu)*Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W1RSR_intrp(nus[0])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W1RSR_intrp(nus[Ntrap_nu-1])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )
# 	else:
# 		#res = []
# 		Nl = len(t) 
# 		res = np.empty(Nl)
# 		for i in range(Nl):
# 			res[i] = (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu)  * (2.0 * np.sum([W1RSR_intrp(nu)*Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W1RSR_intrp(nus[0])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W1RSR_intrp(nus[Ntrap_nu-1])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )  

# 		return res

# def F_W2_fxdR_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nus):
# 	Ntrap_nu = len(nus)
# 	if (type(t) is float or type(t) is np.float64):
# 		return  (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu)  * (2.0 * np.sum([W2RSR_intrp(nu)*Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W2RSR_intrp(nus[0])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W2RSR_intrp(nus[Ntrap_nu-1])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )
# 	else:
# 		#res = []
# 		Nl = len(t) 
# 		res = np.empty(Nl)
# 		for i in range(Nl):
# 			res[i] = (nus[Ntrap_nu-1]-nus[0])/(2.*Ntrap_nu)  * (2.0 * np.sum([W2RSR_intrp(nu)*Fnu_fxdR_UVthick_IRThin_QuadInt(nu, t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) for nu in nus]) - W2RSR_intrp(nus[0])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[0], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) - W2RSR_intrp(nus[Ntrap_nu-1])*Fnu_fxdR_UVthick_IRThin_QuadInt(nus[Ntrap_nu-1], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths) )  

# 		return res










def F_W1_fxdR_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus):
	if (type(t) is float or type(t) is np.float64):
		return np.trapz(np.trapz( np.trapz(W1RSR_intrp(nus[2])*Fnuint_fxdR_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=0), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
	else:
		res = np.empty(len(t))
		for i in range(len(t)):
			res[i] = np.trapz(np.trapz( np.trapz(W1RSR_intrp(nus[2])*Fnuint_fxdR_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=0), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
		return res


def F_W2_fxdR_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nus):
	if (type(t) is float or type(t) is np.float64):
		return np.trapz(np.trapz( np.trapz(W2RSR_intrp(nus[2])*Fnuint_fxdR_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=0), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
	else:
		res = np.empty(len(t))
		for i in range(len(t)):
			res[i] = np.trapz(np.trapz( np.trapz(W2RSR_intrp(nus[2])*Fnuint_fxdR_UVthick_IRThin_Iso(nus[0], nus[1], nus[2], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), dx=2.*ma.pi/(Ntrap_ph), axis=0), dx=ma.pi/(Ntrap_th), axis=0), nus[2][0][0], axis=0)
		return res


	#res = np.empty(len(t))
	#for i in range(len(t)):
# 			res[i] = np.trapz(np.trapz( np.trapz(W2RSR_intrp(nus[2])*Fnuint_fxdR_UVthick_IRThin_Iso(nus[1], nus[0], nus[2], t[i], Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn), axis=0), axis=0), axis=0) 
# 	return res



### tplquad 316/277 times slower than trap above
# def W1_tpl_fxdR(ph, thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp):
# 	return W1RSR_intrp(nu)*Fnuint_fxdR_UVthick_IRThin_Iso(ph, thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn)


# def F_W1_fxdR_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus):
# 	Ntrap_nu = len(nus)
# 	return intg.tplquad(W1_tpl_fxdR, nus[0], nus[Ntrap_nu-1], lambda nu: 0.0, lambda nu: ma.pi, lambda nu, theta: 0.0, lambda nu, theta: 2.*ma.pi, args =(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp), epsabs=myabs, epsrel=myrel )[0]



# def W2_tpl_fxdR(ph, thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp):
# 	return W2RSR_intrp(nu)*Fnuint_fxdR_UVthick_IRThin_Iso(ph, thet, nu, t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn)


# def F_W2_fxdR_ShTorOptThin_Iso_QuadInt(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nus):
# 	Ntrap_nu = len(nus)
# 	return intg.tplquad(W2_tpl_fxdR, nus[0], nus[Ntrap_nu-1], lambda nu: 0.0, lambda nu: ma.pi, lambda nu, theta: 0.0, lambda nu, theta: 2.*ma.pi, args =(t, Dist, Aargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp), epsabs=myabs, epsrel=myrel )[0]











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






















