import numpy as np
from numpy import *
from FluxFuncs_TDEs import *

#PRIORS FROM MEASUREMENTS!
def ln_IR_prior(params):
	#etaR, cosTT, sinJJ, FQfac = params
	etaR, cosTT, sinJJ, nu0, Lav = params
					
	if sinJJ < -1.0 or sinJJ > 1.0:
		return -np.inf

	if cosTT < -1.0 or cosTT > 1.0:
		return -np.inf
					
	if etaR <= 0.0 :
		return -np.inf

	if nu0 <= 0.0 :
	 	return -np.inf

	if Lav <= 0.0001 or Lav > 100.0 :
	 	return -np.inf
			
	return 0.


def ln_V_prior(params):
	#sinJJ, cosTT, Rin, alpha = params
	L0, t0, tfb, gam = params
					
	if L0 <= 0.0 or L0 > 10.0:  ##in units of 10^45 erg/s  bol
		return -np.inf

	if t0 < 0.0 or t0 > 2.0: ##dont shift more than 4 years
		return -np.inf
					
	if tfb < 0.5 or tfb > 5.:
		return -np.inf

	if gam < 0.1 or gam > 4.0:
		return -np.inf

	# if FQfac <= 90.0 :
	# 	return -np.inf
			
	return 0.


##LIKELIHOODS
def ln_IR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	return -(IRTDE_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		

def ln_IR_fxdR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	return -(IRTDE_fxdR_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		


def IRLC_point(p, t, args, RHStable, Ttable):
	etaR, cosT, sinJJt, nu0, Lav = p
	FQfac = 100.0
	nu0=nu0*numicron
	Lav = Lav*10.**(45)
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)

#[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FIR_gal, gam = args
	## Rde doesnt matter here
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_ShTorOptThin_Iso_QuadInt(numn, numx, t, Dst, IRargs, RHStable, Ttable) + FIR_gal)/FRel)





def IRLC_fxdR_point(p, t, args, RHStable, Ttable):
	#etaR, cosT, JJt, FQfac = p
	Rde, cosT, sinJJt, nu0, Lav = p
	FQfac = 100.0
	nu0=nu0*numicron
	Lav = Lav*10.**(45)
	Rde = Rde*pc2cm
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)

 #[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, etaR, FIR_gal, gam = args
	## etaR doesnt matter here
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_fxdR_ShTorOptThin_Iso_QuadInt(numn, numx, t, Dst, IRargs, RHStable, Ttable) + FIR_gal)/FRel)




def IRTDE_Err2(p, t, argW1, argW2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	
	chiW1 = ( y1 - np.minimum(IRLC_point(p, t, argW1, RHStable, Ttable), 12.9)  )/ dy1

	chiW2 = ( y2 - np.minimum(IRLC_point(p, t, argW2, RHStable, Ttable), 11.26) )/ dy2

	print "chiW1^2 = ", sum(chiW1*chiW1)

	print "chiW1 = ", chiW1

	print "chiW2^2 = ", sum(chiW2*chiW2) 

	print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) 
	print(chi2)
	return chi2



def IRTDE_Err2_fmin(p, t, argW1, argW2, RHStable, Ttable, y1, dy1, y2, dy2):
	ln_p = -ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_Err2(p, t, argW1, argW2, RHStable, Ttable, y1, dy1, y2, dy2)
	return ln_l + ln_p



def IRTDE_fxdR_Err2(p, t, argW1, argW2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	
	chiW1 = (  y1 - np.minimum(IRLC_fxdR_point(p, t, argW1, RHStable, Ttable), 12.9)   )/ dy1

	chiW2 = (  y2 - np.minimum(IRLC_fxdR_point(p, t, argW2, RHStable, Ttable), 11.26)  )/ dy2

	print "chiW1^2 = ", sum(chiW1*chiW1)

	print "chiW1 = ", chiW1

	print "chiW2^2 = ", sum(chiW2*chiW2) 

	print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) 
	print(chi2)
	return chi2


def IRTDE_fxdR_Err2_fmin(p, t, argW1, argW2, RHStable, Ttable, y1, dy1, y2, dy2):
	ln_p = -ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_fxdR_Err2(p, t, argW1, argW2, RHStable, Ttable, y1, dy1, y2, dy2)
	return ln_l + ln_p




def ln_V_likelihood(p, t, arg, y, dy):
	return -(VTDE_Err2(p, t, arg, y, dy))




def VLC_point(p, t, arg, LavIRfit):
	L0, t0, tfb, gam = p
	FQfac = 100.0
	L0 = L0*10.**45
	t0 = t0*yr2sec
	tfb = tfb*yr2sec
	Dst, Frel, FV_gal = arg
	BC = 1.0#12.*5./2.*LavIRfit ## bol corr to Vband - doesnt matter for fits, jsut makes Vband correct mag in end
	#return -2.5*np.log10(( Fsrc_Anl_Fit(t, Dst, L0, tfb, t0, gam, FQfac)/BC + FV_gal)/Frel)
	return -2.5*np.log10(( Fsrc_Anl(t, Dst, L0/BC, tfb, t0, gam, FQfac) + FV_gal)/Frel)

#Fsrc_Anl(t, r, Lavg, tfb, t0, FQfac)

def VTDE_Err2(p, t, arg, y, dy):
	print "EVAL", p
	chi = (y - np.minimum(VLC_point(p, t, arg, 1.0),17.6) )/ dy
	chi2 = sum(chi*chi) 
	print(chi2)
	return chi2







##POSTERIORS
def ln_IR_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	ln_p = ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
	return ln_l + ln_p

def ln_IR_fxdR_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	ln_p = ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_fxdR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
	return ln_l + ln_p




def ln_V_posterior(p, t, arg, y, dy):
	ln_p = ln_V_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_V_likelihood(p, t, arg, y, dy)
	return ln_l + ln_p





