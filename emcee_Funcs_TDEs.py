import numpy as np
from numpy import *
from FluxFuncs_TDEs import *







############################
## MAX LIK FNS for Simul FIT of W1 W2 and V0bands
############################
#PRIORS



##ERR2 FUNCS (-LogLik now)
def IRTDE_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	print "EVAL", p

	pIR = [p[0], p[1], p[2], p[3], p[4], p[5]]
	pVb = [p[6], p[7], p[8], p[9]]

	nn = 2.*len(t) + len(tVb)
	
	chiW1 = 0.5*( y1 - np.minimum(IRLC_ML_point(pIR, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn), 12.9)  )/ np.sqrt( dy1*dy1 + p[5]*p[5]  )

	chiW2 = 0.5*( y2 - np.minimum(IRLC_ML_point(pIR, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn), 11.26) )/ np.sqrt( dy2*dy2 + p[5]*p[5]  )

	chiVb = 0.5*(yVb - np.minimum(VLC_point(pVb, tVb, argVb, 1.0),17.6) )/  np.sqrt( dyVb*dyVb + p[5]*p[5]  )



	# print "chiW1^2 = ", sum(chiW1*chiW1)

	# print "chiW1 = ", chiW1

	# print "chiW2^2 = ", sum(chiW2*chiW2) 

	# print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) + sum(chiVb*chiVb)

	penlty = 0.5 * ( sum( np.log(dy1*dy1 + p[5]*p[5]) ) + sum( np.log(dy2*dy2 + p[5]*p[5]) ) + sum( np.log(dyVb*dyVb + p[5]*p[5]) )  )

	LogLik = chi2/2. + nn/2. * np.log(2.*ma.pi) + penlty


	##LogLik i made negative in liklihood function
	print(-LogLik)
	return LogLik

def IRTDE_fxdR_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	print "EVAL", p

	pIR = [p[0], p[1], p[2], p[3], p[4], p[5]]
	pVb = [p[6], p[7], p[8], p[9]]

	nn = 2.*len(t) + len(tVb)
	
	chiW1 = 0.5*( y1 - np.minimum(IRLC_fxdR_ML_point(pIR, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn), 12.9)  )/ np.sqrt( dy1*dy1 + p[5]*p[5]  )

	chiW2 = 0.5*( y2 - np.minimum(IRLC_fxdR_ML_point(pIR, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn), 11.26) )/ np.sqrt( dy2*dy2 + p[5]*p[5]  )

	chiVb = 0.5*(yVb - np.minimum(VLC_point(pVb, tVb, argVb, 1.0),17.6) )/  np.sqrt( dyVb*dyVb + p[5]*p[5]  )




	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) + sum(chiVb*chiVb)

	penlty = 0.5 * ( sum( np.log(dy1*dy1 + p[5]*p[5]) ) + sum( np.log(dy2*dy2 + p[5]*p[5]) ) + sum( np.log(dyVb*dyVb + p[5]*p[5]) )  )

	LogLik = chi2/2. + nn/2. * np.log(2.*ma.pi) + penlty


	##LogLik i made negative in liklihood function
	print(-LogLik)
	return LogLik



##likliehoods (- Err2 Funcs)
def ln_IR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	return -(IRTDE_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb))
		

def ln_IR_fxdR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	return -(IRTDE_fxdR_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb))
		





##POSTERIORS
def ln_IR_ALL_posterior(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	pIR = [p[0], p[1], p[2], p[3], p[4], p[5]]
	pVb = [p[6], p[7], p[8], p[9]]
	ln_p = ln_IR_ML_prior(pIR) + ln_V_prior(pVb)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p

def ln_IR_fxdR_ALL_posterior(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	pIR = [p[0], p[1], p[2], p[3], p[4], p[5]]
	pVb = [p[6], p[7], p[8], p[9]]
	ln_p = ln_IR_ML_prior(pIR) + ln_V_prior(pVb)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_fxdR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p







#Fmin
def IRTDE_fxdR_FitAll_Err2_fmin(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	pIR = [p[0], p[1], p[2], p[3], p[4], p[5]]
	pVb = [p[6], p[7], p[8], p[9]]
	ln_p = ln_IR_ML_prior(pIR) + ln_V_prior(pVb)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_fxdR_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p


def IRTDE_sblR_FitAll_Err2_fmin(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb):
	pIR = [p[0], p[1], p[2], p[3], p[4], p[5]]
	pVb = [p[6], p[7], p[8], p[9]]
	ln_p = ln_IR_ML_prior(pIR) + ln_V_prior(pVb)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p





















############################
##Chi^2 Min Liks
############################
#PRIOR
def ln_IR_prior(params):
	#etaR, cosTT, sinJJ, FQfac = params
	etaR, cosTT, sinJJ, nu0, Lav, sigML = params
					
	if (sinJJ < -1.0 or sinJJ > 1.0):
		return -np.inf

	if (cosTT < -1.0 or cosTT > 1.0):
		return -np.inf
					
	if (etaR <= 0.0 or etaR >=100.0):
		return -np.inf

	if nu0 <= 0.0:
	 	return -np.inf

	if (Lav > 100.0 or Lav <= 0.0001):
	 	return -np.inf
			
	return 0.








### point calcs to go to Err funcs
def IRLC_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn):
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
	return -2.5*np.log10( (F_ShTorOptThin_Iso_QuadInt(numn, numx, t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn) + FIR_gal)/FRel)




def IRLC_fxdR_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn):
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
	return -2.5*np.log10( (F_fxdR_ShTorOptThin_Iso_QuadInt(numn, numx, t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn) + FIR_gal)/FRel)





### chi^2 funcs
def IRTDE_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	print "EVAL", p
	
	chiW1 = 0.5*( y1 - np.minimum(IRLC_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn), 12.9)  )/ dy1

	chiW2 = 0.5*( y2 - np.minimum(IRLC_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn), 11.26) )/ dy2

	# print "chiW1^2 = ", sum(chiW1*chiW1)

	# print "chiW1 = ", chiW1

	# print "chiW2^2 = ", sum(chiW2*chiW2) 

	# print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) 
	print(chi2)
	return chi2





def IRTDE_fxdR_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	print "EVAL", p
	
	chiW1 = 0.5*(  y1 - np.minimum(IRLC_fxdR_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn), 12.9)   )/ dy1

	chiW2 = 0.5*(  y2 - np.minimum(IRLC_fxdR_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn), 11.26)  )/ dy2 

	# print "chiW1^2 = ", sum(chiW1*chiW1)

	# print "chiW1 = ", chiW1

	# print "chiW2^2 = ", sum(chiW2*chiW2) 

	# print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) 
	print(chi2)
	return chi2






##LIKELIHOODS (- Err2 funcs)
def ln_IR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	return -(IRTDE_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2))
		

def ln_IR_fxdR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	return -(IRTDE_fxdR_ML_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2))
		


#POSTERIOS
def ln_IR_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p

def ln_IR_fxdR_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_fxdR_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p





#Fmin
def IRTDE_fxdR_Err2_fmin(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = -ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_fxdR_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p


def IRTDE_Err2_fmin(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = -ln_IR_prior(p)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p


############################
############################
























############################
## MAX LIK FNS
############################
#PRIORS
def ln_IR_ML_prior(params):
	#etaR, cosTT, sinJJ, FQfac = params
	etaR, cosTT, sinJJ, nu0, Lav, sigML = params
					
	if (sinJJ < -1.0 or sinJJ > 1.0):
		return -np.inf

	if (cosTT < -1.0 or cosTT > 1.0):
		return -np.inf
					
	if (etaR <= 0.0 or etaR >=100.0):
		return -np.inf

	if nu0 <= 0.0:
	 	return -np.inf

	if (Lav > 100.0 or Lav <= 0.0001):
	 	return -np.inf

	if sigML < 0.0:
	 	return -np.inf	
			
	return 0.

### Point calculations
def IRLC_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn):
	etaR, cosT, sinJJt, nu0, Lav, sigML = p
	FQfac = 100.0
	nu0=nu0*numicron
	Lav = Lav*10.**(45)
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)

#[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FIR_gal, gam = args
	## Rde doesnt matter here
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_ShTorOptThin_Iso_QuadInt(numn, numx, t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn) + FIR_gal)/FRel)




def IRLC_fxdR_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn):
	Rde, cosT, sinJJt, nu0, Lav, sigML = p
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
	return -2.5*np.log10( (F_fxdR_ShTorOptThin_Iso_QuadInt(numn, numx, t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn) + FIR_gal)/FRel)





##ERR2 FUNCS (-LogLik now)
def IRTDE_ML_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	print "EVAL", p

	nn = 2.*len(t)
	
	chiW1 = 0.5*( y1 - np.minimum(IRLC_ML_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn), 12.9)  )/ np.sqrt( dy1*dy1 + p[5]*p[5]  )

	chiW2 = 0.5*( y2 - np.minimum(IRLC_ML_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn), 11.26) )/ np.sqrt( dy2*dy2 + p[5]*p[5]  )

	# print "chiW1^2 = ", sum(chiW1*chiW1)

	# print "chiW1 = ", chiW1

	# print "chiW2^2 = ", sum(chiW2*chiW2) 

	# print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) 

	penlty = 0.5 * ( sum( np.log(dy1*dy1 + p[5]*p[5]) ) + sum( np.log(dy2*dy2 + p[5]*p[5]) )  )

	LogLik = chi2/2. + nn/2. * np.log(2.*ma.pi) + penlty


	##LogLik i made negative in liklihood function
	print(-LogLik)
	return LogLik

def IRTDE_fxdR_ML_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	print "EVAL", p

	nn = 2.*len(t)

	chiW1 = 0.5*(  y1 - np.minimum(IRLC_fxdR_ML_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn), 12.9)   )/ np.sqrt( dy1*dy1 + p[5]*p[5] )

	chiW2 = 0.5*(  y2 - np.minimum(IRLC_fxdR_ML_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn), 11.26)  )/ np.sqrt( dy2*dy2 + p[5]*p[5] )

	# print "chiW1^2 = ", sum(chiW1*chiW1)

	# print "chiW1 = ", chiW1

	# print "chiW2^2 = ", sum(chiW2*chiW2) 

	# print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) 

	penlty = 0.5 * ( sum( np.log(dy1*dy1 + p[5]*p[5]) ) + sum( np.log(dy2*dy2 + p[5]*p[5]) )  )

	LogLik = chi2/2. + nn/2. * np.log(2.*ma.pi) + penlty


	##LogLik i made negative in liklihood function
	print(-LogLik)
	return LogLik



##likliehoods (- Err2 Funcs)
def ln_IR_ML_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	return -(IRTDE_ML_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2))
		

def ln_IR_fxdR_ML_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	return -(IRTDE_fxdR_ML_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2))
		





##POSTERIORS
def ln_IR_ML_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = ln_IR_ML_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_ML_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p

def ln_IR_fxdR_ML_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = ln_IR_ML_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_fxdR_ML_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p





## FOR Fmin
def IRTDE_Err2_ML_fmin(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = -ln_IR_ML_prior(p)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_ML_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p	

def IRTDE_fxdR_ML_Err2_fmin(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	ln_p = -ln_IR_ML_prior(p)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_fxdR_ML_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2)
	return ln_l + ln_p
############################
############################














############################
## SOURCE FITTING
############################
#PRIOR
def ln_V_prior(params):
	#sinJJ, cosTT, Rin, alpha = params
	L0, t0, tfb, gam = params
					
	if L0 <= 0.0 or L0 > 10.0:  ##in units of 10^45 erg/s  bol
		return -np.inf

	if t0 < 0.0 or t0 > 2.0: ##dont shift more than 4 years
		return -np.inf
					
	if tfb < 0.0 or tfb > 5.:
		return -np.inf

	if gam < 0.0 or gam > 4.0:
		return -np.inf

	# if FQfac <= 90.0 :
	# 	return -np.inf
			
	return 0.








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
	chi = 0.5*(y - np.minimum(VLC_point(p, t, arg, 1.0),17.6) )/ dy
	chi2 = sum(chi*chi) 
	print(chi2)
	return chi2



def ln_V_likelihood(p, t, arg, y, dy):
	return -(VTDE_Err2(p, t, arg, y, dy))




def ln_V_posterior(p, t, arg, y, dy):
	ln_p = ln_V_prior(p)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_V_likelihood(p, t, arg, y, dy)
	return ln_l + ln_p
############################
############################






