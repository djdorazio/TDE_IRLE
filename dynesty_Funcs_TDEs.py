import numpy as np
from numpy import *
from FluxFuncs_TDEs import *
from scipy.interpolate import interp1d



##
#misc pickles

class interp1d_picklable:
    """ class wrapper for piecewise linear function
    """
    def __init__(self, xi, yi, **kwargs):
        self.xi = xi
        self.yi = yi
        self.args = kwargs
        self.f = interp1d(xi, yi, **kwargs)

    def __call__(self, xnew):
        return self.f(xnew)

    def __getstate__(self):
        return self.xi, self.yi, self.args

    def __setstate__(self, state):
        self.f = interp1d(state[0], state[1], **state[2])


############################
## MAX LIK FNS for Simul FIT of W1 W2 and V0bands
############################
#PRIORS
def ptform(params):
	#sinJJ, cosTT, Rin, alpha = params
	#etaR, cosTT, sinJJ, nu0, kk, Lav, sigMLIR,L0, t0, tfb, gam, sigMLV = params
	etaR, cosTT, sinJJ, nu0, kk, Lav, sigMLIR,L0, t0, tfb, gam = params
					
	#Transform from unit cube to actual param values

	L0out = 10.**( (L0*2.0 - 1.0)*4.0) ##from 10^(-3) to 10^(3)  

	t0out = t0*2.0

	tfbout = 10.**( (tfb*2.0 - 1.0)*1.0)

	gamout = gam*4.0

	etaRout = 10.**( (etaR*2.0 - 1.0)*2.0)

	sinJJout = (sinJJ*2.0 - 1.0)

	cosTTout = (cosTT*2.0 - 1.0)

	nu0out = 10.**( (nu0*2.0 - 1.0)*5.0)  ##from 10^(-5) to 10^(5)

	kkout = 20.0*kk

	Lavout = 10.**( (Lav*2.0 - 1.0)*5.0) ##from 10^(-4) to 10^(4)

	sigMLIRout = sigMLIR*0.5

	#sigMLVout = sigMLV*0.5

	#x = np.array([etaRout, sinJJout, cosTTout, nu0out, kkout, Lavout, sigMLIRout, L0out, t0out, tfbout, gamout, sigMLVout])		
	x = np.array([etaRout, sinJJout, cosTTout, nu0out, kkout, Lavout, sigMLIRout, L0out, t0out, tfbout, gamout])		
	return x














##ERR2 FUNCS (-LogLik now)
def IRTDE_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phi, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	#print "EVAL", p


	ptot = [p[0], p[1], p[2], p[3], p[4], p[5], p[8], p[9], p[10]]  #6 is sgimaIR, 7 is  only for src Lum
	#pVb = [p[6], p[7], p[8], p[9]]
	##wtih fitting kk

	### IF YOU WANT TO VARY nu0 and nnk, have to remake looksup table
	
	#print "Creating look up tables"
	nu0 = numicron*p[3]
	nne = p[4]
	NT = 60
	RHS_table = np.zeros(NT)
	TT_table = np.linspace(1., 1800., NT)
	nus_RHS = np.linspace(1.0, np.log(5.*numicron), N_RHS)
	for i in range(NT):
		RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
	RHS_mx = RHS_table[len(RHS_table)-1]
	RHS_mn = RHS_table[0]


	###MAKE TDUST INTERP
	#print "Making Tdust interp"
	Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table)



	pVb = [p[7], p[8], p[9], p[10]]  ##11 is sigma V band

	nn = 2.*len(t) + len(tVb)
	
	chiW1 = ( y1 - IRLC_W1_ML_point(ptot, t, argW1, RHStable, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phi, ths, nuW1)  ) / np.sqrt( dy1*dy1 + p[6]*p[6]  )

	chiW2 = ( y2 - IRLC_W2_ML_point(ptot, t, argW2, RHStable, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phi, ths, nuW2)  ) / np.sqrt( dy2*dy2 + p[6]*p[6]  )

	##USE same sigML as for IR
	#chiVb = 0.5*(yVb - VLC_point(pVb, tVb, argVb, 1.0) )/  np.sqrt( dyVb*dyVb + p[10]*p[10]  )

	#USE diff sigML
	#chiVb = (yVb - VLC_point(pVb, tVb, argVb, 1.0) ) /  np.sqrt( dyVb*dyVb + p[11]*p[11]  )
	chiVb = (yVb - VLC_point(pVb, tVb, argVb, 1.0) ) /  np.sqrt( dyVb*dyVb  )



	# print "chiW1^2 = ", sum(chiW1*chiW1)

	# # print "chiW1 = ", chiW1

	# print "chiW2^2 = ", sum(chiW2*chiW2) 

	# # print "chiW2 = ", chiW2

	# print "chiV^2 = ", sum(chiVb*chiVb) 

	chi2 = np.sum(chiW1*chiW1) + np.sum(chiW2*chiW2) + np.sum(chiVb*chiVb)

	##USE same sigML as for IR
	#penlty = 0.5 * ( sum( np.log(dy1*dy1 + p[5]*p[5]) ) + sum( np.log(dy2*dy2 + p[5]*p[5]) ) + sum( np.log(dyVb*dyVb + p[5]*p[5]) )  )
	#USE diff sigML
#	penlty = 0.5 * ( np.sum( np.log(dy1*dy1 + p[6]*p[6]) ) + np.sum( np.log(dy2*dy2 + p[6]*p[6]) ) + np.sum( np.log(dyVb*dyVb + p[11]*p[11]) )  )
	penlty = 0.5 * ( np.sum( np.log(dy1*dy1 + p[6]*p[6]) ) + np.sum( np.log(dy2*dy2 + p[6]*p[6]) )  )


	negLogLik = chi2/2. + nn/2. * np.log(2.*ma.pi) + penlty


	##LogLik i made negative in liklihood function
	#print(-negLogLik)
	return negLogLik

def IRTDE_fxdR_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp,  phi, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	print "EVAL", p

	#pIR = [p[0], p[1], p[2], p[3], p[4], p[5]]
	#pVb = [p[6], p[7], p[8], p[9]]


	ptot = [p[0], p[1], p[2], p[3], p[4], p[5], p[8], p[9], p[10]]  #6 is sgimaIR
	#pVb = [p[6], p[7], p[8], p[9]]
	##wtih fitting kk
	nu0 = numicron*p[3]
	nne = p[4]
	NT = 60
	RHS_table = np.zeros(NT)
	TT_table = np.linspace(1., 1800., NT)
	nus_RHS = np.linspace(1.0, np.log(5.*numicron), N_RHS)
	for i in range(NT):
		RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
	RHS_mx = RHS_table[len(RHS_table)-1]
	RHS_mn = RHS_table[0]


	###MAKE TDUST INTERP
	#print "Making Tdust interp"
	Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table)


	pVb = [p[7], p[8], p[9], p[10]]  ##11 is sigma V band

	nn = 2.*len(t) + len(tVb)
	
	chiW1 = ( y1 - IRLC_W1_fxdR_ML_point(ptot, t, argW1, RHStable, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phi, ths, nuW1)  )/ np.sqrt( dy1*dy1 + p[6]*p[6]  )

	chiW2 = ( y2 - IRLC_W2_fxdR_ML_point(ptot, t, argW2, RHStable, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phi, ths, nuW2)  )/ np.sqrt( dy2*dy2 + p[6]*p[6]  )

	##USE same sigML as for IR
	#chiVb = 0.5*(yVb - VLC_point(pVb, tVb, argVb, 1.0) )/  np.sqrt( dyVb*dyVb + p[5]*p[5]  )

	#USE diff sigML
	#chiVb = (yVb - VLC_point(pVb, tVb, argVb, 1.0) )/  np.sqrt( dyVb*dyVb + p[11]*p[11]  )
	chiVb = (yVb - VLC_point(pVb, tVb, argVb, 1.0) )/  np.sqrt( dyVb*dyVb )

	print "chiW1^2 = ", sum(chiW1*chiW1)

	# print "chiW1 = ", chiW1

	print "chiW2^2 = ", sum(chiW2*chiW2) 

	# print "chiW2 = ", chiW2

	print "chiV^2 = ", sum(chiVb*chiVb) 

	chi2 = np.sum(chiW1*chiW1) + np.sum(chiW2*chiW2) + np.sum(chiVb*chiVb)

	##USE same sigML as for IR
	#penlty = 0.5 * ( sum( np.log(dy1*dy1 + p[5]*p[5]) ) + sum( np.log(dy2*dy2 + p[5]*p[5]) ) + sum( np.log(dyVb*dyVb + p[5]*p[5]) )  )
	#USE diff sigML
	#penlty = 0.5 * ( np.sum( np.log(dy1*dy1 + p[6]*p[6]) ) + np.sum( np.log(dy2*dy2 + p[6]*p[6]) ) + np.sum( np.log(dyVb*dyVb + p[11]*p[11]) )  )
	penlty = 0.5 * ( np.sum( np.log(dy1*dy1 + p[6]*p[6]) ) + np.sum( np.log(dy2*dy2 + p[6]*p[6]) )  )#+ np.sum( np.log(dyVb*dyVb + p[11]*p[11]) )  )

	negLogLik = chi2/2. + nn/2. * np.log(2.*ma.pi) + penlty


	##LogLik i made negative in likliehood function
	print(-negLogLik)
	return negLogLik



##likliehoods (- Err2 Funcs)
def ln_IR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	return -(IRTDE_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb))
		

def ln_IR_fxdR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	return -(IRTDE_fxdR_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb))
		





##POSTERIORS
def ln_IR_ALL_posterior(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	pIR_prior = [p[0], p[1], p[2], p[3], p[4], p[5], p[6]]
	#pVb = [p[6], p[7], p[8], p[9]]
	#for fitting kk
	pVb_prior = [p[7], p[8], p[9], p[10], p[11]]  
	pVb_prior = [p[7], p[8], p[9], p[10], p[11]]  
	ln_p = ln_IR_ML_prior(pIR_prior) + ln_V_prior(pVb_prior)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p

def ln_IR_fxdR_ALL_posterior(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	pIR_prior = [p[0], p[1], p[2], p[3], p[4], p[5], p[6]]
	#pVb = [p[6], p[7], p[8], p[9]]
	#for fitting kk
	pVb_prior = [p[7], p[8], p[9], p[10], p[11]]  
	ln_p = ln_IR_ML_prior(pIR_prior) + ln_V_prior(pVb_prior)
	if not np.isfinite(ln_p):
		return -np.inf
	
	ln_l = ln_IR_fxdR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p


##for dynasty with ptform() at top of this file
def ln_IR_fxdR_ALL_posterior_dyn(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	ln_l = ln_IR_fxdR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l 

def ln_IR_ALL_posterior_dyn(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	ln_l = ln_IR_ALL_likelihood(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l







#Fmin
def IRTDE_fxdR_FitAll_Err2_fmin(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	pIR_prior = [p[0], p[1], p[2], p[3], p[4], p[5], p[6]]
	#pVb = [p[6], p[7], p[8], p[9]]
	#for fitting kk
	pVb_prior = [p[7], p[8], p[9], p[10], p[11]]  
	ln_p = ln_IR_ML_prior(pIR_prior) + ln_V_prior(pVb_prior)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_fxdR_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p


def IRTDE_sblR_FitAll_Err2_fmin(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb):
	pIR_prior = [p[0], p[1], p[2], p[3], p[4], p[5], p[6]]
	#pVb = [p[6], p[7], p[8], p[9]]
	#for fitting kk
	pVb_prior = [p[7], p[8], p[9], p[10], p[11]]  
	ln_p = ln_IR_ML_prior(pIR_prior) + ln_V_prior(pVb_prior)
	if not np.isfinite(ln_p):
		return np.inf
	
	ln_l = IRTDE_ALL_Err2(p, t, tVb, argW1, argW2, argVb, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, y1, dy1, y2, dy2, yVb, dyVb)
	return ln_l + ln_p





















############################
##Chi^2 Min Liks
############################
#PRIOR
# def ln_IR_prior(params):
# 	etaR, cosTT, sinJJ, nu0, Lav, sigML = params
					
# 	if (sinJJ < -1.0 or sinJJ > 1.0):
# 		return -np.inf

# 	if (cosTT < -1.0 or cosTT > 1.0):
# 		return -np.inf
					
# 	if (etaR <= 0.0 or etaR >=100.0):
# 		return -np.inf

# 	if nu0 <= 0.0:
# 	 	return -np.inf

# 	if (Lav > 100.0 or Lav <= 0.0001):
# 	 	return -np.inf
			
# 	return 0.


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

	if kk < 0.0:
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
	
	chiW1 = 0.5*( y1 - IRLC_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn)  )/ dy1

	chiW2 = 0.5*( y2 - IRLC_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn) )/ dy2

	# print "chiW1^2 = ", sum(chiW1*chiW1)

	# print "chiW1 = ", chiW1

	# print "chiW2^2 = ", sum(chiW2*chiW2) 

	# print "chiW2 = ", chiW2

	chi2 = sum(chiW1*chiW1) + sum(chiW2*chiW2) 
	print(chi2)
	return chi2





def IRTDE_fxdR_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	print "EVAL", p
	
	chiW1 = 0.5*(  y1 - IRLC_fxdR_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn)   )/ dy1

	chiW2 = 0.5*(  y2 - IRLC_fxdR_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn)  )/ dy2 

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
def IR_ML_ptransform(params):
	#etaR, cosTT, sinJJ, FQfac = params
	etaR, cosTT, sinJJ, nu0, kk, Lav, sigML = params
		
	# if (etaR <= 0.0 or etaR >=100.0):
	# 	return -np.inf

	etaRout = 10.**( (etaR*2.0 - 1.0)*2.0)

	# if (sinJJ < -1.0 or sinJJ > 1.0):
	# 	return -np.inf

	sinJJout = (sinJJ*2.0 - 1.0)

	# if (cosTT < -1.0 or cosTT > 1.0):
	# 	return -np.inf

	cosTTout = (cosTT*2.0 - 1.0)

	# if nu0 <= 0.0:
	#  	return -np.inf

	nu0out = 10.**( (nu0*2.0 - 1.0)*5.0) ##from 10^(-5) to 10^(5)

	# if kk < 0.0:
	#  	return -np.inf	

	kkout = 10.0*kk

	# if (Lav > 1000.0 or Lav <= 0.001):
	#  	return -np.

	Lavout = 10.**( (Lav*2.0 - 1.0)*4.0) ##from 10^(-4) to 10^(4)

	# if (sigML < 0.0 or sigML > 0.05):
	# if sigML < 0.0:
	#  	return -np.inf	
			
	x = np.array([etaRout, sinJJout, cosTTout, nu0out, kkout, Lavout, sigML])		

	return x


# def ln_IR_ML_prior(params):
# 	#etaR, cosTT, sinJJ, FQfac = params
# 	etaR, cosTT, sinJJ, nu0, kk, Lav, sigML = params
					
# 	if (sinJJ < -1.0 or sinJJ > 1.0):
# 		return -np.inf

# 	if (cosTT < -1.0 or cosTT > 1.0):
# 		return -np.inf
					
# 	if (etaR <= 0.0 or etaR >=100.0):
# 		return -np.inf

# 	if nu0 <= 0.0:
# 	 	return -np.inf

# 	if kk < 0.0:
# 	 	return -np.inf	

# 	if (Lav > 1000.0 or Lav <= 0.001):
# 	 	return -np.inf

# 	# if (sigML < 0.0 or sigML > 0.05):
# 	if sigML < 0.0:
# 	 	return -np.inf	
			
# 	return 0.

### Point calculations
def IRLC_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn):
	etaR, cosT, sinJJt, nu0, Lav, sigML = p
	FQfac = 100.0
	nu0=nu0*numicron
	nne = kk
	Lav = Lav*10.**(45)
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)

#[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FIR_gal, gam = args
	## Rde doesnt matter here
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_ShTorOptThin_Iso_QuadInt(numn, numx, t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn) + FIR_gal)/FRel)



def IRLC_W1_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus):
	#[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	#FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FIR_gal, gam = args
	FRel, numn, numx, Dst, n0, pp, aeff, etaR, FIR_gal = args

	#etaR, cosT, sinJJt, nu0, Lav, sigML = p
	etaR, cosT, sinJJt, nu0, nne, Lav, t0, tfb, gam = p
	FQfac = 100.0
	t0 = t0*yr2sec
	tfb = tfb*yr2sec
	nu0 = nu0*numicron
	#nne = kk
	Lav = Lav*10.**(45)
	#LVbnd = LVbnd*10.**(45)  # nt used here
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)


	## Rde doesnt matter here
	Rde = etaR*pc2cm
	#Lavg, tfb, n0, Rd, p, thetT, JJ, aeff, nu0, nn, FQfac, t0, etaR, gam = args
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_W1_ShTorOptThin_Iso_QuadInt(t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nus) + FIR_gal)/FRel)


def IRLC_W2_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nus):
	#[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	#FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FIR_gal, gam = args
	FRel, numn, numx, Dst, n0, pp, aeff, etaR, FIR_gal = args

	#etaR, cosT, sinJJt, nu0, Lav, sigML = p
	etaR, cosT, sinJJt, nu0, nne, Lav, t0, tfb, gam = p
	FQfac = 100.0
	t0 = t0*yr2sec
	tfb = tfb*yr2sec
	nu0 = nu0*numicron
	Lav = Lav*10.**(45)
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)



	## Rde doesnt matter here
	Rde = etaR*pc2cm
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_W2_ShTorOptThin_Iso_QuadInt(t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nus) + FIR_gal)/FRel)





def IRLC_fxdR_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths ,nus):
	
	FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, etaR, FIR_gal, gam = args

	Rde, cosT, sinJJt, nu0, kk, Lav, t0, tfb, gam = p
	t0 = t0*yr2sec
	tfb = tfb*yr2sec
	nu0=nu0*numicron
	nne = kk
	Lav = Lav*10.**(45)
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)
	Rde = Rde*pc2cm

	FQfac = 100.0#DEFUNCT


 #[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	
	## etaR doesnt matter here
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_fxdR_ShTorOptThin_Iso_QuadInt(t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn, phis, ths, nus) + FIR_gal)/FRel)







def IRLC_W1_fxdR_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths ,nus):
	#[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
	#FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, etaR, FIR_gal, gam = args
	FRel, numn, numx, Dst, n0, pp, aeff, etaR, FIR_gal = args


	Rde, cosT, sinJJt, nu0, kk, Lav, t0, tfb, gam = p
	t0 = t0*yr2sec
	tfb = tfb*yr2sec
	nu0 = nu0*numicron
	nne = kk
	Lav = Lav*10.**(45)
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)
	Rde = Rde*pc2cm

	FQfac = 100.0#DEFUNCT


	## etaR doesnt matter here
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_W1_fxdR_ShTorOptThin_Iso_QuadInt(t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths ,nus) + FIR_gal)/FRel)




def IRLC_W2_fxdR_ML_point(p, t, args, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths ,nus):
	#[FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
		#FRel, numn, numx, Dst, tfb, n0, pp, aeff, nne, t0, etaR, FIR_gal, gam = args
	FRel, numn, numx, Dst, n0, pp, aeff, etaR, FIR_gal = args

	Rde, cosT, sinJJt, nu0, kk, Lav, t0, tfb, gam = p
	t0 = t0*yr2sec
	tfb = tfb*yr2sec
	nu0=nu0*numicron
	nne = kk
	Lav = Lav*10.**(45)
	thetTst = np.arccos(cosT)
	JJt = np.arcsin(sinJJt)
	Rde = Rde*pc2cm

	FQfac = 100.0#DEFUNCT

 
	## etaR doesnt matter here
	IRargs = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0, etaR, gam]
	return -2.5*np.log10( (F_W2_fxdR_ShTorOptThin_Iso_QuadInt(t, Dst, IRargs, RHStable, Ttable, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths ,nus) + FIR_gal)/FRel)


















##ERR2 FUNCS (-LogLik now)
def IRTDE_ML_Err2(p, t, argW1, argW2, RHStable, Ttable, RHS_mx, RHS_mn, y1, dy1, y2, dy2):
	print "EVAL", p

	nn = 2.*len(t)
	
	chiW1 = 0.5*( y1 - IRLC_ML_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn) )/ np.sqrt( dy1*dy1 + p[5]*p[5]  )

	chiW2 = 0.5*( y2 - IRLC_ML_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn) )/ np.sqrt( dy2*dy2 + p[5]*p[5]  )

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

	chiW1 = 0.5*(  y1 - IRLC_fxdR_ML_point(p, t, argW1, RHStable, Ttable, RHS_mx, RHS_mn)   )/ np.sqrt( dy1*dy1 + p[5]*p[5] )

	chiW2 = 0.5*(  y2 - IRLC_fxdR_ML_point(p, t, argW2, RHStable, Ttable, RHS_mx, RHS_mn)  )/ np.sqrt( dy2*dy2 + p[5]*p[5] )

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
def V_ptransform(params):
	#sinJJ, cosTT, Rin, alpha = params
	L0, t0, tfb, gam, sigML = params
					
	#Transform from unit cube to actual param values

	# if L0 <= 0.001 or L0 > 100.0:  ##in units of 10^45 erg/s  bol
	# 	return -np.inf

	L0out = 10.**( (L0*2.0 - 1.0)*3.0) ##from 10^(-3) to 10^(3)  


	# if t0 < 0.0 or t0 > 2.0: ##dont shift more than 4 years
	# 	return -np.inf

	t0out = t0*2.0
					
	# if tfb < 0.0 or tfb > 5.0:
	# 	return -np.inf

	tfbout = 10.**( (tfb*2.0 - 1.0)*2.0) ##from 10^(-2) to 10^(2)  # same as >0, <100

	# if gam < 0.0 or gam > 5.0:
	# 	return -np.inf

	gamout = gam*4.0

	# if sigML <= 0.0 :
	# 	return -np.inf	
	#keep sigML 0->1

	# if FQfac <= 90.0 :
	# 	return -np.inf

	x = np.array([L0out, t0out, tfbout, gamout, sigML])		

	return x

# def ln_V_prior(params):
# 	#sinJJ, cosTT, Rin, alpha = params
# 	L0, t0, tfb, gam, sigML = params
					
# 	if L0 <= 0.001 or L0 > 100.0:  ##in units of 10^45 erg/s  bol
# 		return -np.inf

# 	if t0 < 0.0 or t0 > 2.0: ##dont shift more than 4 years
# 		return -np.inf
					
# 	if tfb < 0.0 or tfb > 5.0:
# 		return -np.inf

# 	if gam < 0.0 or gam > 5.0:
# 		return -np.inf

# 	if sigML <= 0.0 :
# 		return -np.inf	

# 	# if FQfac <= 90.0 :
# 	# 	return -np.inf
			
# 	return 0.








def VLC_point(p, t, arg, LavIRfit):
	L0, t0, tfb, gam = p
	FQfac = 100.0 ##defunct
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
	chi = 0.5*(y - VLC_point(p, t, arg, 100.0) )/ dy
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






