import cPickle as pickle
import numpy as np
from numpy import *
import scipy as sc
from scipy import *

import emcee
#from emcee.utils import MPIPool

import matplotlib
#matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt


from FluxFuncs_TDEs import *
from emcee_Funcs_TDEs import *

################################
###############################
### OPTIONS
################################
################################
#Trap_Int = False

MaxL = True
Rstrt = False
#RstrtFile = "Restart/Rstrt_sublR_Trap10_MaxLik__src_longerFB_chain.txt"

Pplot = True
plot_solns = False
USE_RSfns = True ##If FALSE remove W1 and W2 RSR fns from plotting

Fit = False
Fit_fmin = True

Fit_ALL = True
Fit_Src = False
Fit_IR = False


Src_BF = False ## doesnt matter if Fit All (all in how set V_prior)
Rfxd = False
No_ReGro = False ##HAVE to change by hand right now in FLuxFn- this only changes names of files



##multiprocessing
NThread = 40
#pool = MPIPool(loadbalance=True)


##Trim beginning of Vband
TrimB = True

zTDE = 0.11






################################
###############################
### IMPORT DATA
################################
################################
print "Importing Data to fit..."

###############################
### OPTICAL DATA
################################
Tt   =  (loadtxt("dat/Vband_Lums_All.txt", usecols=[5], comments="|", delimiter=',')  - 55000)/(1.+zTDE) * day2sec
V_mag  =  loadtxt("dat/Vband_Lums_All.txt", usecols=[1], comments="|", delimiter=',')
V_sig =  loadtxt("dat/Vband_Lums_All.txt", usecols=[2], comments="|", delimiter=',')

V_Mean = mean(V_mag)
#Lum      = (Lum)# - Lum_Mean)  ## for def of mag - mag0

## put in time order
tS   = zip(Tt,V_mag,V_sig)
tS.sort()
TtS  = transpose(tS)
tV_srt   =  np.array(TtS[0])
V_srt    =  np.array(TtS[1])
V_sigsrt =  np.array(TtS[2])



if (TrimB):
	idelE = np.where(tV_srt > 0.0)[0]
	idelB = np.where(tV_srt <= 0.0)[0]
	tV_srt   =  np.delete(tV_srt,idelB)
	V_srt    =  np.delete(V_srt,idelB)
	V_sigsrt =  np.delete(V_sigsrt,idelB)

	tV_srtE   =  np.delete(tV_srt,idelE)
	V_srtE    =  np.delete(V_srt,idelE)
	V_sigsrtE =  np.delete(V_sigsrt,idelE)

V_Gal = 17.6#mean(V_srtE)
#V_srt = V_srt - V_Gal


# magVnet = -2.5*np.log10(10.**(-V_srtE/2.5) - 10.**(-V_Gal/2.5))
# V_srt = magVnet




###### get average value for each cluster of data points in time
iseg = []
iseg.append(-1)
for i in range(len(tV_srt)-1):
	if (abs((tV_srt[i] - tV_srt[i+1])) > 100*day2sec ):
		iseg.append(i)
iseg.append(len(tV_srt)-1)

tV_avg = []
V_avg = []
V_avsg = []

#VarMn_test_W1 = []
#VarMn_test_W2 = []


for i in range(0 , len(iseg)-1):
	tV_avg.append(np.mean(tV_srt[iseg[i]+1:iseg[i+1]+1]))

	V_avg.append(np.mean(V_srt[iseg[i]+1:iseg[i+1]+1]))
	

	
	Nseg = len(V_sigsrt[iseg[i]+1:iseg[i+1]]) + 1


	V_avsg.append(np.sqrt(sum( (V_sigsrt[iseg[i]+1:iseg[i+1]])**2 )/Nseg  ))


	# ## see if variance of means or mean of variances is larger
	# VarMn_test_W1.append(np.mean(W1_sig[iseg[i]+1:iseg[i+1]+1])/np.std(W1_mag[iseg[i]+1:iseg[i+1]+1]))
	# VarMn_test_W2.append(np.mean(W2_sig[iseg[i]+1:iseg[i+1]+1])/np.std(W2_mag[iseg[i]+1:iseg[i+1]+1]))


tV_avg = np.array(tV_avg)
V_avg = np.array(V_avg)
V_avsg = np.array(V_avsg)









###############################
### IR DATA
################################
### All IR DATA already cleaned and reprinted in "WISE_plotGen.py"
t_MJD = (np.genfromtxt("dat/IR_Lums_All.txt",usecols=1, comments="|"))/(1.+zTDE) * day2sec## put in rest frame

W1_mag = np.genfromtxt("dat/IR_Lums_All.txt",usecols=2, comments="|")
W2_mag = np.genfromtxt("dat/IR_Lums_All.txt",usecols=3, comments="|")

W1_sig = np.genfromtxt("dat/IR_Lums_All.txt",usecols=4, comments="|")
W2_sig = np.genfromtxt("dat/IR_Lums_All.txt",usecols=5, comments="|")

##binned
t_avg = (np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=0, comments="|"))/(1.+zTDE) * day2sec## put in rest frame

W1_avg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=1, comments="|")
W2_avg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=2, comments="|")

W1_avsg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=3, comments="|")
W2_avsg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=4, comments="|")


t_avg = np.array(t_avg)
W1_avg = np.array(W1_avg)
W2_avg = np.array(W2_avg)
W1_avsg = np.array(W1_avsg)
W2_avsg = np.array(W2_avsg)



####TEMP: FORCE QUESCIENT POINTS TO HAVE SMALL ERRORS!
W1_avsg[0] *= 0.1
W1_avsg[1] *= 0.1
W2_avsg[0] *= 0.1
W2_avsg[1] *= 0.1

W1_avsg[8] *= 5.0
W1_avsg[7] *= 5.0
W1_avsg[6] *= 2.0

W2_avsg[8] *= 8.0
W2_avsg[7] *= 5.0
W2_avsg[6] *= 3.0



##TEMP CHECK
# W1_avsg[0] = W1_avsg[0]*10.
# W1_avsg[1] = W1_avsg[1]*10.

# W2_avsg[0] = W2_avsg[0]*10.
# W2_avsg[1] = W2_avsg[1]*10.

## subtract off galaxy background form first two preflare epochs
# magW1net = -2.5*np.log10(10.**(-W1_avg/2.5) - 10.**(-12.9/2.5))
# magW2net = -2.5*np.log10(10.**(-W2_avg/2.5) - 10.**(-11.26/2.5))
# W1_avg = magW1net
# W2_avg = magW2net

##########################
## SOURCE PROPERTIES
##########################
#BEST FITS:
if (Src_BF):
	Lav = 10.**45 #erg/s
	LVbnd = 0.18523108#*10.**45
	Dst = 0.513*10**9*pc2cm
	t0 =  0.75866946  * yr2sec
	tfb = 0.02625009 * yr2sec
	gam = 1.12948904
#^mcmc beest fit [ 0.18523108,  0.75866946,  0.02625009,  1.12948904]
else:
	Lav = 10.**45 #erg/s
	LVbnd = 0.0387673
	Dst = 0.513*10**9*pc2cm
	t0 =  0.49686742  * yr2sec
	tfb = 0.50811559 * yr2sec
	gam = 1.74163261
#^restricting tfb min to be 0.5 [ 0.0387673 ,  0.49686742,  0.50811559,  1.74163261]

## TEST VALUES
### DUST stuff
nne = 0.0#1.6  ##emission/absortption efficeincy power law  --  Draine or 1.6 
nu0 = numicron/3.0  ## emission/absortption efficeincy cutoff freq
Rde = 0.0001*pc2cm #* np.sqrt(Lav/10.**46) * pc2cm     ## inner edge of dust
etaR = 15.0 #factor multiplying sublmation radius when using Rde(T) = subl front
thetTst = 0.01*np.pi/4.  ## opening angle for dust torus (0 is sphere, pi/2 is ring)
JJt = 0.01*np.pi/4.       ## inclination angle of dust torus (0 is edge on pi/2 is face on)
aeff = (c/nu0)/(2.*np.pi) #set dust grain size by cutoff frequency


###OPT THIN DUST
##For dust which is optically thin to opt/UV ##IGNORE FOR NOW
Rrout = 1.0*Rde ## oter edge of dust
pp = 2.0 	    ## dust density power law
n10 = 1.0#/(np.pi*Rde*aeff*aeff) * (pp-1.)
nfac = 10.0 
n0 = nfac*n10 ##6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))
##########################
## SOURCE PROPERTIES
##########################

# ### WISE BAND + Observational STUFF
## edges of WISE bands
####IMPORT relative spectral response functions
nu1_RSR = c/(np.genfromtxt("dat/BandW1.txt",usecols=0, comments="#")*10.**(-8))
W1_RSR = np.genfromtxt("dat/BandW1.txt",usecols=1, comments="#")

nu2_RSR = c/(np.genfromtxt("dat/BandW2.txt",usecols=0, comments="#")*10.**(-8))
W2_RSR = np.genfromtxt("dat/BandW2.txt",usecols=1, comments="#")

W1RSR_intrp = sc.interpolate.interp1d(nu1_RSR, W1_RSR)
W2RSR_intrp = sc.interpolate.interp1d(nu2_RSR, W2_RSR)

###TEST Interp
#nus1 = linspace(nu1_RSR[0], nu1_RSR[len(nu1_RSR)-1], 100)
#nus2 = linspace(nu2_RSR[0], nu2_RSR[len(nu2_RSR)-1], 100)
# plt.figure()
# plt.scatter(nu1_RSR, W1_RSR, color='orange')
# plt.scatter(nu2_RSR, W2_RSR, color='red')
# plt.plot(nus1, W1RSR_intrp(nus1))
# plt.plot(nus2, W2RSR_intrp(nus2))
# plt.savefig('RSR_interps.png')



###
W1mn_eff = numicron/2.8
W1mx_eff = numicron/4.0
W2mn_eff = numicron/3.9
W2mx_eff = numicron/5.3

W1mn = nu1_RSR[2]#numicron/2.8
W1mx = nu1_RSR[len(nu1_RSR)-3]# numicron/4.0
W2mn = nu2_RSR[2]#numicron/3.9
W2mx = nu2_RSR[len(nu2_RSR)-3]#numicron/5.3

Vmx = c/(5.45*10**(-5) - 0.44*10**(-5))
Vmn = c/(5.45*10**(-5) + 0.44*10**(-5))


### zero point fluxes - check these...
nuVbnd = c/(5.45*10**(-5))
FVbndRel = 3.636*10**(-20)*nuVbnd 
#FW1Rel = 3.09540*10**(-21)*8.8560*10**(13)#(W1mn + W1mx)/2
#FW2Rel = 1.71787*10**(-21)*6.4451*10**(13)#(W2mn + W2mx)/2

#FW1Rel = 3.09540*10**(-21)*(W1mx_eff-W1mn_eff)#(W1mn + W1mx)/2 ##IN VEGA and mult by dv that we integrate over for tophat(
#FW2Rel = 1.71787*10**(-21)*(W2mx_eff-W2mn_eff)#(W2mn + W2mx)/2
Nnrm = 10000.0
nus1 = np.linspace(W1mn, W1mx, Nnrm)
nus2 = np.linspace(W2mn, W2mx, Nnrm)
FW1Rel = 3.09540*10**(-21)*(W1mx-W1mn)/(2.*Nnrm) * (2.0 * np.sum([W1RSR_intrp(nu) for nu in nus1]) - W1RSR_intrp(nus1[0]) -   W1RSR_intrp(nus1[Ntrap_nu-1]) )
FW2Rel = 1.71787*10**(-21)*(W2mx-W2mn)/(2.*Nnrm) * (2.0 * np.sum([W2RSR_intrp(nu) for nu in nus2]) - W2RSR_intrp(nus2[0]) -   W2RSR_intrp(nus2[Ntrap_nu-1]) )



##TABULATE T's and RHSs
print "Creating look up tables"
NT = 1800
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 1800., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)
RHS_mx = RHS_table[len(RHS_table)-1]
RHS_mn = RHS_table[0]


plt.figure()
plt.scatter(T_table, RHS_table)
plt.savefig(RHS_Table)



#Temp_Num(t, r, Lavg, tfb, t0, gam, FQfac, nu0, nn)
#Td_intrp = sc.interpolate.interp1d(tts, Temp_Num(tts, Rde, Lav, tfb, t0, gam, FQfac, nu0, nn))



## measured background flux to add to TDE LC
FV_gal = 10.**(-17.6/2.5) * FVbndRel
FW1_gal = 10.**(-12.9/2.5) * FW1Rel
FW2_gal = 10.**(-11.26/2.5) * FW2Rel


###CALCUALTE HOW MUCH VBAND contributes to IR - does FV_gal go atend here?
#argVplus = [FVbndRel, Vmn, Vmx, Dst, Lav, tfb, n0, pp, aeff, nne, t0, Rde, 0.0]
argVplus = [FVbndRel, Vmn, Vmx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, 0.0, gam]

argtst = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, 100., t0, etaR, gam]


# argW1 = [FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal, gam]
# argW2 = [FW2Rel, W2mn, W2mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW2_gal, gam]
argW1 = [FW1Rel, W1mn, W1mx, Dst, n0, pp, aeff, Rde, FW1_gal]
argW2 = [FW2Rel, W2mn, W2mx, Dst, n0, pp, aeff, Rde, FW2_gal]
Varg  = [Dst, FVbndRel, FV_gal]


phis = np.linspace(0.0,2.*ma.pi, Ntrap_ph)
ths = np.linspace(0.0, ma.pi, Ntrap_th)
nuW1s = np.linspace(W1mn, W1mx, Ntrap_nu)
nuW2s = np.linspace(W2mn, W2mx, Ntrap_nu)

nuW1 = np.meshgrid(phis, ths, nuW1s) 
nuW2 = np.meshgrid(phis, ths, nuW2s) 

sigML = 0.1


# class Args2Pass:
# 	RHS_mx = RHS_mx
# 	RHS_mn = RHS_mn
# 	RHS_table = RHS_table
# 	T_table   = T_table
# 	argVplus = argVplus
# 	argtst = argtst 
# 	argW1 = argW1
# 	argW2 = argW2
# 	Varg = Varg






if (Fit_fmin):
	from scipy import optimize
	from scipy.optimize import fmin


	if (Fit_ALL):
		if (Rfxd):
				Shell_File = "_FminML_FitALLIRandOPt_FxdR_Trap%g_" %Ntrap_nu
				param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", "k", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$", r"$\sigma^{\rm{Vbnd}_{\rm{ML}}$"]
				#p0 =[  0.956647,   0.969297,   1.0,  0.302956, 3.33748, 0.107118, 0.0374935, 0.47906, 0.492371, 1.67332]
				#p0 = [1.40165688,  0.25494203,  0.94144883,  0.28295825,  10.0,  0.0875567, LVbnd, t0/yr2sec, tfb/yr2sec, gam ]
				#p0 = [1.00165688,  0.25494203,  0.94144883,  0.28295825,  8.0,  0.0875567, LVbnd, t0/yr2sec, tfb/yr2sec, gam ]
				#p0 = [1.23793989,  0.22661969,  0.75340754,  0.36330853,  7.07261987, 0.04719663,  0.04075957,  0.55869568,  0.6239391 ,  1.5733816 ] 
				#p0 = [1.16647024,  0.22769808,  0.75363487,  0.36921657,  6.97276134, 0.04658234,  0.04162513,  0.57117547,  0.62769134,  1.65204143 ]
				#p0 = [9.23726538e-01,   9.02902839e-02,  -1.39823086e-01, 1.52954915e-01,   6.66873575e+00,   1.31067343e-01, 6.61179768e-02,   8.70825658e-05,   1.24761379e+00, 2.24192681e+00, 0.01]
			
				## WITH extra sigmaML for Vband
				#p0 = [  8.99271639e-01,   9.22691018e-02,  -1.42585549e-01, 1.46374543e-01, 1.6,  6.79000508e+00,   1.32665069e-01, 6.66665189e-02,   9.00623387e-05,   1.26117162e+00, 2.29801380e+00,   1.02734196e-02]
				p0 = [ 9.38084792e-01,   9.18561421e-02, -1.31365501e-01,   1.31487333e-01, 2.52415700e+00,   4.49555437e+00,   1.64801865e-01,   2.19669805e-02, 7.49921946e-05,   7.63580401e-01,   1.08668587e+00,    2.29383233e-02 ]

				IR_p0 = np.array(p0)
				popt  = sc.optimize.fmin(IRTDE_fxdR_FitAll_Err2_fmin,  IR_p0, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg), full_output=1, disp=False, ftol=0.1)[0]
		else:
				if (No_ReGro):
					Shell_File = "_FminML_FitALLIRandOPt_sblR_NOREGROW_Trap%g_" %Ntrap_nu
				else:
					Shell_File = "_FminML_FitALLIRandOPt_sblR_Trap%g_" %Ntrap_nu
				param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", "k", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$", r"$\sigma^{\rm{Vbnd}_{\rm{ML}}$"]
				#p0 = [  14.6678,   0.980941,   1.0, 0.295263,   5.19492, 0.0156566, 0.0421324, 0.483296, 0.291635, 2.25087]
				#p0 = [ 18.0,   0.70941,   0.0, 0.295263,   2.0, 0.0156566, 0.0421324, 0.483296, 0.291635, 2.25087]
				#p0 = [18.0,   7.08844123e-01,  -1.25866437e-03, 5.18330219e-01,  8.67447142e+00,   2.05544796e-02, LVbnd, t0/yr2sec, tfb/yr2sec, gam]
				#p0 = [ 1.59415002e+01,   7.18785035e-01,  -9.62103056e-04, 6.04657602e-01,   1.20797977e+01,   2.66838527e-02, 5.03637857e-02,   3.51993058e-01,   5.51367583e-01, 1.27481570e+00 ]
				#p0 = [ 1.54100504e+01,   3.82739483e-01,   4.66174470e-04, 1.93773202e-01,   5.75805590e+00,   1.74165532e-01, 2.51103284e-03,   2.37980357e-01,   4.64101485e-01, 3.49929669e+00 ]
				#p0 = [ 1.54100504e+01,   3.82739483e-01,   4.66174470e-04, 1.93773202e-01,   5.75805590e+00,   1.74165532e-01, LVbnd, t0/yr2sec, tfb/yr2sec, gam, 0.01]
				
				### with extra sigmaML for Vband
				#p0 = [  1.54555898e+01,   3.81640899e-01,   4.66995006e-04, 1.94302693e-01,   5.78479853e+00,   1.74874406e-01, 3.88715034e-02,   4.98321895e-01,   5.07851677e-01, 1.74972141e+00,   1.00254982e-02]
				#p0 = [  1.54555898e+01,   3.81640899e-01,   4.66995006e-04, 1.94302693e-01,   1.6, 5.78479853e+00,   1.74874406e-01, 3.88715034e-02,   4.98321895e-01,   5.07851677e-01, 1.74972141e+00,   1.00254982e-02]
				#p0  = [  1.53352529e+01,   3.68507850e-01,   5.10341847e-04, 2.01771863e-01,   1.72056978e+00,   7.07650311e+00, 1.52872603e-01,   4.34919585e-02,   2.89131689e-01, 5.50664403e-01,   1.89425164e+00,   9.65336541e-03]
				p0  = [  1.50648220e+01,   3.50311977e-01,   5.04816357e-04, 1.98229977e-01,   1.00099855e+00,   7.68582797e+00, 1.51661512e-01,  5.00930035e-02,  1.1, 0.06,   1.1,   1.03412462e-02]


				IR_p0 = np.array(p0)
				popt  = sc.optimize.fmin(IRTDE_sblR_FitAll_Err2_fmin,  IR_p0, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg), full_output=1, disp=False, ftol=0.1)[0]

				





	else:

		if (MaxL):
		#arg1 = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
			if (Rfxd):
				Shell_File = "_FitALLIR_FxdR_Trap%g_fminML_" %Ntrap_nu
				param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$"]
				if (Src_BF):
					#p0IR = [0.8, np.cos(thetTst), np.sin(JJt), nu0/numicron, Lav/10.**45]
					p0IR =[  8.16625053e-01,   9.83459365e-01,   7.80853708e-03,  3.60872083e-06,   3.44648331e+00, sigML]
				else:
					#longer fallback
					#1024 MCMC fit
					p0IR = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476, sigML]
					#p0IR = [  9.83366786e-01,   7.80424532e-01,   9.05889101e-03, 1.13584010e-05,   1.79596264e+00]
				IR_p0 = np.array(p0IR)
				#p0IR = [0.229, 0.862, 0.05, 0.147, 18.4891]
				popt  = sc.optimize.fmin(IRTDE_fxdR_ML_Err2_fmin,  IR_p0, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False, ftol=0.0001)[0]

				

			else:
				Shell_File = "_FitALLIR_sublR_Trap%g_fminML_" %Ntrap_nu
				param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$"]
				#p0IR = [etaR, np.cos(thetTst), np.sin(JJt), nu0/numicron, Lav/10.**45]
				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					p0IR = [19.371,   0.753721,  0.00895658,  0.0106707,  1.88784, sigML]
				else:
					Shell_File = Shell_File + "_src_longerFB_"
					#longer fallback
					p0IR = [  2.01213715e+01,   9.94873454e-01,   6.76740228e-03, 1.71450267e-06,   1.16555735e+00,   1.99141744e-01]
				IR_p0 = np.array(p0IR)
				popt  = sc.optimize.fmin(IRTDE_Err2_ML_fmin,  IR_p0, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False, ftol=0.0001)[0]

		else:
			if (Rfxd):
				Shell_File = "_FitALLIR_FxdR_Trap%g_fminNOML_" %Ntrap_nu
				param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", r"$L_{45}$"]
				if (Src_BF):
					p0IR =[  8.16625053e-01,   9.83459365e-01,   7.80853708e-03,  3.60872083e-06,   3.44648331e+00]
				else:
					#longer fallback
					p0IR = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476]
				IR_p0 = np.array(p0IR)
				popt  = sc.optimize.fmin(IRTDE_fxdR_Err2_fmin,  IR_p0, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False, ftol=0.0001)[0]
			
			else:
				Shell_File = "_FitALLIR_sublR_Trap%g_fmin_" %Ntrap_nu
				param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", r"$L_{45}$"]
				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					p0IR = [19.371,   0.753721,  0.00895658,  0.0106707,  1.88784]
				else:
					Shell_File = Shell_File + "_src_longerFB_"
					#longer fallback
					p0IR = [  2.01213715e+01,   9.94873454e-01,   6.76740228e-03, 1.71450267e-06,   1.16555735e+00]
				IR_p0 = np.array(p0IR)
				popt  = sc.optimize.fmin(IRTDE_Err2_fmin,  IR_p0, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False, ftol=0.0001)[0]






	if (Fit_ALL):
		All_p_opt = popt
	else:
		IR_p_opt = popt

	filename = "emcee_Results/TDE_results_"+Shell_File+".txt"
	print "Printing Results"
	target = open(filename, 'w')
	target.truncate()

	for i in range(len(popt)):
		target.write(param_names[i]+" = %g" %popt[i])
		target.write("\n")

	target.close()















if (Fit):


	if(Fit_ALL):
		##ALways doing MaxL
		if (Rfxd):
				Shell_File = "_FitIRandOpt_FxdR_Trap%g_MaxLik_" %Ntrap_nu
				param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$"]
				p0 = [ 9.38084792e-01,   9.18561421e-02, -1.31365501e-01,   1.31487333e-01, 2.52415700e+00,   4.49555437e+00,   1.64801865e-01,   2.19669805e-02, 7.49921946e-05,   7.63580401e-01,   1.08668587e+00,    2.29383233e-02 ]

				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					p0 = 
				else:
					Shell_File = Shell_File + "_src_longerFB_"
					##longer fall backs
					#1024 MCMC
					p0 = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476, sigML, LVbnd, t0/yr2sec, tfb/yr2sec, gam]
					perr = [0.0037, -0.0607, 0.3757, -0.1978, -0.0736]
					merr = [0.0003, 0.1323, -0.1746, 0.2659, 0.1329]
					#p0IR = [  9.83366786e-01,   7.80424532e-01,   9.05889101e-03, 1.13584010e-05,   1.79596264e+00]
				ndim = len(p0)
				nwalkers = ndim*4

				All_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_fxdR_ALL_posterior, threads=NThread, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg))


		else:
				Shell_File = "_sublR_Trap%g_MaxLik_" %Ntrap_nu
				param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$"]
				#p0IR = [etaR, np.cos(thetTst), np.sin(JJt), nu0/numicron, Lav/10.**45]
				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					#best fit
					#p0IR = [19.371,   0.753721,  0.00895658,  0.0106707,  1.88784]
					#fmin
					p0 = [  14.6678,   0.980941,   1.0, 0.295263,   5.19492, 0.0156566, 0.0421324, 0.483296, 0.291635, 2.25087]
					#p0 = [  1.62162628e+01,   9.94984173e-01,   9.07876972e-03, 0.3,   2.70076201e+00, sigML, LVbnd, t0/yr2sec, tfb/yr2sec, gam]
				else:
					#longer fallback
					Shell_File = Shell_File + "_src_longerFB_"
					p0 = [  1.92818999e+01,   8.74318091e-01,   1.07231518e-02, 0.3,   1.65353712e+00, sigML, LVbnd, t0/yr2sec, tfb/yr2sec, gam]
				ndim = len(p0)
				nwalkers = ndim*4


				All_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_ALL_posterior, threads=NThread, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg))

		All_p0 = np.array(p0)

		All_walker_p0 = np.random.normal(All_p0, np.abs(All_p0)*1E-4, size=(nwalkers, ndim))

					
		clen = 512
		#for ():
		#run as iterable
		#acor function
		All_pos,_,_ = All_sampler.run_mcmc(All_walker_p0, clen)
		#manipulte (replace) fidn outliers 2std form median in beginning


		print "SAVING THE PICKLE mmmmm"
		with open("emcee_data/Pickles/TDE_All_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((All_sampler.chain, All_sampler.lnprobability), f1)


				


		### OPEN OUTPUT DATA
		with open("emcee_data/Pickles/TDE_All_%iwalkers.pickle" %clen) as f1:
			All_chain,All_lnprobs = pickle.load(f1)



		#V_acor = V_sampler.acor

		All_flatchain = np.vstack(All_chain[:,clen/4:])
		All_flatlnprobs = np.vstack(All_lnprobs[:,clen/4:])
				

		All_p_opt = All_flatchain[All_flatlnprobs.argmax()]		



		##record final state for restart
		f_rstrt = open("Restart/Rstrt"+Shell_File+"chain.txt", "w")
		f_rstrt.close()

		for result in All_sampler.sample(All_pos, iterations=1, storechain=False):
		    position = result[0]
		    f_rstrt = open("Restart/Rstrt"+Shell_File+"chain.txt", "a")
		    for k in range(position.shape[0]):
		    	#print k
		    	f_rstrt.write("%i  %g %g %g %g" %(k, position[k][0], position[k][1], position[k][2], position[k][3]))
		    	f_rstrt.write("\n")
		        #f_rstrt.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
		    f_rstrt.close()



		print "ANALYSING MCMC (!)..."
		with open("emcee_data/Pickles/TDE_All_%iwalkers.pickle" %clen) as f1:
			All_chain,All_lnprobs = pickle.load(f1)



		##PLOT dem WALKERS
		for k in range(All_chain.shape[2]):
			plt.figure()
			#plt.figure()
			for i in range(All_chain.shape[0]):
				plt.plot(All_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(param_names[k])
				plt.xlabel('steps')
				plt.tight_layout()
			plt.savefig('emcee_data/src_'+Shell_File+'_TDE_%s_%iwalkers.png' %(param_names[k],clen))
			plt.clf()




		###CORNER PLOT	
		All_flatchain = np.vstack(All_chain[:,clen/4:])
		All_flatlnprobs = np.vstack(All_lnprobs[:,clen/4:])



				

		#import triangle
		import corner as triangle
		All_fig = triangle.corner(All_flatchain,  labels=param_names, quantiles=[0.15, 0.5, 0.85],show_titles=True, title_kwargs={"fontsize": 14},label_kwargs={"fontsize": 18})			
		All_fig.savefig('emcee_data/src_'+Shell_File+'_TDE_Corner_Plot_%iwalkers.png' %clen)




		## Do some stats on the walkers
		from scipy.stats import scoreatpercentile as scoretpercentile



		## max posterior + percentiles
		All_MAP_vals = All_flatchain[All_flatlnprobs.argmax()]
		All_perc = scoretpercentile(All_flatchain, [15,85], axis=0)




		filename = "emcee_Results/TDE_results"+Shell_File+"%iwalkers.txt" %clen
		print "Printing Results"
		target = open(filename, 'w')
		target.truncate()

		All_diff_minus = np.zeros(len(param_names))
		All_diff_plus  = np.zeros(len(param_names))
		for i,name in enumerate(param_names):
			All_diff_minus[i] = All_MAP_vals[i] - All_perc[0,i]
			All_diff_plus[i] = All_perc[1,i] - All_MAP_vals[i]
			target.write("TDE_All: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(All_MAP_vals[i], All_diff_plus[i], All_diff_minus[i], name=name))
			target.write("\n")





		All_mxprbs = zeros(nwalkers)

					

		for i in range(nwalkers):
			All_mxprbs[i] = max(All_lnprobs[i])


		chi2_pdf_All = -max(All_mxprbs)/(len(tV_avg) + 2*len(t_avg) - len(param_names) - 1)/2 ## add 2 not in liklihood func
		target.write("\n")		
		target.write("Max Lik =  %04g" %max(All_mxprbs))
		target.write("\n")
		target.write("\n")		
		target.write("Shell TDE Vband source fit, reduced chi2 =  %04g" %chi2_pdf_All)
		target.write("\n")
		target.close()


















































	if (Fit_Src):
		##########################
		## FIT SOURCE PROPERTIES
		##########################

		Shell_File = "_Fit_Vband_"
		param_names = [r"$L_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$"]
		p0V = [LVbnd, t0/yr2sec, tfb/yr2sec, gam]
		ndim = len(p0V)
		nwalkers = ndim*4
		#V_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_V_posterior, threads=NThread,args=(tV_srt, Varg, V_srt, V_sigsrt))
		V_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_V_posterior, threads=NThread, args=(tV_avg, Varg, V_avg, V_avsg))




		V_p0 = np.array(p0V)

		V_walker_p0 = np.random.normal(V_p0, np.abs(V_p0)*1E-3, size=(nwalkers, ndim))

					
		clen = 512
		#for ():
		#run as iterable
		#acor function
		V_pos,_,_ = V_sampler.run_mcmc(V_walker_p0, clen)
		#manipulte (replace) fidn outliers 2std form median in beginning


		print "SAVING THE PICKLE mmmmm"
		with open("emcee_data/Pickles/TDE_V_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((V_sampler.chain, V_sampler.lnprobability), f1)


				


		### OPEN OUTPUT DATA
		with open("emcee_data/Pickles/TDE_V_%iwalkers.pickle" %clen) as f1:
			V_chain,V_lnprobs = pickle.load(f1)



		#V_acor = V_sampler.acor

		V_flatchain = np.vstack(V_chain[:,clen/4:])
		V_flatlnprobs = np.vstack(V_lnprobs[:,clen/4:])
				

		V_p_opt = V_flatchain[V_flatlnprobs.argmax()]		



		##record final state for restart
		f_rstrt = open("Restart/Rstrt"+Shell_File+"chain.txt", "w")
		f_rstrt.close()

		for result in V_sampler.sample(V_pos, iterations=1, storechain=False):
		    position = result[0]
		    f_rstrt = open("Restart/Rstrt"+Shell_File+"chain.txt", "a")
		    for k in range(position.shape[0]):
		    	#print k
		    	f_rstrt.write("%i  %g %g %g %g" %(k, position[k][0], position[k][1], position[k][2], position[k][3]))
		    	f_rstrt.write("\n")
		        #f_rstrt.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
		    f_rstrt.close()



		print "ANALYSING MCMC (!)..."
		with open("emcee_data/Pickles/TDE_V_%iwalkers.pickle" %clen) as f1:
			V_chain,V_lnprobs = pickle.load(f1)



		##PLOT dem WALKERS
		for k in range(V_chain.shape[2]):
			plt.figure()
			#plt.figure()
			for i in range(V_chain.shape[0]):
				plt.plot(V_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(param_names[k])
				plt.xlabel('steps')
				plt.tight_layout()
			plt.savefig('emcee_data/src_'+Shell_File+'_TDE_%s_%iwalkers.png' %(param_names[k],clen))
			plt.clf()




		###CORNER PLOT	
		V_flatchain = np.vstack(V_chain[:,clen/4:])
		V_flatlnprobs = np.vstack(V_lnprobs[:,clen/4:])



				

		#import triangle
		import corner as triangle
		V_fig = triangle.corner(V_flatchain,  labels=param_names, quantiles=[0.15, 0.5, 0.85],show_titles=True, title_kwargs={"fontsize": 14},label_kwargs={"fontsize": 18})			
		V_fig.savefig('emcee_data/src_'+Shell_File+'_TDE_Corner_Plot_%iwalkers.png' %clen)




		## Do some stats on the walkers
		from scipy.stats import scoreatpercentile as scoretpercentile



		## max posterior + percentiles
		V_MAP_vals = V_flatchain[V_flatlnprobs.argmax()]
		V_perc = scoretpercentile(V_flatchain, [15,85], axis=0)




		filename = "emcee_Results/TDE_results_"+Shell_File+"%iwalkers.txt" %clen
		print "Printing Results"
		target = open(filename, 'w')
		target.truncate()


		for i,name in enumerate(param_names):
			V_diff_minus = V_MAP_vals[i] - V_perc[0,i]
			V_diff_plus = V_perc[1,i] - V_MAP_vals[i]
			target.write("TDE_V: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(V_MAP_vals[i], V_diff_plus, V_diff_minus, name=name))
			target.write("\n")





		V_mxprbs = zeros(nwalkers)

					

		for i in range(nwalkers):
			V_mxprbs[i] = max(V_lnprobs[i])


		chi2_pdf_V = -max(V_mxprbs)/(len(tV_srt) - len(param_names) - 1)/2 ## add 2 not in liklihood func
		target.write("\n")		
		target.write("Shell TDE Vband source fit, reduced chi2 =  %04g" %chi2_pdf_V)
		target.write("\n")
		target.close()




























	if (Fit_IR): 
		##########################
		## FIT DUST PROPERTIES
		##########################
		##fit for Lav, tfb with optical data
		#arg1 = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]

		if (MaxL):
			if (Rfxd):
				Shell_File = "_FxdR_Trap%g_MaxLik_" %Ntrap_nu
				param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$"]
				#p0IR = [0.8, np.cos(thetTst), np.sin(JJt), nu0/numicron, Lav/10.**45]
				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					#fmin best
					#p0IR = [0.8149, 1.0, 0.00886, 1.105e-5, 2.13 ]
					p0IR =[  8.16625053e-01,   9.83459365e-01,   7.80853708e-03,  3.60872083e-06,   3.44648331e+00, sigML]
				else:
					Shell_File = Shell_File + "_src_longerFB_"
					##longer fall back
					#1024 MCMC
					p0IR = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476, sigML]
					perr = [0.0037, -0.0607, 0.3757, -0.1978, -0.0736]
					merr = [0.0003, 0.1323, -0.1746, 0.2659, 0.1329]
					#p0IR = [  9.83366786e-01,   7.80424532e-01,   9.05889101e-03, 1.13584010e-05,   1.79596264e+00]
				ndim = len(p0IR)
				nwalkers = ndim*5

				IR_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_fxdR_ML_posterior, threads=NThread, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg))


			else:
				Shell_File = "_sublR_Trap%g_MaxLik_" %Ntrap_nu
				param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$"]
				#p0IR = [etaR, np.cos(thetTst), np.sin(JJt), nu0/numicron, Lav/10.**45]
				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					#best fit
					#p0IR = [19.371,   0.753721,  0.00895658,  0.0106707,  1.88784]
					p0IR = [  1.62162628e+01,   9.94984173e-01,   9.07876972e-03, 9.44932700e-03,   2.70076201e+00, sigML]
				else:
					#longer fallback
					Shell_File = Shell_File + "_src_longerFB_"
					p0IR = [  1.92818999e+01,   8.74318091e-01,   1.07231518e-02, 9.81932600e-03,   1.65353712e+00, sigML]
				ndim = len(p0IR)
				nwalkers = ndim*5


				IR_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_ML_posterior, threads=NThread, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg))



		else:
			if (Rfxd):
				Shell_File = "_FxdR_Trap%g_NoML_" %Ntrap_nu
				param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", r"$L_{45}$"]
				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					p0IR =[  8.16625053e-01,   9.83459365e-01,   7.80853708e-03,  3.60872083e-06,   3.44648331e+00]
				else:
					Shell_File = Shell_File + "_src_longerFB_"
					p0IR = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476]
				ndim = len(p0IR)
				nwalkers = ndim*5

				IR_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_fxdR_posterior, threads=NThread, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg))


			else:
				Shell_File = "_sublR_Trap%g_NoML_" %Ntrap_nu
				param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", r"$L_{45}$"]
				if (Src_BF):
					Shell_File = Shell_File + "_src_BF_"
					p0IR = [  1.62162628e+01,   9.94984173e-01,   9.07876972e-03, 9.44932700e-03,   2.70076201e+00]
				else:
					#longer fallback
					Shell_File = Shell_File + "_src_longerFB_"
					p0IR = [  1.92818999e+01,   8.74318091e-01,   1.07231518e-02, 9.81932600e-03,   1.65353712e+00]
				ndim = len(p0IR)
				nwalkers = ndim*5



				IR_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_posterior, threads=NThread, args=(t_avg, argW1, argW2, RHS_table, T_table, RHS_mx, RHS_mn, W1_avg, W1_avsg, W2_avg, W2_avsg))




		if (Rstrt):
			##define walker init from file
			RstrtFile = "Restart/Rstrt"+Shell_File+"chain.txt"
			print "Restarting from File "+RstrtFile 
			in0 = np.zeros([ndim, nwalkers])
			for i in range(0, len(param_names)):
				in0[i] = np.array(np.genfromtxt(RstrtFile, usecols=i+1, comments='$'))

			IR_walker_p0 = np.transpose(in0)
		else:
			IR_p0 = np.array(p0IR)

			IR_walker_p0 = np.random.normal(IR_p0, np.abs(IR_p0)*1E-4, size=(nwalkers, ndim))

					
		clen = 1024#4048
		IR_pos,_,_ = IR_sampler.run_mcmc(IR_walker_p0 , clen)



		print "SAVING THE PICKLE mmmmm"
		with open("emcee_data/Pickles/TDE_IR_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((IR_sampler.chain, IR_sampler.lnprobability), f1)


				
		### SAAVE END STATE to RESTART (IR_pos)
		f_rstrt = open("Restart/Rstrt"+Shell_File+"chain.txt", "w")
		f_rstrt.close()

		for result in IR_sampler.sample(IR_pos, iterations=1, storechain=False):
		    position = result[0]
		    f_rstrt = open("Restart/Rstrt"+Shell_File+"chain.txt", "a")
		    f_rstrt.write(" ".join(param_names))
		    f_rstrt.write("\n")
		    for k in range(position.shape[0]):
		    	#print k
		    	f_rstrt.write("%i " %k)
		    	for i in range(len(param_names)):
		    		f_rstrt.write(" %g " %position[k][i])
		    	f_rstrt.write("\n")
		        #f_rstrt.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
		    f_rstrt.close()



		### OPEN OUTPUT DATA
		with open("emcee_data/Pickles/TDE_IR_%iwalkers.pickle" %clen) as f1:
			IR_chain,IR_lnprobs = pickle.load(f1)



		#IR_acor  = IR_sampler.acor

		IR_flatchain = np.vstack(IR_chain[:,clen/2:])
		IR_flatlnprobs = np.vstack(IR_lnprobs[:,clen/2:])
				

		IR_p_opt = IR_flatchain[IR_flatlnprobs.argmax()]		














		print "ANALYSING MCMC (!)..."
		with open("emcee_data/Pickles/TDE_IR_%iwalkers.pickle" %clen) as f1:
			IR_chain,IR_lnprobs = pickle.load(f1)



		##PLOT dem WALKERS
		for k in range(IR_chain.shape[2]):
			plt.figure()
			#plt.figure()
			for i in range(IR_chain.shape[0]):
				plt.plot(IR_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(param_names[k])
				plt.xlabel('steps')
				plt.tight_layout()
			plt.savefig('emcee_data/'+Shell_File+'_TDE_%s_%iwalkers.png' %(param_names[k],clen))
			plt.clf()




		###CORNER PLOT	
		IR_flatchain = np.vstack(IR_chain[:,clen/2:])
		IR_flatlnprobs = np.vstack(IR_lnprobs[:,clen/2:])



				

		#import triangle
		import corner as triangle
		IR_fig = triangle.corner(IR_flatchain,  labels=param_names, quantiles=[0.15, 0.5, 0.85],show_titles=True, title_kwargs={"fontsize": 14},label_kwargs={"fontsize": 18})			
		IR_fig.savefig('emcee_data/src_'+Shell_File+'_TDE_Corner_Plot_%iwalkers.png' %clen)




		## Do some stats on the walkers
		from scipy.stats import scoreatpercentile as scoretpercentile



		## max posterior + percentiles
		IR_MAP_vals = IR_flatchain[IR_flatlnprobs.argmax()]
		IR_perc = scoretpercentile(IR_flatchain, [15,85], axis=0)




		filename = "emcee_Results/TDE_results_"+Shell_File+"%iwalkers.txt" %clen
		print "Printing Results"
		target = open(filename, 'w')
		target.truncate()

		IR_diff_minus = np.zeros(len(param_names))
		IR_diff_plus  = np.zeros(len(param_names))
		for i,name in enumerate(param_names):
			IR_diff_minus[i] = IR_MAP_vals[i] - IR_perc[0,i]
			IR_diff_plus[i] = IR_perc[1,i] - IR_MAP_vals[i]
			target.write("TDE_IR: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(IR_MAP_vals[i], IR_diff_plus[i], IR_diff_minus[i], name=name))
			target.write("\n")





		IR_mxprbs = zeros(nwalkers)

					

		for i in range(nwalkers):
			IR_mxprbs[i] = max(IR_lnprobs[i])


		chi2_pdf_IR = -max(IR_mxprbs)/(2.*len(t_avg) - len(param_names) - 1)
		target.write("\n")		
		target.write("Shell TDE W1 and W2 fit, reduced chi2 =  %04g" %chi2_pdf_IR)
		target.write("\n")
		target.close()


		##extra analaysis and plot best fit here:
		###FOR NOW PLOT BEST FITS:













if (Pplot):
	### PLOT POINTS
	print "PLOTTING"
	Nt=40
	tt = np.linspace(0.00, 12.,       Nt)*yr2sec




	if (Fit_ALL):
		if (Fit==False and Fit_fmin == False):
			Shell_File = '_No_Fit_'
			param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$"]
			All_diff_plus = np.zeros(11)
			All_diff_minus = np.zeros(11)
			if (Rfxd):
				Shell_File += "Rfxd_"
				#All_p_opt = [0.910846, 0.848817, 0.00366485, 0.269825, 5.38774, 0.0200015, LVbnd, t0/yr2sec, tfb/yr2sec, gam]
				#All_p_opt = [1.02466091e+00,   3.45764549e-01,  -5.71892301e-01, 7.45045151e-02,   2.71663490e+00,   1.09507827e-01, 2.46264503e-02,   1.01923360e-03,   1.35243703e+00, 3.47157797e+00]
				#All_p_opt = [1.00165688,  0.25494203,  0.94144883,  0.28295825,  8.0,  0.0875567, LVbnd, t0/yr2sec, tfb/yr2sec, gam ]
			    #All_p_opt = [1.23793989,  0.22661969,  0.75340754,  0.36330853,  7.07261987, 0.04719663,  0.04075957,  0.55869568,  0.6239391 ,  1.5733816 ] 
			   # All_p_opt = [9.23726538e-01,   9.02902839e-02,  -1.39823086e-01, 1.52954915e-01,   6.66873575e+00,   1.31067343e-01, 6.61179768e-02,   8.70825658e-05,   1.24761379e+00, 2.24192681e+00]
			    #All_p_opt = [9.23726538e-01,   9.02902839e-02,  -1.39823086e-01, 1.52954915e-01,   6.66873575e+00,   1.31067343e-01, LVbnd, t0/yr2sec, tfb/yr2sec, gam]
				#
				#All_p_opt = [  8.99271639e-01,   9.22691018e-02,  -1.42585549e-01, 1.46374543e-01,   6.79000508e+00,   1.32665069e-01, 6.66665189e-02,   9.00623387e-05,   1.26117162e+00, 2.29801380e+00,   1.02734196e-02]

				#All_p_opt = [  9.38966595e-01,   1.00046904e-01,  -1.23664774e-01, 1.43521177e-01,   8.20451155e+00,   1.81756034e-01, 3.78754285e-02,   7.60534611e-05,   8.48637449e-01, 1.90787326e+00,   1.38884649e-02]
				#All_p_opt = [  0.7,   1.00046904e-01,  1.23664774e-01, 1.43521177e-01, 1.6,   10.20451155e+00,   1.81756034e-01, 3.78754285e-02,  1.0,   0.5, 1.90787326e+00,   1.38884649e-02]
				
				#April 24
				#All_p_opt = [ 9.38084792e-01,   9.18561421e-02, -1.31365501e-01,   1.31487333e-01, 2.52415700e+00,   4.49555437e+00,   1.64801865e-01,   2.19669805e-02, 7.49921946e-05,   7.63580401e-01,   1.08668587e+00,    2.29383233e-02 ]
				#All_p_opt = [ 9.26611608e-01,   8.84820679e-02,  -1.31123118e-01, 1.39123818e-01,   2.55517364e+00,   3.98501867e+00, 1.69738012e-01,   2.22059297e-02,   7.67266183e-05, 7.99562234e-01,   1.07145682e+00,   2.33878895e-02]
				All_p_opt = [ 9.26611608e-01,   8.84820679e-02,  -1.31123118e-01, 1.39123818e-01,   2.55517364e+00,   3.98501867e+00, 1.69738012e-01,   0.04,   1.0, 0.1,   1.3,   2.33878895e-02]


			else:
				Shell_File += "Rsbl_"
				if (No_ReGro):
					Shell_File += "NoReGro_"
					##FOR no reformation
					#All_p_opt = [5.88926854e+00,   5.07748657e-01,   3.99091677e-04, 2.14734914e-01,   6.57536803e+00,   1.83177674e-01, 3.76384062e-02,   6.58923433e-01,   5.29856296e-01, 1.66304499e+00]
					## WITH DUST REFORMATION
					#All_p_opt = [ 1.54100504e+01,   3.82739483e-01,   4.66174470e-04, 1.93773202e-01,   5.75805590e+00,   1.74165532e-01, 2.51103284e-03,   2.37980357e-01,   4.64101485e-01, 3.49929669e+00 ]
					#All_p_opt = [ 1.54100504e+01,   3.82739483e-01,   4.66174470e-04, 1.93773202e-01,   5.75805590e+00,   1.74165532e-01, LVbnd, t0/yr2sec, tfb/yr2sec, gam]
					
					#All_p_opt  = [  1.54555898e+01,   3.81640899e-01,   4.66995006e-04, 1.94302693e-01,   5.78479853e+00,   1.74874406e-01, 3.88715034e-02,   4.98321895e-01,   5.07851677e-01, 1.74972141e+00,   1.00254982e-02]
					All_p_opt  = [  4.85345348e+00,   3.93356583e-01,   4.33738530e-04, 2.65052050e-01,  1.6,  6.76041754e+00,   2.21162887e-01, 2.43621429e-02,   5.58555871e-01,   5.32053290e-01, 1.34841137e+00,   1.43193852e-02]

				
				else:
					Shell_File += "ReGro_"
					#All_p_opt  = [  1.25344184e+01,   3.48208739e-01,   4.86479966e-04, 1.89160955e-01, 1.6,   6.10055565e+00,   1.99350048e-01, 4.07340721e-02,   4.85946899e-01,   5.40749964e-01, 1.61745957e+00,   1.03397755e-02]
					#All_p_opt  = [  1.53352529e+01,   3.68507850e-01,   5.10341847e-04, 2.01771863e-01,   1.72056978e+00,   7.07650311e+00, 1.52872603e-01,  4.34919585e-02,  2.89131689e-01, 5.50664403e-01,   1.89425164e+00,   9.65336541e-03]
					
					#All_p_opt  = [  1.50648220e+01,   3.50311977e-01,   5.04816357e-04, 1.98229977e-01,   1.90099855e+00,   7.68582797e+00, 1.51661512e-01,  5.00930035e-02,  2.65080426e-01, 5.00964890e-01,   1.72216616e+00,   1.03412462e-02]
					#All_p_opt  = [  1.50648220e+01,   3.50311977e-01,   5.04816357e-04, 1.98229977e-01,   1.00099855e+00,   7.68582797e+00, 1.51661512e-01,  5.00930035e-02,  1.1, 0.06,   1.1,   1.03412462e-02]
					All_p_opt  = [  1.58034282e+01,   4.30724394e-01,   4.67063482e-04, 1.79919506e-01,   8.99745907e-01,   1.02496577e+01, 2.42810823e-01,   3.06004833e-02,   7.48266626e-01, 5.98355497e-02,   1.14594283e+00,   7.06828862e-03]

		IR_p_opt = [All_p_opt[0], All_p_opt[1], All_p_opt[2], All_p_opt[3], All_p_opt[4], All_p_opt[5], All_p_opt[8], All_p_opt[9],  All_p_opt[10]]
		V_p_opt = [All_p_opt[7], All_p_opt[8], All_p_opt[9], All_p_opt[10]]

		#if (plot_solns):
		IR_perr = [All_diff_plus[0], All_diff_plus[1], All_diff_plus[2], All_diff_plus[3], All_diff_plus[4], All_diff_plus[5]]
		IR_merr = [All_diff_minus[0], All_diff_minus[1], All_diff_minus[2], All_diff_minus[3], All_diff_minus[4], All_diff_minus[5]]

		V_perr = [All_diff_plus[6], All_diff_plus[7], All_diff_plus[8], All_diff_plus[9]]
		V_merr = [All_diff_minus[6], All_diff_minus[7], All_diff_minus[8], All_diff_minus[9]]


		Nsolns = len(param_names)*2 + 1
		FI1_slns = np.zeros([Nsolns, Nt])
		FI2_slns = np.zeros([Nsolns, Nt])
		Fsrc_slns = np.zeros([Nsolns, Nt])
		All_p_slns = np.zeros([Nsolns,12])
		IR_p_slns = np.zeros([Nsolns,9])
		V_p_slns = np.zeros([Nsolns,4])
		for i in range(Nsolns):
			All_p_slns[i] = All_p_opt
			IR_p_slns[i] = IR_p_opt
			V_p_slns[i]  = V_p_opt


		k  = 0
		for j in range(10):
			for i in range(2):
				if (i==0):
					err = -All_diff_minus[j] 
				if (i==1):
					err = All_diff_plus[j] 
				All_p_slns[k][j] = All_p_opt[j] + err
				k = k+1

		for k in range(Nsolns):
			for i in range(6):
				IR_p_slns[k][i] = All_p_slns[k][i]
			for j in range(6, 10):
				V_p_slns[k][j-6]  = All_p_slns[k][j]




			## Plot all allowed solutions
		if (plot_solns):
			for j in range(Nsolns):
				for i in range(Nt):
					if (Rfxd):
						FI1_slns[j][i] = IRLC_fxdR_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn)
						FI2_slns[j][i] = IRLC_fxdR_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn)

					else:
						FI1_slns[j][i] = IRLC_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn)
						FI2_slns[j][i] = IRLC_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn)

				Fsrc_slns[j] = VLC_point(V_p_slns[j], tt, Varg, All_p_slns[j][4])





	else:
		if (Fit_fmin == False):
			if (Fit==False or Fit_IR==False or Fit_ALL==False):
				Shell_File = '_No_Fit_'
				if (MaxL):
					if (Rfxd):
						if (Src_BF):
							IR_p_opt = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476, sigML]
							perr = [0.0037, -0.0607, 0.3757, -0.1978, -0.0736, 0.0]
							merr = [0.0003, 0.1323, -0.1746, 0.2659, 0.1329, 0.0]
						else:
							#CONVERGED
							IR_p_opt = [0.7752, 0.9070, -0.4418, 0.2750, 1.1061, 0.0034]
							perr = [0.2399, 0.0156, 0.9733, 0.0182, 0.7653, 0.0]
							merr = [0.0059, 0.1848, 0.2036, 0.0094, 0.0349, 0.0]
					else:
						if (Src_BF):
							IR_p_opt = [14.4698, 0.9961, 0.2437, 0.2632, 3.5439, sigML]
							perr = [0.8689, -0.0127, -0.0034, -0.0035, 0.0035, 0.0]
							merr = [0.0093, 0.1769, 0.2438, 0.2617, 0.2791, 0.0]
						else:
							##CONVERGED
							IR_p_opt = [18.2154, 0.9518, -0.1839, 0.2839, 1.5222, 0.0014]
							perr = [0.5836, 0.0201, 0.9066, 0.0181, 0.0889, 0.0]
							merr = [0.6233, 0.1352, 0.4791, 0.0140, 0.2516, 0.0]

				else:
					if (Rfxd):
						if (Src_BF):
							IR_p_opt = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476]
							perr = [0.0037, -0.0607, 0.3757, -0.1978, -0.0736]
							merr = [0.0003, 0.1323, -0.1746, 0.2659, 0.1329]
						else:
							#Converged
							IR_p_opt = [0.7752, 0.9070, -0.4418, 0.2750, 1.1061]
							perr = [0.2399, 0.0156, 0.9733, 0.0182, 0.7653]
							merr = [0.0059, 0.1848, 0.2036, 0.0094, 0.0349]
					else:
						if (Src_BF):
							IR_p_opt = [14.4698, 0.9961, 0.2437, 0.2632, 3.5439]
							perr = [0.8689, -0.0127, -0.0034, -0.0035, 0.0035]
							merr = [0.0093, 0.1769, 0.2438, 0.2617, 0.2791]
						else:
							##CONVERGED
							IR_p_opt = [18.3361, 0.9596, -0.1268, 0.2817, 1.5030]
							perr = [0.6522, 0.0128, 0.7840, 0.0276, 0.2518, 0.0]
							merr = [0.8675, 0.1091, 0.3847, 0.0113, 0.1721, 0.0]
			
				if (Fit_Src == False):
					V_p_opt = [LVbnd, t0/yr2sec, tfb/yr2sec, gam]
		else:
			V_p_opt = [LVbnd, t0/yr2sec, tfb/yr2sec, gam]







	FsrcI1 = np.zeros(Nt)
	FVplus = np.zeros(Nt)
	FI1 = np.zeros(Nt)
	FI2 = np.zeros(Nt)
	F1chk = np.zeros(len(t_avg))


	
	LIRfit = IR_p_opt[5] ##DEFUNCT thsi is overwritten by V-p_opt value
	FsrcI1 = VLC_point(V_p_opt, tt, Varg, LIRfit)


	##Plot best fit solutions
	if (MaxL):
		for i in range(Nt):
			#FsrcI1[i] = -2.5*np.log10(Fsrc_Anl_Fit(tt[i], Dst, V_p_opt[0], V_p_opt[2], V_p_opt[1], V_p_opt[3], V_p_opt[4]))
			if (Rfxd):
				if (USE_RSfns):
					FI1[i]    = IRLC_W1_fxdR_ML_point(IR_p_opt, tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
					FI2[i]    = IRLC_W2_fxdR_ML_point(IR_p_opt, tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)
					#FVplus[i] = IRLC_fxdR_ML_point(IR_p_opt, tt[i], argVplus, RHS_table, T_table, RHS_mx, RHS_mn)
				else:
					FI1[i]    = IRLC_fxdR_ML_point(IR_p_opt, tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn)
					FI2[i]    = IRLC_fxdR_ML_point(IR_p_opt, tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn)
			else:
				if (USE_RSfns):
					FI1[i]    = IRLC_W1_ML_point(IR_p_opt, tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
					FI2[i]    = IRLC_W2_ML_point(IR_p_opt, tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)
					#FVplus[i] = IRLC_ML_point(IR_p_opt, tt[i], argVplus, RHS_table, T_table, RHS_mx, RHS_mn)
				else:
					FI1[i]    = IRLC_ML_point(IR_p_opt, tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn)
					FI2[i]    = IRLC_ML_point(IR_p_opt, tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn)

	else:
		for i in range(Nt):
			#FsrcI1[i] = -2.5*np.log10(Fsrc_Anl_Fit(tt[i], Dst, V_p_opt[0], V_p_opt[2], V_p_opt[1], V_p_opt[3], V_p_opt[4]))
			if (Rfxd):
				FI1[i]    = IRLC_W1_fxdR_point(IR_p_opt, tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
				FI2[i]    = IRLC_W2_fxdR_point(IR_p_opt, tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)
				#FVplus[i] = IRLC_fxdR_point(IR_p_opt, tt[i], argVplus, RHS_table, T_table, RHS_mx, RHS_mn)
			else:
				FI1[i]    = IRLC_W1_point(IR_p_opt, tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
				FI2[i]    = IRLC_W2_point(IR_p_opt, tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)
				#FVplus[i] = IRLC_point(IR_p_opt, tt[i], argVplus, RHS_table, T_table, RHS_mx, RHS_mn)

	# for i in range(len(t_avg)):
	# 	if (Rfxd):
	# 		F1chk[i]  = min(IRLC_fxdR_point(IR_p_opt, t_avg[i], argW1, RHS_table, T_table), 12.9)
	# 	else:
	# 		F1chk[i]  = min(IRLC_point(IR_p_opt, t_avg[i], argW1, RHS_table, T_table), 12.9)



	## Plot all allowed solutions
	if (plot_solns and Fit_ALL==False):
		IR_p_opt = np.array(IR_p_opt)
		perr = np.array(perr)
		merr = np.array(merr)


		Nsolns = 11
		FI1_slns = np.zeros([Nsolns, Nt])
		FI2_slns = np.zeros([Nsolns, Nt])
		IR_p_slns = np.zeros([Nsolns,6])
		for i in range(Nsolns):
			IR_p_slns[i] = IR_p_opt

		k  = 0
		for j in range(5):
			for i in range(2):
				if (i==0):
					err = -merr[j] 
				if (i==1):
					err = perr[j] 
				IR_p_slns[k][j] = IR_p_opt[j] + err
				k = k+1

		for j in range(Nsolns):
			for i in range(Nt):
				if (Rfxd):
					FI1_slns[j][i] = IRLC_W1_fxdR_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
					FI2_slns[j][i] = IRLC_W2_fxdR_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)

				else:
					FI1_slns[j][i] = IRLC_W1_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, T_table, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
					FI2_slns[j][i] = IRLC_W2_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, T_table, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)





	#FVtot = -2.5*np.log10(10.**(-FsrcI1/2.5) + 10.**(-FVplus/2.5))


	###PLOT###
	plt.figure()

	if (Rfxd):
		plt.title(r"$R = %g$ pc" %IR_p_opt[0])
	else:
		R0s = 0.05*IR_p_opt[0]
		plt.title(r"$R = %g L^{1/2}_{44}$ pc" %R0s)

	#plt.scatter(t_avg/day2sec, F1chk, color='orange')

	#st = plt.plot(tt/(yr2sec)*365., FVtot-3.5, linestyle = '-', color='blue', linewidth=2)
	s1 = plt.plot(tt/(yr2sec)*365., FsrcI1-3.5, linestyle = '-', color='blue', linewidth=3)


	Vdat   = plt.errorbar(tV_srt/day2sec, V_srt-3.5, yerr=V_sigsrt, linestyle="none", color='blue', alpha=0.2, elinewidth=1.5)
	Vsct   = plt.scatter(tV_srt/day2sec, V_srt-3.5, color='blue', alpha=0.2)

	Vavg     = plt.errorbar(tV_avg/day2sec, V_avg-3.5, yerr=V_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
	Vavsct   = plt.scatter(tV_avg/day2sec, V_avg-3.5, color='black', alpha=1.)


	#Vav   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)

	IR1 = plt.plot(tt/yr2sec*365., FI1, color='orange', linewidth=3)#color='#1b9e77', linewidth=3)
	IR2 = plt.plot(tt/yr2sec*365., FI2, color='red', linewidth=3)#color='#d95f02', linewidth=3)

	if (plot_solns):
		for i in range(Nsolns):
			plt.plot(tt/yr2sec*365., FI1_slns[i], color='orange', linewidth=2, alpha=0.5)#color='#1b9e77', linewidth=3)
			plt.plot(tt/yr2sec*365., FI2_slns[i], color='red', linewidth=2, alpha=0.5)#color='#d95f02', linewidth=3)

			plt.plot(tt/(yr2sec)*365., Fsrc_slns[i]-3.5, linestyle = '-', color='blue', linewidth=2, alpha=0.5)







		
	W1dat   = plt.errorbar(t_MJD/day2sec, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=0.5, elinewidth=1.5)
	W1sct   = plt.scatter(t_MJD/day2sec, W1_mag,   color='orange', alpha=0.5)

	W1av   = plt.errorbar(t_avg/day2sec, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1.0, elinewidth=1.5)
	W1as   = plt.scatter(t_avg/day2sec, W1_avg,   color='black', alpha=1.0)




	W2dat   = plt.errorbar(t_MJD/day2sec, W2_mag, yerr=W2_sig, linestyle="none", color='red', alpha=0.5, elinewidth=1.5)
	W2sct   = plt.scatter(t_MJD/day2sec, W2_mag,  color='red', alpha=0.5)

	W2av   = plt.errorbar(t_avg/day2sec, W2_avg, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
	W2as   = plt.scatter(t_avg/day2sec, W2_avg,   color='black', alpha=1.)

	plt.figtext(0.75, 0.25, "V-3.5", color='blue')
	plt.figtext(0.75, 0.45, "W1", color='orange')
	plt.figtext(0.75, 0.75, "W2", color='red')


	plt.grid(b=True, which='both')
	#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$R_d = R_0$',   r'$R_d = 0.8R_0$',   r'$R_d = 1.\bar{33}R_0$'), loc='upper right', fontsize=18)

	plt.xlabel(r"$t/t_{\rm{fb}}$")
	plt.ylabel("mag")
	#plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

	plt.ylim(plt.ylim(9.5, 14.5)[::-1])

	plt.tight_layout()

	if (Fit):
		if (Rfxd):
			Savename = "plots/BestFits"+Shell_File+"TDE_Rfxd_AnalySrc_Ntrapphi%g_NTrapnu%g_tfb%g_clen%g_%gwalkers_nne%g.png" %(Ntrap_ph, Ntrap_nu,tfb,clen,nwalkers, nne)
		else:
			Savename = "plots/BestFits"+Shell_File+"TDE_Rsubl_AnalySrc_Ntrapphi%g_NTrapnu%g_tfb%g_clen%g_%gwalkers_nne%g.png" %(Ntrap_ph, Ntrap_nu,tfb,clen,nwalkers, nne)
	else:
		if (Rfxd):
			Savename = "plots/BestFits"+Shell_File+"TDE_Rfxd_AnalySrc_Ntrapphi%g_NTrapnu%g__tfb%g_nne%g.png" %(Ntrap_ph, Ntrap_nu, tfb, nne)
		else:
			Savename = "plots/BestFits"+Shell_File+"TDE_Rsubl_AnalySrc_Ntrapphi%g_NTrapnu%g_tfb%g_nne%g.png" %(Ntrap_ph, Ntrap_nu, tfb, nne)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)


