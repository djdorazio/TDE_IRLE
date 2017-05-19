import cPickle as pickle
import numpy as np
from numpy import *
import scipy as sc
from scipy import *


#from emcee.utils import MPIPool

import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
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

Rstrt = 2
#RstrtFile = "Restart/Rstrt_sublR_Trap10_MaxLik__src_longerFB_chain.txt"

Pplot = True
plot_solns = False
USE_RSfns = True ##If FALSE remove W1 and W2 RSR fns from plotting

Fit_fmin = False
Fit_MC = True

Src_BF = False ## doesnt matter if Fit All (all in how set V_prior)
Rfxd = False
No_ReGro = True ##HAVE to change by hand right now in FLuxFn- this only changes names of files



##multiprocessing
NThread = 1
NTdust = 60
N_RHS = 100.
#pool = MPIPool(loadbalance=True)


##Trim beginning of Vband
TrimB = True
TrimE = True

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
Tt    =  (loadtxt("dat/Vband_Lums_All.txt",  usecols=[5], comments="|", delimiter=',')  - 55000)/(1.+zTDE) * day2sec
V_mag =  loadtxt("dat/Vband_Lums_All.txt", usecols=[1], comments="|", delimiter=',')
V_sig =  loadtxt("dat/Vband_Lums_All.txt",  usecols=[2], comments="|", delimiter=',')

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
	idelB = np.where(tV_srt <= 0.0)[0]
	tV_srt   =  np.delete(tV_srt,idelB)
	V_srt    =  np.delete(V_srt,idelB)
	V_sigsrt =  np.delete(V_sigsrt,idelB)



if (TrimE):
	idelE = np.where(tV_srt > 2000.0* day2sec)[0]
	tV_srt   =  np.delete(tV_srt,idelE)
	V_srt    =  np.delete(V_srt,idelE)
	V_sigsrt =  np.delete(V_sigsrt,idelE)



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
# W1_avsg[0] *= 0.1
# W1_avsg[1] *= 0.1
# W2_avsg[0] *= 0.1
# W2_avsg[1] *= 0.1

# W1_avsg[8] *= 5.0
# W1_avsg[7] *= 5.0
# W1_avsg[6] *= 2.0

# W2_avsg[8] *= 8.0
# W2_avsg[7] *= 5.0
# W2_avsg[6] *= 3.0



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
nne = 1.6  ##emission/absortption efficeincy power law  --  Draine or 1.6 
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




##TABULATE T's and RHSs
print "Creating look up tables"
NT = NTdust
RHS_table = np.zeros(NT)
TT_table = np.linspace(1., 1800., NT)
nus_RHS = np.linspace(1.0, np.log(5.*numicron), N_RHS)
for i in range(NT):
#	RHS_table[i] = T_RHS(TT_table[i], nu0, nne)
	RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
RHS_mx = RHS_table[len(RHS_table)-1]
RHS_mn = RHS_table[0]


###MAKE TDUST INTERP
print "Making Tdust interp"
#temp rename

Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table,)

plt.figure()
plt.scatter(RHS_table, Td_intrp(RHS_table))
plt.savefig('T_Interp.png')

plt.figure()
plt.scatter(TT_table, RHS_table)
plt.savefig('RHSTABLE.png')




if (Fit_fmin):
	from scipy import optimize
	from scipy.optimize import fmin


	if (Rfxd):
			if (TrimE):
				Shell_File = "_TrimE_FminML_FitALLIRandOPt_FxdR_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
			else:
				Shell_File = "_FminML_FitALLIRandOPt_FxdR_Trap%g_" %Ntrap_nu
			param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", "k", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$", r"$\sigma^{\rm{Vbnd}}_{\rm{ML}}$"]
		
			
			#p0 = [ 9.38084792e-01,   9.18561421e-02, -1.31365501e-01,   1.31487333e-01, 2.52415700e+00,   4.49555437e+00,   1.64801865e-01,   2.19669805e-02, 7.49921946e-05,   7.63580401e-01,   1.08668587e+00,    2.29383233e-02 ]
			p0 = [ 0.839886,  0.457649,  0.491424,  2.82186e-07,  5.88129, 2.54017,  0.0535668,  0.0100001,  1.19606,  0.127709, 0.772826,  0.00662339]
			

			##TABULATE T's and RHSs
			print "Creating look up tables"
			nu0 = numicron*p0[3]
			nne = p0[4]
			NT = NTdust
			RHS_table = np.zeros(NT)
			TT_table = np.linspace(1., 1800., NT)
			nus_RHS = np.linspace(0.0, 3.0*numicron, N_RHS)
			for i in range(NT):
			#	RHS_table[i] = T_RHS(TT_table[i], nu0, nne)
				RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
			RHS_mx = RHS_table[len(RHS_table)-1]
			RHS_mn = RHS_table[0]


			###MAKE TDUST INTERP
			print "Making Tdust interp"
			#temp rename

			Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table,)

			plt.figure()
			plt.scatter(RHS_table, Td_intrp(RHS_table))
			plt.savefig('T_Interp.png')

			plt.figure()
			plt.scatter(TT_table, RHS_table)
			plt.savefig('RHSTABLE.png')

			IR_p0 = np.array(p0)
			popt  = sc.optimize.fmin(IRTDE_fxdR_FitAll_Err2_fmin,  IR_p0, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg), full_output=1, disp=False, ftol=0.1)[0]
	else:
			if (TrimE):
				if (No_ReGro):
					Shell_File = "_TrimeE_FminML_FitALLIRandOPt_sblR_NOREGROW_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
				else:
					Shell_File = "_TrimeE_FminML_FitALLIRandOPt_sblR_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
			else:
				if (No_ReGro):
					Shell_File = "_FminML_FitALLIRandOPt_sblR_NOREGROW_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
				else:
					Shell_File = "_FminML_FitALLIRandOPt_sblR_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)

			param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", "k", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$", r"$\sigma^{\rm{Vbnd}}_{\rm{ML}}$"]
			
		
			p0  = [  3.63219272e+00,   8.10934666e-01,   8.65639639e-01, 2.43747544e-01,   5.70042738e+00,   5.07455278e+00, 3.45666610e-02,   4.11066719e-02,   1.00817126e+00, 6.65924733e-02,   1.00789213e+00,   4.28704440e-03]

			##TABULATE T's and RHSs
			print "Creating look up tables"
			nu0 = numicron*p0[3]
			nne = p0[4]
			NT = NTdust
			RHS_table = np.zeros(NT)
			TT_table = np.linspace(1., 1800., NT)
			for i in range(NT):
			#	RHS_table[i] = T_RHS(TT_table[i], nu0, nne)
				RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
			RHS_mx = RHS_table[len(RHS_table)-1]
			RHS_mn = RHS_table[0]


			###MAKE TDUST INTERP
			print "Making Tdust interp"
			#temp rename

			Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table,)

			plt.figure()
			plt.scatter(RHS_table, Td_intrp(RHS_table))
			plt.savefig('T_Interp.png')

			plt.figure()
			plt.scatter(TT_table, RHS_table)
			plt.savefig('RHSTABLE.png')


			IR_p0 = np.array(p0)
			popt  = sc.optimize.fmin(IRTDE_sblR_FitAll_Err2_fmin,  IR_p0, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg), full_output=1, disp=False, ftol=0.1)[0]

				



	All_p_opt = popt


	filename = "emcee_Results/TDE_results_"+Shell_File+".txt"
	print "Printing Results"
	target = open(filename, 'w')
	target.truncate()

	for i in range(len(popt)):
		target.write(param_names[i]+" = %g" %popt[i])
		target.write("\n")

	target.close()















if (Fit_MC):
	import emcee

	if (Rfxd):
			if (TrimE):
				Shell_File = "_TrimE_FitIRandOpt_FxdR_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
			else:
				Shell_File = "_FitIRandOpt_FxdR_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
			
			param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", "k", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$", r"$\sigma^{\rm{Vbnd}}_{\rm{ML}}$"]
			p0 = [ 9.38084792e-01,   9.18561421e-02, -1.31365501e-01,   1.31487333e-01, 2.52415700e+00,   4.49555437e+00,   1.64801865e-01,   2.19669805e-02, 7.49921946e-05,   7.63580401e-01,   1.08668587e+00,    2.29383233e-02 ]


						##TABULATE T's and RHSs
			print "Creating look up tables"
			nu0 = numicron*p0[3]
			nne = p0[4]
			NT = NTdust
			RHS_table = np.zeros(NT)
			TT_table = np.linspace(1., 1800., NT)
			for i in range(NT):
			#	RHS_table[i] = T_RHS(TT_table[i], nu0, nne)
				RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
			RHS_mx = RHS_table[len(RHS_table)-1]
			RHS_mn = RHS_table[0]


			###MAKE TDUST INTERP
			print "Making Tdust interp"
			#temp rename

			Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table,)

			plt.figure()
			plt.scatter(RHS_table, Td_intrp(RHS_table))
			plt.savefig('T_Interp.png')

			plt.figure()
			plt.scatter(TT_table, RHS_table)
			plt.savefig('RHSTABLE.png')

			
			ndim = len(p0)
			nwalkers = ndim*4
																											
			All_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_fxdR_ALL_posterior, threads=NThread, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg))


	else:
		if (TrimE):
			if (No_ReGro):
				Shell_File = "_TrimE_sublR_NoRegro_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
			else:
				Shell_File = "_TrimE_sublR_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
		else:
			if (No_ReGro):
				Shell_File = "_sublR_NoRegro_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)
			else:
				Shell_File = "_sublR_NTrapnu%g_NTrapph%g_NTrapth%g_" %(Ntrap_nu, Ntrap_ph, Ntrap_th)




		param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", "k", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$", r"$\sigma^{\rm{Vbnd}}_{\rm{ML}}$"]
		p0  = [ 3.63219272e+00,   8.10934666e-01,   8.65639639e-01, 2.43747544e-01,   5.70042738e+00,   5.07455278e+00, 3.45666610e-02,   4.11066719e-02,   1.00817126e+00, 6.65924733e-02,   1.00789213e+00,   4.28704440e-03]


		##TABULATE T's and RHSs
		print "Creating look up tables"
		nu0 = numicron*p0[3]
		nne = p0[4]
		NT = NTdust
		RHS_table = np.zeros(NT)
		TT_table = np.linspace(1., 1800., NT)
		for i in range(NT):
		#	RHS_table[i] = T_RHS(TT_table[i], nu0, nne)
			RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
		RHS_mx = RHS_table[len(RHS_table)-1]
		RHS_mn = RHS_table[0]


		###MAKE TDUST INTERP
		print "Making Tdust interp"
		#temp rename

		Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table,)

		plt.figure()
		plt.scatter(RHS_table, Td_intrp(RHS_table))
		plt.savefig('T_Interp.png')

		plt.figure()
		plt.scatter(TT_table, RHS_table)
		plt.savefig('RHSTABLE.png')


		ndim = len(p0)
		nwalkers = ndim*6

		All_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_ALL_posterior, threads=NThread, args=(t_avg, tV_avg, argW1, argW2, Varg, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, W2RSR_intrp, phis, ths, nuW1, nuW2, W1_avg, W1_avsg, W2_avg, W2_avsg, V_avg, V_avsg))


	

	if (Rstrt>0):		
		RstrtFN = Rstrt - 1
		RstrtFile = "Restart/Restart%i_"%RstrtFN+Shell_File+"chain.txt"
		##define walker init from file
		#RstrtFile = "Restart/Rstrt"+Shell_File+"chain.txt"
		print "Restarting from File "+RstrtFile 
		in0 = np.zeros([ndim, nwalkers])
		for i in range(0, len(param_names)):
			in0[i] = np.array(np.genfromtxt(RstrtFile, usecols=i+1, comments='$'))

		All_walker_p0  = np.transpose(in0)
	else:
		All_p0 = np.array(p0)
		All_walker_p0 = np.random.normal(All_p0, np.abs(All_p0)*1E-4, size=(nwalkers, ndim))


	Shell_File = "Restart%i_"%Rstrt +Shell_File 
	
	clen = 256
	#for ():
	#run as iterable
	#acor function
	All_pos,_,_ = All_sampler.run_mcmc(All_walker_p0, clen)
	#manipulte (replace) fidn outliers 2std form median in beginning
	print "emcee DONE RUNNING"



	print "SAVING THE PICKLE mmmmm"
	with open("emcee_data/Pickles/TDE_All_%iwalkers.pickle" %clen, "w") as f1:
		pickle.dump((All_sampler.chain, All_sampler.lnprobability), f1)


			


	### OPEN OUTPUT DATA
	print "ANALYSING MCMC (!)..."
	with open("emcee_data/Pickles/TDE_All_%iwalkers.pickle" %clen) as f1:
		All_chain,All_lnprobs = pickle.load(f1)



	#V_acor = V_sampler.acor

	All_flatchain = np.vstack(All_chain[:,clen/4:])
	All_flatlnprobs = np.vstack(All_lnprobs[:,clen/4:])
			

	All_p_opt = All_flatchain[All_flatlnprobs.argmax()]		



	##record final state for restart
	print  "WRITING RESTART FILE"
	f_rstrt = open("Restart/"+Shell_File+"chain.txt", "w")
	f_rstrt.close()

	for result in All_sampler.sample(All_pos, iterations=1, storechain=False):
	    position = result[0]
	    f_rstrt = open("Restart/"+Shell_File+"chain.txt", "a")
	    for k in range(position.shape[0]):
	    	#print k
	    	f_rstrt.write("%i  %g %g %g %g %g %g %g %g %g %g %g %g" %(k, position[k][0], position[k][1], position[k][2], position[k][3], position[k][4], position[k][5], position[k][6], position[k][7], position[k][8], position[k][9], position[k][10], position[k][11]))
	    	f_rstrt.write("\n")
	        #f_rstrt.write("{0:4d} {1:s}\n".format(k, " ".join(position[k])))
	    f_rstrt.close()





	##PLOT dem WALKERS
	print  "PLOT dem WALKERS"
	for k in range(All_chain.shape[2]):
		plt.figure()
		#plt.figure()
		for i in range(All_chain.shape[0]):
			plt.plot(All_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
			plt.ylabel(param_names[k])
			plt.xlabel('steps')
			plt.tight_layout()
			
			plt.savefig('emcee_data/'+Shell_File+'_TDE_%s_%iwalkers.png' %(param_names[k],clen))

		plt.clf()




	###CORNER PLOT	
	All_flatchain = np.vstack(All_chain[:,clen/4:])
	All_flatlnprobs = np.vstack(All_lnprobs[:,clen/4:])



			

	#import triangle
	print  "PLOT dat CORNER"
	import corner as triangle
	All_fig = triangle.corner(All_flatchain,  labels=param_names, quantiles=[0.15, 0.5, 0.85],show_titles=True, title_kwargs={"fontsize": 14},label_kwargs={"fontsize": 18})			
	All_fig.savefig('emcee_data/'+Shell_File+'_TDE_Corner_Plot_%iwalkers.png' %clen)




	## Do some stats on the walkers
	from scipy.stats import scoreatpercentile as scoretpercentile



	## max posterior + percentiles
	All_MAP_vals = All_flatchain[All_flatlnprobs.argmax()]
	All_perc = scoretpercentile(All_flatchain, [15,85], axis=0)



	print "PRINTING RESULTS TO FILE"
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





























if (Pplot):
	### PLOT POINTS
	print "PLOTTING"
	Nt=120
	tt = np.linspace(0.00, 50.,       Nt)*yr2sec



	if (Fit_MC==False and Fit_fmin == False):
		Shell_File = '_No_Fit_'
		param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} \mu\rm{m}\nu_0$", r"$L_{45}$", r"$\sigma_{\rm{ML}}$", r"$L^{V}_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$"]
		All_diff_plus = np.zeros(11)
		All_diff_minus = np.zeros(11)
		if (Rfxd):
			Shell_File += "Rfxd_"

			#All_p_opt = [9.35220121e-01,   3.46767905e-01,   7.96086756e-01,   9.97171726e-08, 5.18182796e-01,   1.41984700e+00,   1.92290030e-01,   1.34571892e-02, 3.48190849e-01,   2.00052884e-01,   5.93820751e-01,   5.83430920e-05 ]
			All_p_opt = [7.28841671e-01,   6.62484534e-01,   8.19048508e-01,   1.05024448e-07, 8.54690827e+00,   3.01452951e+00,   1.63695778e-01,   9.50280404e-03, 1.41175127e+00,  4.08711266e-04,   3.32723164e-01,   3.10365790e-03 ]
		else:
			Shell_File += "Rsbl_"
			if (No_ReGro):
				Shell_File += "NoReGro_"
				
				All_p_opt  = [  3.63219272e+00,   8.10934666e-01,   8.65639639e-01, 2.43747544e-01,   5.70042738e+00,   5.07455278e+00, 3.45666610e-02,   4.11066719e-02,   1.00817126e+00, 6.65924733e-02,   1.00789213e+00,   4.28704440e-03]
			
			else:
				Shell_File += "ReGro_"
			
				All_p_opt  = [1.53069822e+00,   6.66538508e-01,   9.23748358e-01,   6.36056796e-02, 1.23992169e-01,   6.09638524e+00,   2.53046481e-01,   6.41176692e-02, 7.97110106e-01,   8.10336449e-02,   1.07679385e+00,   3.65140372e-03 ]



	##TABULATE T's and RHSs
	print "Creating look up tables"
	nu0 = numicron*All_p_opt[3]
	nne = All_p_opt[4]
	NT = NTdust
	RHS_table = np.zeros(NT)
	TT_table = np.linspace(1., 1800., NT)
	for i in range(NT):
	#	RHS_table[i] = T_RHS(TT_table[i], nu0, nne)
		RHS_table[i] = T_RHS(nus_RHS, TT_table[i], nu0, nne)
	RHS_mx = RHS_table[len(RHS_table)-1]
	RHS_mn = RHS_table[0]


	###MAKE TDUST INTERP
	print "Making Tdust interp"
	#temp rename

	Td_intrp = sc.interpolate.interp1d(RHS_table, TT_table,)

	plt.figure()
	plt.scatter(RHS_table, Td_intrp(RHS_table))
	plt.savefig('T_Interp.png')

	plt.figure()
	plt.scatter(TT_table, RHS_table)
	plt.savefig('RHSTABLE.png')

	IR_p_opt = [All_p_opt[0], All_p_opt[1], All_p_opt[2], All_p_opt[3], All_p_opt[4], All_p_opt[5], All_p_opt[8], All_p_opt[9],  All_p_opt[10]]
	V_p_opt = [All_p_opt[7], All_p_opt[8], All_p_opt[9], All_p_opt[10]]

	# #if (plot_solns):
	# IR_perr = [All_diff_plus[0], All_diff_plus[1], All_diff_plus[2], All_diff_plus[3], All_diff_plus[4], All_diff_plus[5], All_diff_plus[8], All_diff_plus[9], All_diff_plus[10]]
	# IR_merr = [All_diff_minus[0], All_diff_minus[1], All_diff_minus[2], All_diff_minus[3], All_diff_minus[4], All_diff_minus[5], All_diff_minus[8], All_diff_minus[9], All_diff_minus[10]]

	# V_perr = [All_diff_plus[7], All_diff_plus[8], All_diff_plus[9], All_diff_plus[10]]
	# V_merr = [All_diff_minus[7], All_diff_minus[8], All_diff_minus[9], All_diff_minus[10]]


	# Nsolns = len(param_names)*2 + 1
	# FI1_slns = np.zeros([Nsolns, Nt])
	# FI2_slns = np.zeros([Nsolns, Nt])
	# Fsrc_slns = np.zeros([Nsolns, Nt])
	# All_p_slns = np.zeros([Nsolns,len(All_p_opt)-1])
	# IR_p_slns = np.zeros([Nsolns,len(IR_p_opt)-1])
	# V_p_slns = np.zeros([Nsolns,len(V_p_opt)-1])
	# for i in range(Nsolns):
	# 	All_p_slns[i] = All_p_opt
	# 	IR_p_slns[i] = IR_p_opt
	# 	V_p_slns[i]  = V_p_opt


	# k  = 0
	# for j in range(12):
	# 	for i in range(2):
	# 		if (i==0):
	# 			err = -All_diff_minus[j] 
	# 		if (i==1):
	# 			err = All_diff_plus[j] 
	# 		All_p_slns[k][j] = All_p_opt[j] + err
	# 		k = k+1

	# for k in range(Nsolns):
	# 	for i in range(9):
	# 		IR_p_slns[k][i] = All_p_slns[k][i]
	# 	for j in range(9, 12):
	# 		V_p_slns[k][j-6]  = All_p_slns[k][j]




	# 	## Plot all allowed solutions
	# if (plot_solns):
	# 	for j in range(Nsolns):
	# 		for i in range(Nt):
	# 			if (Rfxd):
	# 				FI1_slns[j][i] = IRLC_W1_fxdR_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
	# 				FI2_slns[j][i] = IRLC_W2_fxdR_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)

	# 			else:
	# 				FI1_slns[j][i] = IRLC_W1_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
	# 				FI2_slns[j][i] = IRLC_W2_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)

	# 		Fsrc_slns[j] = VLC_point(V_p_slns[j], tt, Varg, 1.0)





	# else:
	# 	if (Fit_fmin == False):
	# 		if (Fit_MC==False or Fit_IR==False or Fit_ALL==False):
	# 			Shell_File = '_No_Fit_'
	# 			if (MaxL):
	# 				if (Rfxd):
	# 					if (Src_BF):
	# 						IR_p_opt = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476, sigML]
	# 						perr = [0.0037, -0.0607, 0.3757, -0.1978, -0.0736, 0.0]
	# 						merr = [0.0003, 0.1323, -0.1746, 0.2659, 0.1329, 0.0]
	# 					else:
	# 						#CONVERGED
	# 						IR_p_opt = [0.7752, 0.9070, -0.4418, 0.2750, 1.1061, 0.0034]
	# 						perr = [0.2399, 0.0156, 0.9733, 0.0182, 0.7653, 0.0]
	# 						merr = [0.0059, 0.1848, 0.2036, 0.0094, 0.0349, 0.0]
	# 				else:
	# 					if (Src_BF):
	# 						IR_p_opt = [14.4698, 0.9961, 0.2437, 0.2632, 3.5439, sigML]
	# 						perr = [0.8689, -0.0127, -0.0034, -0.0035, 0.0035, 0.0]
	# 						merr = [0.0093, 0.1769, 0.2438, 0.2617, 0.2791, 0.0]
	# 					else:
	# 						##CONVERGED
	# 						IR_p_opt = [18.2154, 0.9518, -0.1839, 0.2839, 1.5222, 0.0014]
	# 						perr = [0.5836, 0.0201, 0.9066, 0.0181, 0.0889, 0.0]
	# 						merr = [0.6233, 0.1352, 0.4791, 0.0140, 0.2516, 0.0]

	# 			else:
	# 				if (Rfxd):
	# 					if (Src_BF):
	# 						IR_p_opt = [0.9696, 0.8149, 0.0232, 0.2672, 1.7476]
	# 						perr = [0.0037, -0.0607, 0.3757, -0.1978, -0.0736]
	# 						merr = [0.0003, 0.1323, -0.1746, 0.2659, 0.1329]
	# 					else:
	# 						#Converged
	# 						IR_p_opt = [0.7752, 0.9070, -0.4418, 0.2750, 1.1061]
	# 						perr = [0.2399, 0.0156, 0.9733, 0.0182, 0.7653]
	# 						merr = [0.0059, 0.1848, 0.2036, 0.0094, 0.0349]
	# 				else:
	# 					if (Src_BF):
	# 						IR_p_opt = [14.4698, 0.9961, 0.2437, 0.2632, 3.5439]
	# 						perr = [0.8689, -0.0127, -0.0034, -0.0035, 0.0035]
	# 						merr = [0.0093, 0.1769, 0.2438, 0.2617, 0.2791]
	# 					else:
	# 						##CONVERGED
	# 						IR_p_opt = [18.3361, 0.9596, -0.1268, 0.2817, 1.5030]
	# 						perr = [0.6522, 0.0128, 0.7840, 0.0276, 0.2518, 0.0]
	# 						merr = [0.8675, 0.1091, 0.3847, 0.0113, 0.1721, 0.0]
			
	# 			if (Fit_Src == False):
	# 				V_p_opt = [LVbnd, t0/yr2sec, tfb/yr2sec, gam]
	# 	else:
	# 		V_p_opt = [LVbnd, t0/yr2sec, tfb/yr2sec, gam]







	FsrcI1 = np.zeros(Nt)
	FVplus = np.zeros(Nt)
	FI1 = np.zeros(Nt)
	FI2 = np.zeros(Nt)
	F1chk = np.zeros(len(t_avg))


	
	FsrcI1 = VLC_point(V_p_opt, tt, Varg, 1.0) ##DEFUNCT thsi is overwritten by V-p_opt value


	##Plot best fit solutions
	for i in range(Nt):
		#FsrcI1[i] = -2.5*np.log10(Fsrc_Anl_Fit(tt[i], Dst, V_p_opt[0], V_p_opt[2], V_p_opt[1], V_p_opt[3], V_p_opt[4]))
		if (Rfxd):
			if (USE_RSfns):
				FI1[i]    = IRLC_W1_fxdR_ML_point(IR_p_opt, tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
				FI2[i]    = IRLC_W2_fxdR_ML_point(IR_p_opt, tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)
				#FVplus[i] = IRLC_fxdR_ML_point(IR_p_opt, tt[i], argVplus, RHS_table, Td_intrp, RHS_mx, RHS_mn)
			else:
				FI1[i]    = IRLC_fxdR_ML_point(IR_p_opt, tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn)
				FI2[i]    = IRLC_fxdR_ML_point(IR_p_opt, tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn)
		else:
			if (USE_RSfns):
				FI1[i]    = IRLC_W1_ML_point(IR_p_opt, tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
				FI2[i]    = IRLC_W2_ML_point(IR_p_opt, tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)
				#FVplus[i] = IRLC_ML_point(IR_p_opt, tt[i], argVplus, RHS_table, Td_intrp, RHS_mx, RHS_mn)
			else:
				FI1[i]    = IRLC_ML_point(IR_p_opt, tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn)
				FI2[i]    = IRLC_ML_point(IR_p_opt, tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn)




	# ## Plot all allowed solutions
	# if (plot_solns and Fit_ALL==False):
	# 	IR_p_opt = np.array(IR_p_opt)
	# 	perr = np.array(perr)
	# 	merr = np.array(merr)


	# 	Nsolns = 11
	# 	FI1_slns = np.zeros([Nsolns, Nt])
	# 	FI2_slns = np.zeros([Nsolns, Nt])
	# 	IR_p_slns = np.zeros([Nsolns,6])
	# 	for i in range(Nsolns):
	# 		IR_p_slns[i] = IR_p_opt

	# 	k  = 0
	# 	for j in range(5):
	# 		for i in range(2):
	# 			if (i==0):
	# 				err = -merr[j] 
	# 			if (i==1):
	# 				err = perr[j] 
	# 			IR_p_slns[k][j] = IR_p_opt[j] + err
	# 			k = k+1

	# 	for j in range(Nsolns):
	# 		for i in range(Nt):
	# 			if (Rfxd):
	# 				FI1_slns[j][i] = IRLC_W1_fxdR_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
	# 				FI2_slns[j][i] = IRLC_W2_fxdR_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)

	# 			else:
	# 				FI1_slns[j][i] = IRLC_W1_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
	# 				FI2_slns[j][i] = IRLC_W2_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)





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

	plt.xlabel(r"$t [day]$")
	plt.ylabel("mag")
	#plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

	plt.ylim(plt.ylim(9.5, 14.5)[::-1])

	plt.tight_layout()


	if (Fit_MC):
		if (Rfxd):
			Savename = "plots/BestFits"+Shell_File+"TDE_Rfxd_AnalySrc_Ntrapth%g_Ntrapphi%g_NTrapnu%g_tfb%g_clen%g_%gwalkers_nne%g.png" %(Ntrap_th, Ntrap_ph, Ntrap_nu,tfb,clen,nwalkers, nne)
		else:
			Savename = "plots/BestFits"+Shell_File+"TDE_Rsubl_AnalySrc_Ntrapth%g_Ntrapphi%g_NTrapnu%g_tfb%g_clen%g_%gwalkers_nne%g.png" %(Ntrap_th, Ntrap_ph, Ntrap_nu,tfb,clen,nwalkers, nne)
	else:
		if (Rfxd):
			Savename = "plots/BestFits"+Shell_File+"TDE_Rfxd_AnalySrc_Ntrapth%g_Ntrapphi%g_NTrapnu%g__tfb%g_nne%g.png" %(Ntrap_th,Ntrap_ph, Ntrap_nu, tfb, nne)
		else:
			Savename = "plots/BestFits"+Shell_File+"TDE_Rsubl_AnalySrc_Ntrapth%g_Ntrapphi%g_NTrapnu%g_tfb%g_nne%g.png" %(Ntrap_th,Ntrap_ph, Ntrap_nu, tfb, nne)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)


