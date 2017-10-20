import cPickle as pickle
import numpy as np
from numpy import *
import scipy as sc
from scipy import *


#from emcee.utils import MPIPool

import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm

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
TrimE = True
plot_MCMC = True
plot_ensemble = False 


perc_rng = 15


##OPTIONS THAT YOU DONT CHANGE
#Strt_Rstrt = 4
Src_BF = False
Rfxd = False
No_ReGro = True ##HAVE to change by hand right now in FLuxFn- this only changes names of files
##Trim beginning of Vband
TrimB = True
zTDE = 0.11
NTdust = 60





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


clen = 512
ndim = len(p0)
nwalkers = ndim*6
	


### OPEN OUTPUT DATA
print "ANALYSING MCMC (!)..."
#Strt_Rstrt = 3

#Shell_File_s = "Restart%i_"%Strt_Rstrt +Shell_File
Shell_File_3 = "Restart3_" +Shell_File
Shell_File_4 = "Restart4_" +Shell_File
Shell_File_5 = "Restart5_" +Shell_File
Shell_File_6 = "Restart6_" +Shell_File
Shell_File_7 = "Restart7_" +Shell_File
Shell_File_8 = "Restart8_" +Shell_File
Shell_File_9 = "Restart9_" +Shell_File
Shell_File_10 = "Restart10_" +Shell_File
Shell_File_11 = "Restart11_" +Shell_File
Shell_File_12 = "Restart12_" +Shell_File
Shell_File_13 = "Restart13_" +Shell_File
Strt_Rstrt = 4
if (TrimE):
	print "Trimming End point!"
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_4+"All_%iwalkers.pickle" %clen) as f1:
		All_chain4,All_lnprobs4 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_5+"All_%iwalkers.pickle" %clen) as f1:
		All_chain5,All_lnprobs5 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_6+"All_%iwalkers.pickle" %clen) as f1:
		All_chain6,All_lnprobs6 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_7+"All_%iwalkers.pickle" %clen) as f1:
		All_chain7,All_lnprobs7 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_8+"All_%iwalkers.pickle" %clen) as f1:
		All_chain8,All_lnprobs8 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_9+"All_%iwalkers.pickle" %clen) as f1:
		All_chain9,All_lnprobs9 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_10+"All_%iwalkers.pickle" %clen) as f1:
		All_chain10,All_lnprobs10 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_11+"All_%iwalkers.pickle" %clen) as f1:
		All_chain11,All_lnprobs11 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_12+"All_%iwalkers.pickle" %clen) as f1:
		All_chain12,All_lnprobs12 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_13+"All_%iwalkers.pickle" %clen) as f1:
		All_chain13,All_lnprobs13 = pickle.load(f1)
	End_Rstrt = 13
	lenRstrt = End_Rstrt-Strt_Rstrt

	All_chain45 = np.append(All_chain4, All_chain5, axis=1)
	All_chain456 = np.append(All_chain45, All_chain6, axis=1)
	All_chain4567 = np.append(All_chain456, All_chain7, axis=1)
	All_chain48 = np.append(All_chain4567, All_chain8, axis=1)
	All_chain49 = np.array(np.append(All_chain48, All_chain9, axis=1))
	All_chain410 = np.array(np.append(All_chain49, All_chain10, axis=1))
	All_chain = np.array(np.append(All_chain410, All_chain11, axis=1))

	All_lnprobs45 = np.append(All_lnprobs4, All_lnprobs5, axis=1)
	All_lnprobs456 = np.append(All_lnprobs45, All_lnprobs6, axis=1)
	All_lnprobs4567 = np.append(All_lnprobs456, All_lnprobs7, axis=1)
	All_lnprobs48 = np.append(All_lnprobs4567, All_lnprobs8, axis=1)
	All_lnprobs49 = np.array(np.append(All_lnprobs48, All_lnprobs9, axis=1))
	All_lnprobs410 = np.array(np.append(All_lnprobs49, All_lnprobs10, axis=1))
	All_lnprobs = np.array(np.append(All_lnprobs410, All_lnprobs11, axis=1))



else:
	print "USING ALL DATA!"
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_4+"All_%iwalkers.pickle" %clen) as f1:
		All_chain4,All_lnprobs4 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_5+"All_%iwalkers.pickle" %clen) as f1:
		All_chain5,All_lnprobs5 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_6+"All_%iwalkers.pickle" %clen) as f1:
		All_chain6,All_lnprobs6 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_7+"All_%iwalkers.pickle" %clen) as f1:
		All_chain7,All_lnprobs7 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_8+"All_%iwalkers.pickle" %clen) as f1:
		All_chain8,All_lnprobs8 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_9+"All_%iwalkers.pickle" %clen) as f1:
		All_chain9,All_lnprobs9 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_10+"All_%iwalkers.pickle" %clen) as f1:
		All_chain10,All_lnprobs10 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_11+"All_%iwalkers.pickle" %clen) as f1:
		All_chain11,All_lnprobs11 = pickle.load(f1)
	with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File_12+"All_%iwalkers.pickle" %clen) as f1:
		All_chain12,All_lnprobs12 = pickle.load(f1)
	End_Rstrt = 12
	lenRstrt = End_Rstrt-Strt_Rstrt


	All_chain45 = np.append(All_chain4, All_chain5, axis=1)
	All_chain456 = np.append(All_chain45, All_chain6, axis=1)
	All_chain4567 = np.append(All_chain456, All_chain7, axis=1)
	All_chain48 = np.append(All_chain4567, All_chain8, axis=1)
	All_chain49 = np.append(All_chain48, All_chain9, axis=1)
	All_chain = np.array(np.append(All_chain49, All_chain10, axis=1))

	All_lnprobs45 = np.append(All_lnprobs4, All_lnprobs5, axis=1)
	All_lnprobs456 = np.append(All_lnprobs45, All_lnprobs6, axis=1)
	All_lnprobs4567 = np.append(All_lnprobs456, All_lnprobs7, axis=1)
	All_lnprobs48 = np.append(All_lnprobs4567, All_lnprobs8, axis=1)
	All_lnprobs49 = np.array(np.append(All_lnprobs48, All_lnprobs9, axis=1))
	All_lnprobs = np.array(np.append(All_lnprobs49, All_lnprobs10, axis=1))



#if (TrimE):
# for Ri in range(Strt_Rstrt+1, Rstrt-1):
# 	print Ri
# 	Shell_File_it = "Restart%i_"%Ri +Shell_File 
# 	with open("plots/yeti/obscured/subl_NoRegro_TrimE/Pickles/TDE"+Shell_File_it+"All_%iwalkers.pickle" %clen) as f1:
# 		All_chains ,All_lnprob = pickle.load(f1)
# 	All_chain = np.append(All_chain, All_chains, axis=0)
# 	All_lnprobs = np.append(All_lnprob, All_lnprobs, axis=0)
#else:	
	# with open("plots/yeti/obscured/subl_NoRegro/Pickles/TDE"+Shell_File+"All_%iwalkers.pickle" %clen) as f1:
	# 	All_chain,All_lnprobs = pickle.load(f1)


clen = lenRstrt*clen

#V_acor = V_sampler.acor
# All_flatchain = np.vstack(All_chain[:,clen/4:])
All_flatchain = np.vstack(All_chain[:,clen*3/4:])
All_flatlnprobs = np.vstack(All_lnprobs[:,clen*3/4:])
		

All_p_opt = All_flatchain[All_flatlnprobs.argmax()]		




###Gelman-Rubin Statistic for convergence
#Mean for each parameter

#N parameters
#m chains per parameter
#i steps per chain

Ptot = All_chain.shape[2] ##number fo params
Ntot = All_chain.shape[1] ##number of steps per chain
Mtot = All_chain.shape[0] ##number of chains
dof = Ntot*Mtot-1.0

thav_m = np.zeros([Mtot,Ptot])
var_m = np.zeros([Mtot,Ptot])

thav = np.zeros(Ptot)
Btwn = np.zeros(Ptot)
Wthn = np.zeros(Ptot)
Vpool = np.zeros(Ptot)
PSRF = np.zeros(Ptot)


##GEllman-Rubin Statistic From: http://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/

for Pp in range(Ptot):
	for Mm in range(Mtot):
		thav_m[Mm, Pp] = np.mean(All_chain[Mm, :, Pp])
		var_m[Mm, Pp]  = (np.std(All_chain[Mm, :, Pp]))**2
	thav[Pp] = np.mean(thav_m[:,Pp])
	Btwn[Pp] = Ntot/(Mtot-1.0) * np.sum( (thav_m[:, Pp] - thav[Pp])**2 )
	Wthn[Pp] = 1.0/Mtot * np.sum(var_m[:, Pp])

	Vpool[Pp] = (Ntot-1.0)/Ntot *Wthn[Pp] + (Mtot+1.0)/(Mtot * Ntot) * Btwn[Pp]


	PSRF[Pp] = np.sqrt((dof+3.0)/(dof+1.0) * Vpool[Pp]/Wthn[Pp])

	print "PSRF for "+param_names[Pp]+ "= %g" %PSRF[Pp]


if (plot_MCMC):
	##PLOT dem WALKERS
	print  "PLOT dem WALKERS"
	colors = cm.tab20(np.linspace(0, 1, All_chain.shape[0]))
	for k in range(All_chain.shape[2]):
		plt.figure()
		#plt.figure()
		for i in range(All_chain.shape[0]):
			plt.plot(All_chain[i,:,k], drawstyle='steps', color=colors[i], marker=None, alpha=0.7, linewidth=1)
			plt.ylabel(param_names[k])
			plt.xlabel('steps')
			plt.tight_layout()
			if (TrimE):
				plt.savefig("plots/yeti/Obscured/subl_NoRegro_TrimE/"+Shell_File+'_TDE_%s_%iwalkers.png' %(param_names[k],clen))
			else:
				plt.savefig("plots/yeti/Obscured/subl_NoRegro/"+Shell_File+'_TDE_%s_%iwalkers.png' %(param_names[k],clen))
		plt.clf()





	###CORNER PLOT	
	#import triangle
	print  "PLOT dat CORNER"
	import corner as triangle
	All_fig = triangle.corner(All_flatchain,  labels=param_names, quantiles=[0.59, 0.5, 0.51],show_titles=True, title_kwargs={"fontsize": 18},label_kwargs={"fontsize": 18})			
	if (TrimE):
		plt.savefig("plots/yeti/Obscured/subl_NoRegro_TrimE/Strt%i_End%i"%(Strt_Rstrt, End_Rstrt)+Shell_File+'_TDE_Corner_Plot_%iwalkers.png' %clen )
	else:
		plt.savefig("plots/yeti/Obscured/subl_NoRegro/Strt%i_End%i_"%(Strt_Rstrt, End_Rstrt)+Shell_File+'_TDE_Corner_Plot_%iwalkers.png' %clen )




## Do some stats on the walkers
from scipy.stats import scoreatpercentile as scoretpercentile



## max posterior + percentiles
All_MAP_vals = All_flatchain[All_flatlnprobs.argmax()]
All_perc = scoretpercentile(All_flatchain, [50-perc_rng,50+perc_rng], axis=0)



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






if (plot_ensemble):
	print "PLOTTING SOUTION ENSEMBLE"
	Nt=80
	tt = np.linspace(0.00, 5000./365.,       Nt)*yr2sec

	##TABULATE T's and RHSs
	print "Creating look up tables"
	nu0 = numicron*All_p_opt[3]
	nne = All_p_opt[4]
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



	FsrcI1 = np.zeros(Nt)
	FVplus = np.zeros(Nt)
	FI1 = np.zeros(Nt)
	FI2 = np.zeros(Nt)
	F1chk = np.zeros(len(t_avg))


	
	
	#All_p_opt  =[ 4.77589089,  0.85436462,  0.86116253,  0.23955337,  7.60425767, 1.83650851,  0.03721246,  0.03448687,  0.64344651,  0.0301423 , 0.63908395,  0.01995151]



	#create array of solutions for above ranges, (12 parameters * >3 values for each parameter - ?overlap = ? solutions)
	IR_p_opt = np.array([All_p_opt[0], All_p_opt[1], All_p_opt[2], All_p_opt[3], All_p_opt[4], All_p_opt[5], All_p_opt[8], All_p_opt[9],  All_p_opt[10]])
	V_p_opt = np.array([All_p_opt[7], All_p_opt[8], All_p_opt[9], All_p_opt[10]])

	#if (plot_solns):
	IR_perr = np.array([All_diff_plus[0], All_diff_plus[1], All_diff_plus[2], All_diff_plus[3], All_diff_plus[4], All_diff_plus[5], All_diff_plus[8], All_diff_plus[9], All_diff_plus[10]])
	IR_merr = np.array([All_diff_minus[0], All_diff_minus[1], All_diff_minus[2], All_diff_minus[3], All_diff_minus[4], All_diff_minus[5], All_diff_minus[8], All_diff_minus[9], All_diff_minus[10]])

	V_perr = np.array([All_diff_plus[7], All_diff_plus[8], All_diff_plus[9], All_diff_plus[10]])
	V_merr = np.array([All_diff_minus[7], All_diff_minus[8], All_diff_minus[9], All_diff_minus[10]])




	Nsolns = len(param_names)*3#(len(param_names)-1)*2 +1
	FI1_slns = np.zeros([Nsolns, Nt])
	FI2_slns = np.zeros([Nsolns, Nt])
	Fsrc_slns = np.zeros([Nsolns, Nt])
	All_p_slns = np.zeros([Nsolns,len(All_p_opt)])
	IR_p_slns = np.zeros([Nsolns,len(IR_p_opt)])
	V_p_slns = np.zeros([Nsolns,len(V_p_opt)])
	for i in range(Nsolns):
		All_p_slns[i] = All_p_opt
		IR_p_slns[i] = IR_p_opt
		V_p_slns[i]  = V_p_opt


	for i in range(Nsolns):
		All_p_slns[i] = (All_p_opt - All_diff_minus) + (All_diff_minus + All_diff_plus)*np.random.rand(len(All_p_slns[0]))
		IR_p_slns[i] = (IR_p_opt - IR_merr) + (IR_perr + IR_merr) * np.random.rand(len(IR_p_slns[0]))
		V_p_slns[i]  = (V_p_opt - V_merr) + (V_perr + V_merr)*np.random.rand(len(V_p_slns[0]))

	# k  = 0
	# #All_p_slns[0] = All_p_opt
	# for j in range(11):
	# 	#print j
	# 	for i in range(2):
	# 		if (i==0):
	# 			err = -All_diff_minus[j] 
	# 		if (i==1):
	# 			err = All_diff_plus[j] 
	# 		All_p_slns[k][j] = All_p_opt[j] + err
	# 		k = k+1

	# for kk in range(Nsolns):
	# 	for ii in range(9):
	# 		IR_p_slns[kk][ii] = All_p_slns[kk][ii]
	# 	for jj in range(7, 11):
	# 		#print jj
	# 		V_p_slns[kk][jj-7]  = All_p_slns[kk][jj]


#Note that All_p_solns has a diff value of the kth parmater for each k, so since V-P_slns is parametrs 7,8,9,10, the first 18 rows will be identical
	#plot
	for j in range(Nsolns):
		for i in range(Nt):
			FI1_slns[j][i] = IRLC_W1_ML_point(IR_p_slns[j], tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
			FI2_slns[j][i] = IRLC_W2_ML_point(IR_p_slns[j], tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)
		Fsrc_slns[j] = VLC_point(V_p_slns[j], tt, Varg, 1.0)

	FsrcI1 = VLC_point(V_p_opt, tt, Varg, 1.0) ##DEFUNCT thsi is overwritten by V-p_opt value
	for i in range(Nt):
		FI1[i]    = IRLC_W1_ML_point(IR_p_opt, tt[i], argW1, RHS_table, Td_intrp, RHS_mx, RHS_mn, W1RSR_intrp, phis, ths, nuW1)
		FI2[i]    = IRLC_W2_ML_point(IR_p_opt, tt[i], argW2, RHS_table, Td_intrp, RHS_mx, RHS_mn, W2RSR_intrp, phis, ths, nuW2)


###PLOT###
	plt.figure()

	R0s = 0.05*IR_p_opt[0]
	plt.title(r"$R = %g L^{1/2}_{44}$ pc" %R0s)
	

	#plt.scatter(t_avg/day2sec, F1chk, color='orange')

	#st = plt.plot(tt/(yr2sec)*365., FVtot-3.5, linestyle = '-', color='blue', linewidth=2)
	s1 = plt.plot(tt/(yr2sec)*365., FsrcI1-3.5, linestyle = '-', color='blue', linewidth=3, zorder=0)


	Vdat   = plt.errorbar(tV_srt/day2sec, V_srt-3.5, yerr=V_sigsrt, linestyle="none", color='blue', alpha=0.2, elinewidth=1.5, zorder=5)
	Vsct   = plt.scatter(tV_srt/day2sec, V_srt-3.5, color='blue', alpha=0.2, zorder=5)

	Vavg     = plt.errorbar(tV_avg/day2sec, V_avg-3.5, yerr=V_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5, zorder=10)
	Vavsct   = plt.scatter(tV_avg/day2sec, V_avg-3.5, color='black', alpha=1., zorder=10)


	#Vav   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)

	IR1 = plt.plot(tt/yr2sec*365., FI1, color='orange', linewidth=3, zorder=1)#color='#1b9e77', linewidth=3)
	IR2 = plt.plot(tt/yr2sec*365., FI2, color='red', linewidth=3, zorder=1)#color='#d95f02', linewidth=3)


	for i in range(Nsolns):
		plt.plot(tt/yr2sec*365., FI1_slns[i], color='orange', linewidth=2, alpha=0.25, zorder=0)#color='#1b9e77', linewidth=3)
		plt.plot(tt/yr2sec*365., FI2_slns[i], color='red', linewidth=2, alpha=0.25, zorder=0)#color='#d95f02', linewidth=3)

		plt.plot(tt/(yr2sec)*365., Fsrc_slns[i]-3.5, linestyle = '-', color='blue', linewidth=2, alpha=0.25, zorder=0)







		
	W1dat   = plt.errorbar(t_MJD/day2sec, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=0.5, elinewidth=1.5, zorder=5)
	W1sct   = plt.scatter(t_MJD/day2sec, W1_mag,  color='orange', alpha=0.5, zorder=5)

	W1av   = plt.errorbar(t_avg/day2sec, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1.0, elinewidth=1.5, zorder=10)
	W1as   = plt.scatter(t_avg/day2sec, W1_avg,   color='black', alpha=1.0, zorder=10)




	W2dat   = plt.errorbar(t_MJD/day2sec, W2_mag, yerr=W2_sig, linestyle="none", color='red', alpha=0.5, elinewidth=1.5, zorder=5)
	W2sct   = plt.scatter(t_MJD/day2sec, W2_mag,  color='red', alpha=0.5, zorder=5)

	W2av   = plt.errorbar(t_avg/day2sec, W2_avg, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5, zorder=10)
	W2as   = plt.scatter(t_avg/day2sec, W2_avg,   color='black', alpha=1., zorder=10)

	plt.figtext(0.75, 0.25, "V-3.5", color='blue', zorder=10)
	plt.figtext(0.75, 0.45, "W1", color='orange', zorder=10)
	plt.figtext(0.75, 0.75, "W2", color='red', zorder=10)


	plt.grid(b=True, which='both', zorder=15)
	#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$R_d = R_0$',   r'$R_d = 0.8R_0$',   r'$R_d = 1.\bar{33}R_0$'), loc='upper right', fontsize=18)

	plt.xlabel(r"$t [mjd-55000]$")
	plt.ylabel("mag")
	#plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

	plt.ylim(plt.ylim(9.0, 14.5)[::-1])

	plt.tight_layout()


	if (TrimE):
		Savename = "plots/yeti/Obscured/subl_NoRegro_TrimE/"+Shell_File+"TDE_Rfxd_AnalySrc_Ntrapth%g_Ntrapphi%g_NTrapnu%g_tfb%g_clen%g_%gwalkers_nne%g.png" %(Ntrap_th, Ntrap_ph, Ntrap_nu,tfb,clen,nwalkers, nne)
	else:
		Savename = "plots/yeti/Obscured/subl_NoRegro/"+Shell_File+"TDE_Rfxd_AnalySrc_Ntrapth%g_Ntrapphi%g_NTrapnu%g_tfb%g_clen%g_%gwalkers_nne%g.png" %(Ntrap_th, Ntrap_ph, Ntrap_nu,tfb,clen,nwalkers, nne)

	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)




