import cPickle as pickle
import numpy as np
from numpy import *
import scipy as sc
from scipy import *

import emcee

import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

from scipy import optimize
from scipy.optimize import fmin

from FluxFuncs_TDEs import *
from emcee_Funcs_TDEs import *

################################
###############################
### OPTIONS
################################
################################
Pplot = True
Fit = True
Fit_Src = False
Fit_IR = True
QvFit = False
Rfxd = True

##multiprocessing
NThread = 32


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
Tt   =  (loadtxt("dat/Vband_Lums_All.txt", usecols=[5], comments="|", delimiter=',')  - 55000)/(1.+zTDE) 
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
	if (abs((tV_srt[i] - tV_srt[i+1])) > 100 ):
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
t_MJD = (np.genfromtxt("dat/IR_Lums_All.txt",usecols=1, comments="|"))/(1.+zTDE) ## put in rest frame

W1_mag = np.genfromtxt("dat/IR_Lums_All.txt",usecols=2, comments="|")
W2_mag = np.genfromtxt("dat/IR_Lums_All.txt",usecols=3, comments="|")

W1_sig = np.genfromtxt("dat/IR_Lums_All.txt",usecols=4, comments="|")
W2_sig = np.genfromtxt("dat/IR_Lums_All.txt",usecols=5, comments="|")

##binned
t_avg = (np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=0, comments="|"))/(1.+zTDE) ## put in rest frame

W1_avg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=1, comments="|")
W2_avg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=2, comments="|")

W1_avsg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=3, comments="|")
W2_avsg = np.genfromtxt("dat/IR_Lums_Avg.txt",usecols=4, comments="|")


t_avg = np.array(t_avg)
W1_avg = np.array(W1_avg)
W2_avg = np.array(W2_avg)
W1_avsg = np.array(W1_avsg)
W2_avsg = np.array(W2_avsg)


## subtract off galaxy background form first two preflare epochs
# magW1net = -2.5*np.log10(10.**(-W1_avg/2.5) - 10.**(-12.9/2.5))
# magW2net = -2.5*np.log10(10.**(-W2_avg/2.5) - 10.**(-11.26/2.5))
# W1_avg = magW1net
# W2_avg = magW2net

##########################
## SOURCE PROPERTIES
##########################
Lav = 10.**45 #erg/s
Dst = 0.513*10**9*pc2cm
tfb = 1.0 * yr2sec
t0 =  175./365. * yr2sec

## TEST VALUES
### DUST stuff
nne = 1.6       ##emission/absortption efficeincy power law  --  Draine or 1.6 
nu0 = numicron/3.0  ## emission/absortption efficeincy cutoff freq
Rde = 0.001*pc2cm #* np.sqrt(Lav/10.**46) * pc2cm     ## inner edge of dust
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

###
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3

Vmx = c/(5.45*10**(-5) - 0.44*10**(-5))
Vmn = c/(5.45*10**(-5) + 0.44*10**(-5))


### zero point fluxes - check these...
nuVbnd = c/(5.45*10**(-5))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-21)*8.8560*10**(13)#(W1mn + W1mx)/2
FW2Rel = 1.71787*10**(-21)*6.4451*10**(13)#(W2mn + W2mx)/2


##TABULATE T's and RHSs
print "Creating look up tables"
NT = 1800
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 1800., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)


## measured background flux to add to TDE LC
FV_gal = 10.**(-17.6/2.5) * FVbndRel
FW1_gal = 10.**(-12.9/2.5) * FW1Rel
FW2_gal = 10.**(-11.26/2.5) * FW2Rel


###CALCUALTE HOW MUCH VBAND contributes to IR - does FV_gal go atend here?
#argVplus = [FVbndRel, Vmn, Vmx, Dst, Lav, tfb, n0, pp, aeff, nne, t0, Rde, 0.0]
argVplus = [FVbndRel, Vmn, Vmx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, 0.0]

argtst = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, 100., t0, etaR]

argW1 = [FW1Rel, W1mn, W1mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW1_gal]
argW2 = [FW2Rel, W2mn, W2mx, Dst, tfb, n0, pp, aeff, nne, t0, Rde, FW2_gal]
Varg  = [Dst, FVbndRel, FV_gal]


if (Fit):
	if (Fit_Src):
		##########################
		## FIT SOURCE PROPERTIES
		##########################

		Shell_File = "_Fit_Vband_"
		param_names = [r"$L_0$",r"$t_0$",r"$t_{fb}$", r"$\gamma$"]
		p0V = [0.7, t0/yr2sec, tfb/yr2sec, 2.3]
		ndim = len(p0V)
		nwalkers = ndim*4
		#V_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_V_posterior, threads=NThread,args=(tV_srt, Varg, V_srt, V_sigsrt))
		V_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_V_posterior, threads=NThread,args=(tV_avg, Varg, V_avg, V_avsg))




		V_p0 = np.array(p0V)

		V_walker_p0 = np.random.normal(V_p0, np.abs(V_p0)*1E-3, size=(nwalkers, ndim))

					
		clen = 1024#4048
		V_pos,_,_ = V_sampler.run_mcmc(V_walker_p0, clen)



		print "SAVING THE PICKLE mmmmm"
		with open("emcee_data/Pickles/TDE_V_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((V_sampler.chain, V_sampler.lnprobability), f1)


				


		### OPEN OUTPUT DATA
		with open("emcee_data/Pickles/TDE_V_%iwalkers.pickle" %clen) as f1:
			V_chain,V_lnprobs = pickle.load(f1)





		V_flatchain = np.vstack(V_chain[:,clen/4:])
		V_flatlnprobs = np.vstack(V_lnprobs[:,clen/4:])
				

		V_p_opt = V_flatchain[V_flatlnprobs.argmax()]		







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
		V_flatchain = np.vstack(V_chain[:,clen/2:])
		V_flatlnprobs = np.vstack(V_lnprobs[:,clen/2:])



				

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


		chi2_pdf_V = -max(V_mxprbs)/(len(tV_srt) - len(param_names) - 1)
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
		if (Rfxd):
			Shell_File = "_FitALLIR_FxdR_Trap80_"
			param_names = [r"$R_d$ [pc]",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$c^{-1} mu\rm{m}\nu_0$", r"$L_{45}$"]
			p0IR = [0.8, np.cos(thetTst), np.sin(JJt), nu0/numicron, Lav/10.**45]
			ndim = len(p0IR)
			nwalkers = ndim*6

			IR_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_fxdR_posterior, threads=NThread,args=(t_avg, argW1, argW2, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg))


		else:
			Shell_File = "_FitALLIR_sublR_Trap80_"
			param_names = [r"$\eta_R$",r"$\cos{\theta_T}$",r"$\sin(J)$", r"$\mu\rm{m}\nu_0c^{-1}$", r"$L_{45}$"]
			p0IR = [etaR, np.cos(thetTst), np.sin(JJt), nu0/numicron, Lav/10.**45]
			ndim = len(p0IR)
			nwalkers = ndim*6

			IR_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_IR_posterior, threads=NThread,args=(t_avg, argW1, argW2, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg))



		IR_p0 = np.array(p0IR)

		IR_walker_p0 = np.random.normal(IR_p0, np.abs(IR_p0)*1E-3, size=(nwalkers, ndim))

					
		clen = 512#4048
		IR_pos,_,_ = IR_sampler.run_mcmc(IR_walker_p0 , clen)



		print "SAVING THE PICKLE mmmmm"
		with open("emcee_data/Pickles/TDE_IR_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((IR_sampler.chain, IR_sampler.lnprobability), f1)


				


		### OPEN OUTPUT DATA
		with open("emcee_data/Pickles/TDE_IR_%iwalkers.pickle" %clen) as f1:
			IR_chain,IR_lnprobs = pickle.load(f1)





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


		for i,name in enumerate(param_names):
			IR_diff_minus = IR_MAP_vals[i] - IR_perc[0,i]
			IR_diff_plus = IR_perc[1,i] - IR_MAP_vals[i]
			target.write("TDE_IR: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(IR_MAP_vals[i], IR_diff_plus, IR_diff_minus, name=name))
			target.write("\n")





		IR_mxprbs = zeros(nwalkers)

					

		for i in range(nwalkers):
			IR_mxprbs[i] = max(IR_lnprobs[i])


		chi2_pdf_IR = -max(IR_mxprbs)/(len(t_avg) - len(param_names) - 1)
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
	tt = np.linspace(0.00, 12.,       Nt)*tfb

	if (Fit==False):
		#IR_p_opt = [etaR, np.cos(thetTst), np.sin(JJt), 100.0]
		if (Rfxd):
			IR_p_opt = [0.8, np.cos(thetTst), np.sin(JJt), nu0/numicron, 1.0]
		else:
			IR_p_opt = [etaR, np.cos(thetTst), np.sin(JJt), nu0/numicron, 1.0]
		V_p_opt = [0.7, t0/yr2sec, tfb/yr2sec, 2.3]
	if (Fit_IR==False):
		#IR_p_opt = [etaR, np.cos(thetTst), np.sin(JJt), 100.0]	
		if (Rfxd):
			IR_p_opt = [0.8, np.cos(thetTst), np.sin(JJt), nu0/numicron, 1.0]
		else:
			IR_p_opt = [etaR, np.cos(thetTst), np.sin(JJt), nu0/numicron, 1.0]
	if (Fit_Src==False):
		V_p_opt = [0.7, t0/yr2sec, tfb/yr2sec, 2.3]

	FsrcI1 = np.zeros(Nt)
	FVplus = np.zeros(Nt)
	FI1 = np.zeros(Nt)
	FI2 = np.zeros(Nt)


	#FsrcI1 = -2.5*np.log10(Fsrc_Anl(tt, Dst, Lav, tfb)/FVbndRel)
						   #Fsrc_Anl_Fit(t, r,      Lavg,       tfb,       t0,       gam,         FQfac)
	LIRfit = IR_p_opt[4]
	FsrcI1 = VLC_point(V_p_opt, tt, Varg, LIRfit)

	for i in range(Nt):
		#FsrcI1[i] = -2.5*np.log10(Fsrc_Anl_Fit(tt[i], Dst, V_p_opt[0], V_p_opt[2], V_p_opt[1], V_p_opt[3], V_p_opt[4]))
		if (Rfxd):
			FI1[i]    = min(IRLC_fxdR_point(IR_p_opt, tt[i], argW1, RHS_table, T_table), 12.9)
			FI2[i]    = min(IRLC_fxdR_point(IR_p_opt, tt[i], argW2, RHS_table, T_table), 11.26)
			FVplus[i] = IRLC_fxdR_point(IR_p_opt, tt[i], argVplus, RHS_table, T_table)
		else:
			FI1[i]    = min(IRLC_point(IR_p_opt, tt[i], argW1, RHS_table, T_table), 12.9)
			FI2[i]    = min(IRLC_point(IR_p_opt, tt[i], argW2, RHS_table, T_table), 11.26)
			FVplus[i] = IRLC_point(IR_p_opt, tt[i], argVplus, RHS_table, T_table)


	FVtot = -2.5*np.log10(10.**(-FsrcI1/2.5) + 10.**(-FVplus/2.5))


	###PLOT###
	plt.figure()

	if (Rfxd):
		plt.title(r"$R = %g pc$" %IR_p_opt[0])
	else:
		R0s = 0.5*IR_p_opt[0]
		plt.title(r"$R = %g L^{1/2}_{46} pc$" %R0s)

	st = plt.plot(tt/(tfb)*365., FVtot-3.5, linestyle = '-', color='blue', linewidth=2)
	s1 = plt.plot(tt/(tfb)*365., FsrcI1-3.5, linestyle = '--', color='black', linewidth=3)


	Vdat   = plt.errorbar(tV_srt, V_srt-3.5, yerr=V_sigsrt, linestyle="none", color='blue', alpha=0.2, elinewidth=1.5)
	Vsct   = plt.scatter(tV_srt, V_srt-3.5, color='blue', alpha=0.2)

	Vavg     = plt.errorbar(tV_avg, V_avg-3.5, yerr=V_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
	Vavsct   = plt.scatter(tV_avg, V_avg-3.5, color='black', alpha=1.)


	#Vav   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)



	IR1 = plt.plot(tt/(tfb)*365., FI1, color='orange', linewidth=3)#color='#1b9e77', linewidth=3)

	W1dat   = plt.errorbar(t_MJD, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=0.5, elinewidth=1.5)
	W1sct   = plt.scatter(t_MJD, W1_mag,   color='orange', alpha=0.5)

	W1av   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1.0, elinewidth=1.5)
	W1as   = plt.scatter(t_avg, W1_avg,   color='black', alpha=1.0)



	IR2 = plt.plot(tt/(tfb)*365., FI2, color='red', linewidth=3)#color='#d95f02', linewidth=3)

	W2dat   = plt.errorbar(t_MJD, W2_mag, yerr=W2_sig, linestyle="none", color='red', alpha=0.5, elinewidth=1.5)
	W2sct   = plt.scatter(t_MJD, W2_mag,  color='red', alpha=0.5)

	W2av   = plt.errorbar(t_avg, W2_avg, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
	W2as   = plt.scatter(t_avg, W2_avg,   color='black', alpha=1.)



	plt.grid(b=True, which='both')
	#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$R_d = R_0$',   r'$R_d = 0.8R_0$',   r'$R_d = 1.\bar{33}R_0$'), loc='upper right', fontsize=18)

	plt.xlabel(r"$t/t_{\rm{fb}}$")
	plt.ylabel("mag")
	#plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

	plt.ylim(plt.ylim(10.2, 14.5)[::-1])

	plt.tight_layout()

	if (Fit):
		if (Rfxd):
			Savename = "plots/TDE_Rfxd_AnalySrc_BestFits_tfb%g_clen%g_%gwalkers.png" %(tfb,clen,nwalkers)
		else:
			Savename = "plots/TDE_Rsubl_AnalySrc_BestFits_tfb%g_clen%g_%gwalkers.png" %(tfb,clen,nwalkers)
	else:
		if (Rfxd):
			Savename = "plots/TDE_Rfxd_AnalySrc_BestFits_tfb%g.png" %tfb
		else:
			Savename = "plots/TDE_Rsubl_AnalySrc_BestFits_tfb%g.png" %tfb
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)
	#plt.show()





