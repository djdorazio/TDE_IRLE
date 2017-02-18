import matplotlib
#matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 20})


from FluxFuncs_TDEs import *






###OPTIONS###OPTIONS
#All_nu = False



## SOURCE PROPERTIES
Lav = 1.*10.**45 #erg/s
Dst = 0.513*10**9*pc2cm
tfb = 1.0 * yr2sec
FQfac = 1.7  ## background flux
t0 = 2.0*yr2sec #500./365. * yr2sec


## TEST VALUES
### DUST stuff
nne = 1.6       ##emission/absortption efficeincy power law  --  Draine or 1.6 
nu0 = numicron/4.  ## emission/absortption efficeincy cutoff freq
Rde = 0.3*pc2cm #* np.sqrt(Lav/10.**46) * pc2cm     ## inner edge of dust
thetTst = 0.0   ## opening angle for dust torus (0 is sphere, pi/2 is ring)
JJt = 0.0       ## inclination angle of dust torus (0 is edge on pi/2 is face on)
aeff = (c/nu0)/(2.*ma.pi) #set dust grain size by cutoff frequency


###OPT THIN DUST
##For dust which is optically thin to opt/UV ##IGNORE FOR NOW
Rrout = 1.0*Rde ## oter edge of dust
pp = 2.0 	    ## dust density power law
n10 = 1.0/(ma.pi*Rde*aeff*aeff) * (pp-1.)
nfac = 10.0 
n0 = nfac*n10 ##6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))






# ### WISE BAND + Observational STUFF
## edges of WISE bands
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3


## zero point fluxes - check these...
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


### PLOT POINTS
Nt=20
tt = np.linspace(0.00, 6.,       Nt)*tfb


## INTEGRATION LIMTS FOR ALL nu
# if (All_nu):
Nnumn = 0.0
Nnumx = 5.0
numn = Nnumn*numicron
numx = Nnumx*numicron
# elif (WISE_nu):
numn1 = W1mn
numx1 = W1mx
numn2 = W2mn
numx2 = W2mx
# else:
# 	print "Define limits of integration over nu"


####################################################
###   Plot echoes in W1 and W2 bands
####################################################


arg1 = [Lav, tfb, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne, FQfac, t0]

FsrcI1 = np.zeros(Nt)
FI1 = np.zeros(Nt)
FI2 = np.zeros(Nt)


#FsrcI1 = -2.5*np.log10(Fsrc_Anl(tt, Dst, Lav, tfb)/FVbndRel)

for i in range(Nt):
	FsrcI1[i] = -2.5*np.log10(Fsrc_Anl(tt[i], Dst, Lav, tfb, t0, 10.0)/10/FVbndRel)
	FI1[i]    = -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt(numn1, numx1, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
	FI2[i] = -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt(numn2, numx2, tt[i], Dst, arg1, RHS_table, T_table)/FW2Rel)

nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)

###PLOT###
plt.figure()

#plt.title(r"Isotropic, Sphere, $P =  R_0/c$")

s1 = plt.plot(tt/(tfb), FsrcI1, linestyle = '--', color='blue', linewidth=3)

IR1 = plt.plot(tt/(tfb), FI1+nrm, color='orange', linewidth=3)#color='#1b9e77', linewidth=3)

IR2 = plt.plot(tt/(tfb), FI2+nrm, color='red', linewidth=3)#color='#d95f02', linewidth=3)


plt.grid(b=True, which='both')
#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$R_d = R_0$',   r'$R_d = 0.8R_0$',   r'$R_d = 1.\bar{33}R_0$'), loc='upper right', fontsize=18)

plt.xlabel(r"$t/t_{\rm{fb}}$")
plt.ylabel("mag")
#plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

plt.ylim(plt.ylim(10.2, 20.)[::-1])

plt.tight_layout()

Savename = "plots/TDE_analytic_tfb%g.png" %tfb
Savename = Savename.replace('.', 'p')
Savename = Savename.replace('ppng', '.png')
plt.savefig(Savename)
#plt.show()



#0.003 cnts/sec in 2003 Chandra 0.2kev  3* 10^39 erg/s

#0.008 in 2009


		
