import numpy as np
import random
import matplotlib.pyplot as plt

#output file to write cuts file
ofname = 'cafe_systematics_cuts_file.csv' 
ofile = open(ofname, 'w+')
ofile.write('# CaFe Systematics Cuts File\n')
ofile.write('# \n'
            '# Header Definitions: \n'
            '# entry           : entry to identify each of the gaussian random-generated cuts below\n'
            '# Q2_min,max      : 4-momentum transfer (MF,SRC) [GeV^2] \n'
            '# Em_min,max_mf   : Missing Energy (MF) [GeV] \n'
            '# Pm_min,max_mf   : Missing Momentum (MF) [GeV/c] \n'
            '# xbj_min,max     : x-bjorken (SRC) \n'
            '# thrq_min,max    : recoil nucleon angle relative to q-vector (SRC) [deg] \n'
            '# Pm_min,max_src  : Missing Momentum (SRC) [GeV/c] \n'

            
            )
ofile.write('entry,Q2_min,Q2_max,Em_min_mf,Em_max_mf,Pm_min_mf,Pm_max_mf,xbj_min,xbj_max,thrq_min,thrq_max,Pm_min_src,Pm_max_src\n') 


#----------------------------------------------------------------
# Set standard CaFe cuts (mu) +/-  sigma to randomly generate
# a value around the about the central cut

# Q2 (same random cut for MF, SRC)
Q2_min_mu = 1.8
Q2_min_sig = 0.1

# Emiss lower bound (MF) :  -20 +/- 5 MeV
Em_min_mu_mf = -0.02
Em_min_sig_mf = 0.005

# Emiss upper bound (MF) : 100 +/- 10
Em_max_mu_mf = 0.1
Em_max_sig_mf = 0.01

# cut variation only needed for Pmax_mf (MF)
Pm_max_mf_mu = 250
Pm_max_mf_sig = 10.

# x-Bjorken (SRC) : 1.2 +/- 0.1
xbj_min_mu = 1.2
xbj_min_sig = 0.1

# th_rq (SRC) : 40 +/- 3
thrq_max_mu = 40
thrq_max_sig = 3.

# should a systematic variation be made on Pmiss ?



# cut variation for both min/max Pm (SRC)
Pm_min_src_mu = 300
Pm_min_src_sig = 10.

Pm_max_src_mu = 700
Pm_max_src_sig = 10.

# ---------------------------------------------------------------


# generate N random entries of cut ranges extracted from a gaussian sample
entries = 1000

# get random cut from randoom gaussian sample
Q2_min       = np.random.normal(Q2_min_mu,    Q2_min_sig, entries)
Q2_max       = 100.  # fixed upper limit

Em_min_mf    = np.random.normal(Em_min_mu_mf, Em_min_sig_mf, entries)
Em_max_mf    = np.random.normal(Em_max_mu_mf, Em_max_sig_mf, entries)

xbj_min_src  = np.random.normal(xbj_min_mu,   xbj_min_sig, entries)
xbj_max_src  = 100.

thrq_min_src = 0.
thrq_max_src = np.random.normal(thrq_max_mu,  thrq_max_sig, entries)

Pm_min_mf    = 0.
Pm_max_mf    = np.random.normal(Pm_max_mf_mu, Pm_max_mf_sig, entries)

Pm_min_src   = np.random.normal(Pm_min_src_mu, Pm_min_src_sig, entries)
Pm_max_src   = np.random.normal(Pm_max_src_mu, Pm_max_src_sig, entries)

for i in range(entries):
    # write random cuts to a .txt file
    ofile.write("%i,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f \n" % (i, Q2_min[i], Q2_max,Em_min_mf[i],Em_max_mf[i],Pm_min_mf,Pm_max_mf[i],xbj_min_src[i],xbj_max_src,thrq_min_src,thrq_max_src[i],Pm_min_src[i],Pm_max_src[i] ))
ofile.close()

#--------------------------------------------------
# prepare subplots for histogramming random cuts
fig, axs = plt.subplots(nrows=2, ncols=4)

axs[0, 0].hist(Q2_min, 100, density=True, histtype='stepfilled', facecolor='violet', alpha=0.75)
axs[0,0].set_title(r'4-Momentum Transfer, $Q^{2}_{min}$')
axs[0,0].set_xlabel(r'$Q^{2}_{min}$ [GeV$^{2}$]')

axs[0, 1].hist(Em_min_mf, 100, density=True, histtype='stepfilled', facecolor='lightskyblue', alpha=0.75)
axs[0, 1].set_title(r'Missing Energy, $E^{MF}_{miss, min}$')
axs[0, 1].set_xlabel(r'$E_{m}$ [GeV]')


axs[0, 2].hist(Em_max_mf, 100, density=True, histtype='stepfilled', facecolor='dodgerblue', alpha=0.75)
axs[0, 2].set_title(r'Missing Energy, $E^{MF}_{miss, max}$')
axs[0, 2].set_xlabel(r'$E_{miss}$ [GeV]')

axs[0, 3].hist(Pm_max_mf, 100, density=True, histtype='stepfilled', facecolor='tomato', alpha=0.75)
axs[0, 3].set_title(r'Missing Momentum, $P^{MF}_{miss, max}$')
axs[0, 3].set_xlabel(r'$P_{miss, max}$ [GeV/c]')


axs[1, 0].hist(xbj_min_src, 100, density=True, histtype='stepfilled', facecolor='gold', alpha=0.75)
axs[1, 0].set_title(r'x-Bjorken, $x_{Bj,min}$')
axs[1, 0].set_xlabel(r'$x_{Bj, min}$')

axs[1, 1].hist(thrq_max_src, 100, density=True, histtype='stepfilled', facecolor='crimson', alpha=0.75)
axs[1, 1].set_title(r'recoil angle, $\theta_{rq,min}$')
axs[1, 1].set_xlabel(r'$\theta_{rq,min} [deg]$')

axs[1, 2].hist(Pm_min_src, 100, density=True, histtype='stepfilled', facecolor='lime', alpha=0.75)
axs[1, 2].set_title(r'Missing Momentum, $P^{SRC}_{miss, min}$')
axs[1, 2].set_xlabel(r'$P^{SRC}_{miss, min}$ [GeV/c]')

axs[1, 3].hist(Pm_max_src, 100, density=True, histtype='stepfilled', facecolor='forestgreen', alpha=0.75)
axs[1, 3].set_title(r'Missing Momentum, $P^{SRC}_{miss, max}$')
axs[1, 3].set_xlabel(r'$P^{SRC}_{miss, max}$ [GeV/c]')

fig.tight_layout()
plt.show()
#--------------------------------------------------

 
'''
# loop over all entries
for ientry in range(entries):

    # get random cut from randoom gaussian sample
    Q2_min       = random.gauss(Q2_min_mu,    Q2_min_sig)
    Q2_max       = 100.  # fixed upper limit

    Em_min_mf    = random.gauss(Em_min_mu_mf, Em_min_sig_mf)
    Em_max_mf    = random.gauss(Em_max_mu_mf, Em_max_sig_mf)

    xbj_min_src  = random.gauss(xbj_min_mu,   xbj_min_sig)
    xbj_max_src  = 100.

    thrq_min_src = 0.
    thrq_max_src = random.gauss(thrq_max_mu,  thrq_max_sig)

    Pm_min_mf    = 0.
    Pm_max_mf    = random.gauss(Pm_max_mf_mu, Pm_max_mf_sig)

    Pm_min_src   = random.gauss(Pm_min_src_mu, Pm_min_src_sig)
    Pm_max_src   = random.gauss(Pm_max_src_mu, Pm_max_src_sig)
    
    debug=False

    if(debug):
        
        print(
            '-----------------------------------------------\n',
            'Analysis Cuts Randomly-Generated from Gaussian \n',
            '-----------------------------------------------\n',
            '\n',
            'entry: %i\n' % ientry, 
            'mean-field (MF): \n',
            'Q2_min           : %.3f GeV^2 \n' % (Q2_min),
            '(Em_min, Em_max) : (%.3f, %.3f) GeV \n' % (Em_min_mf, Em_max_mf),
            '(Pm_min, Pm_max) : (0, %.3f) GeV/c \n' % (Pm_max_mf),
            '\n',
            'short-range correlations (SRC): \n',
            'Q2_min               : %.3f GeV^2 \n' % (Q2_min),
            'xbj_min              : %.3f \n' % (xbj_min_src),
            '(thrq_min, thrq_max) : (0, %.3f) deg \n' % (thrq_max_src),
            '(Pm_min, Pm_max)     : (%.3f, %.3f) GeV/c \n' % (Pm_min_src, Pm_max_src),
            
        )


    # write random cuts to a .txt file
    ofile.write("%i,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f \n" % (ientry, Q2_min, Q2_max,Em_min_mf,Em_max_mf,Pm_min_mf,Pm_max_mf,xbj_min_src,xbj_max_src,thrq_min_src,thrq_max_src,Pm_min_src,Pm_max_src ))
ofile.close()
'''
