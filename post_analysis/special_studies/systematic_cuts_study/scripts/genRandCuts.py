import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd


#output file to write cuts file
ofname = 'cafe_systematics_cuts_file_noQ2var_10k.csv' 
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
            '# hms/shms_coll   : HMS/SHMS collimator scale factor cut variation\n'
            '#'
            )
ofile.write('# entry,Q2_min,Q2_max,Em_min_mf,Em_max_mf,Pm_min_mf,Pm_max_mf,xbj_min,xbj_max,thrq_min,thrq_max,Pm_min_src,Pm_max_src, hms_coll, shms_coll\n') 


#----------------------------------------------------------------
# Set standard CaFe cuts (mu) +/-  2sigma to randomly generate
# a value around the about the central cut

# Q2 (same random cut for MF, SRC) : 1.8 +/- 0.05 GeV^2
Q2_min_mu = 1.8
Q2_min_sig = 0.0 # 0.05

# Emiss lower bound (MF) :  -20 MeV (fixed)
Em_min_mu_mf = -0.02

# Emiss upper bound (MF) : 90 +/- 2.5 MeV
Em_max_mu_mf = 0.09
Em_max_sig_mf = 0.0025

# cut variation only needed for Pmax_mf (MF) : 270 +/- 10 MeV
Pm_max_mf_mu = 0.270
Pm_max_mf_sig = 0.01

# hms/shms collimator (vary scalee factor by sig=4%
hms_coll_mu = 1    # scale factor of 1 (geometry cut on collimator)
hms_coll_sig = 0.04  # 4% variation in the scaler factor 

shms_coll_mu = 1    # scale factor of 1 (geometry cut on collimator)
shms_coll_sig = 0.04  # 4% variation in the scaler factor 

# x-Bjorken (SRC) : 1.2 +/- 0.05
xbj_min_mu = 1.2
xbj_min_sig = 0.05

# th_rq (SRC) : 40 +/- 2 deg
thrq_max_mu = 40
thrq_max_sig = 2.

# cut variation for  min Pm (SRC): 375 +/- 12.5 MeV/c
Pm_min_src_mu = 0.375
Pm_min_src_sig = 0.0125

# cut variation for max Pm (SRC): 700 +/- 50 MeV/c 
Pm_max_src_mu = 0.70
Pm_max_src_sig = 0.05

# ---------------------------------------------------------------


# generate N random entries of cut ranges extracted from a gaussian sample
entries = 10000

# get random cut from randoom gaussian sample
hms_coll       = np.random.normal(hms_coll_mu, hms_coll_sig, entries)
shms_coll      = np.random.normal(shms_coll_mu, shms_coll_sig, entries)

Q2_min       = np.random.normal(Q2_min_mu,    Q2_min_sig, entries)
Q2_max       = 100.  # fixed upper limit

Em_min_mf    = Em_min_mu_mf  # fixed lower limit 
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
    ofile.write("%i,%.3f, %.3f, %.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f \n" % (i, Q2_min[i], Q2_max,Em_min_mf,Em_max_mf[i],Pm_min_mf,Pm_max_mf[i],xbj_min_src[i],xbj_max_src,thrq_min_src,thrq_max_src[i],Pm_min_src[i],Pm_max_src[i], hms_coll[i], shms_coll[i] ))
ofile.close()


#--------------------------------------------------
fig, axs = plt.subplots(nrows=3, ncols=3)
nbins = 20

axs[0,0].hist(Q2_min, nbins, density=True, histtype='stepfilled', facecolor='violet', alpha=0.75)
axs[0,0].set_title(r'4-Momentum Transfer, $Q^{2}_{min}$')
axs[0,0].set_xlabel(r'$Q^{2}_{min}$ [GeV$^{2}$]')


axs[0, 1].hist(Em_max_mf, nbins, density=True, histtype='stepfilled', facecolor='dodgerblue', alpha=0.75)
axs[0, 1].set_title(r'Missing Energy, $E^{MF}_{miss, max}$')
axs[0, 1].set_xlabel(r'$E_{miss}$ [GeV]')

axs[0, 2].hist(Pm_max_mf, nbins, density=True, histtype='stepfilled', facecolor='tomato', alpha=0.75)
axs[0, 2].set_title(r'Missing Momentum, $P^{MF}_{miss, max}$')
axs[0, 2].set_xlabel(r'$P_{miss, max}$ [GeV/c]')


axs[1, 0].hist(xbj_min_src, nbins, density=True, histtype='stepfilled', facecolor='gold', alpha=0.75)
axs[1, 0].set_title(r'x-Bjorken, $x_{Bj,min}$')
axs[1, 0].set_xlabel(r'$x_{Bj, min}$')

axs[1, 1].hist(thrq_max_src, nbins, density=True, histtype='stepfilled', facecolor='crimson', alpha=0.75)
axs[1, 1].set_title(r'recoil angle, $\theta_{rq,min}$')
axs[1, 1].set_xlabel(r'$\theta_{rq,min} [deg]$')

axs[1, 2].hist(Pm_min_src, nbins, density=True, histtype='stepfilled', facecolor='lime', alpha=0.75)
axs[1, 2].set_title(r'Missing Momentum, $P^{SRC}_{miss, min}$')
axs[1, 2].set_xlabel(r'$P^{SRC}_{miss, min}$ [GeV/c]')

axs[2, 0].hist(Pm_max_src, nbins, density=True, histtype='stepfilled', facecolor='forestgreen', alpha=0.75)
axs[2, 0].set_title(r'Missing Momentum, $P^{SRC}_{miss, max}$')
axs[2, 0].set_xlabel(r'$P^{SRC}_{miss, max}$ [GeV/c]')

axs[2, 1].hist(hms_coll, nbins, density=True, histtype='stepfilled', facecolor='darkorange', alpha=0.75)
axs[2, 1].set_title(r'HMS Collimator')
axs[2, 1].set_xlabel(r'Scale Factor')

axs[2, 2].hist(shms_coll, nbins, density=True, histtype='stepfilled', facecolor='plum', alpha=0.75)
axs[2, 2].set_title(r'SHMS Collimator')
axs[2, 2].set_xlabel(r'Scale Factor')


fig.tight_layout()
plt.show()


'''
#---------------------------------------

# this part of the code reads directly the input cuts file to make the plots
# NOTE: This requires the column header to be read without the '#' symbol (without comments)
# in order to be able to read the columns using pandas dataframe

# alternatively read the input cuts file
#fname='../input/cafe_systematics_cuts_file_noQ2var_10k.csv'
fname='./cafe_systematics_cuts_file_noQ2var_10k.csv'

df = pd.read_csv(fname, comment='#')
df.to_numpy()
    
# prepare subplots for histogramming random cuts
fig, axs = plt.subplots(nrows=3, ncols=3)

nbins = 20

axs[0,0].hist(df.Q2_min, nbins, density=True, histtype='stepfilled', facecolor='violet', alpha=0.75)
axs[0,0].set_title(r'4-Momentum Transfer, $Q^{2}_{min}$')
axs[0,0].set_xlabel(r'$Q^{2}_{min}$ [GeV$^{2}$]')


axs[0, 1].hist(df.Em_max_mf, nbins, density=True, histtype='stepfilled', facecolor='dodgerblue', alpha=0.75)
axs[0, 1].set_title(r'Missing Energy, $E^{MF}_{miss, max}$')
axs[0, 1].set_xlabel(r'$E_{miss}$ [GeV]')

axs[0, 2].hist(df.Pm_max_mf, nbins, density=True, histtype='stepfilled', facecolor='tomato', alpha=0.75)
axs[0, 2].set_title(r'Missing Momentum, $P^{MF}_{miss, max}$')
axs[0, 2].set_xlabel(r'$P_{miss, max}$ [GeV/c]')


axs[1, 0].hist(df.xbj_min, nbins, density=True, histtype='stepfilled', facecolor='gold', alpha=0.75)
axs[1, 0].set_title(r'x-Bjorken, $x_{Bj,min}$')
axs[1, 0].set_xlabel(r'$x_{Bj, min}$')

axs[1, 1].hist(df.thrq_max, nbins, density=True, histtype='stepfilled', facecolor='crimson', alpha=0.75)
axs[1, 1].set_title(r'recoil angle, $\theta_{rq,min}$')
axs[1, 1].set_xlabel(r'$\theta_{rq,min} [deg]$')

axs[1, 2].hist(df.Pm_min_src, nbins, density=True, histtype='stepfilled', facecolor='lime', alpha=0.75)
axs[1, 2].set_title(r'Missing Momentum, $P^{SRC}_{miss, min}$')
axs[1, 2].set_xlabel(r'$P^{SRC}_{miss, min}$ [GeV/c]')

axs[2, 0].hist(df.Pm_max_src, nbins, density=True, histtype='stepfilled', facecolor='forestgreen', alpha=0.75)
axs[2, 0].set_title(r'Missing Momentum, $P^{SRC}_{miss, max}$')
axs[2, 0].set_xlabel(r'$P^{SRC}_{miss, max}$ [GeV/c]')

axs[2, 1].hist(df.hms_coll, nbins, density=True, histtype='stepfilled', facecolor='darkorange', alpha=0.75)
axs[2, 1].set_title(r'HMS Collimator')
axs[2, 1].set_xlabel(r'Scale Factor')

axs[2, 2].hist(df.shms_coll, nbins, density=True, histtype='stepfilled', facecolor='plum', alpha=0.75)
axs[2, 2].set_title(r'SHMS Collimator')
axs[2, 2].set_xlabel(r'Scale Factor')


#--------------------------------------------------

 
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
