import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import asarray as ar


# description: this script reads the data yield
# over the systematically-varied cuts for each 
# (target, kinematics) as well as the ratios of
# SRC/MF for target A to a reference nucleus, etc.
# to study the cut sensitivity on ratios as well
# NOTE: When studying cut sensitivity on ratios,
# the ratio of yields should be taken on an entry-by-entry basis
# i.e., ratio of numerical arrays

fname='cafe_systematics_Be9_SRC.csv'
df_data = pd.read_csv(fname, comment='#')
df_data.to_numpy()


#Calculating the Gaussian PDF values given Gaussian parameters and random variable X
def gaus(X,C,X_mean,sigma):
    return C*np.exp(-(X-X_mean)**2/(2*sigma**2))



fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True)

nbins = 20

# define systematic cut arrays to loop over
syst_name=['syst_total_real', 'syst_dPm_min_real', 'syst_dPm_max_real', 'syst_dQ2_real', 'syst_dXbj_real', 'syst_dthrq_real', 'syst_dHcoll_real', 'syst_dScoll_real']
clr=['silver', 'lime', 'forestgreen', 'violet', 'gold', 'crimson', 'lime', 'plum']
title=[r'Integrated Yield (total)', r'Integrated Yield ($\delta Pm_{min}$)', r'Integrated Yield ($\delta Pm_{max}$)', r'Integrated Yield ($\delta Q^{2}$)', r'Integrated Yield ($\delta X^{2}_{bj}$)', r'Integrated Yield ($\delta\theta_{rq}$)', r'Integrated Yield ($\delta HMS_{coll}$)', r'Integrated Yield ($\delta SHMS_{coll}$)']



# loop oved subplot indices (ith-row, jth-col)
row=[0,1,2]
col=[0,1,2]
idx = 0  # index counter
for i in row:
    for j in col:

        #-------------------------
        # plot the data histogram
        #-------------------------
        # yhist:   is the number of counts in each bin of the histogram
        # xedge:   is the left hand edge of each bin
        # patches: is the individual patches used to create the histogram, e.g a collection of rectangles
        
        yhist, xedge, patches = axs[i,j].hist(df_data[syst_name[idx]], nbins, density=False, histtype='stepfilled', facecolor=clr[idx], alpha=0.75)
        xhist = (xedge[:-1] + xedge[1:])/2  # bin center
        
        #---------------------------------------
        # Gaussian least-square fitting process
        #---------------------------------------
        
        # get data info (peak, mean, sigma) to be used as input fit param
        data_peak = yhist.max()
        data_mu   = df_data[syst_name[idx]].mean()
        data_sig  = df_data[syst_name[idx]].std()
        
        # perform fit (provided fit function-> gaus, data-> xhist,yhist, and fit param p0)
        popt,pcov = curve_fit(gaus,xhist,yhist ,p0=[data_peak, data_mu, data_sig], maxfev=5000) 
        
        # get fit parameters
        mu_fit,  mu_fit_err  = popt[1], np.sqrt(pcov[1,1])
        sig_fit, sig_fit_err = popt[2], np.sqrt(pcov[2,2])

        # create evenly-spaced points for plotting fit function
        xmin=df_data[syst_name[idx]].min()
        xmax=df_data[syst_name[idx]].max()
        x_fit = np.linspace(xmin, xmax, 500) 
        
        # plot fit function and write fit parameters as label for legend
        axs[i,j].plot(x_fit,gaus(x_fit,*popt), color='r', label=r'$\mu:{0:.0f}\pm{1:.0f}$''\n''$\sigma:{2:.0f}\pm{3:.0f}$'.format(mu_fit, mu_fit_err, sig_fit, sig_fit_err))
        axs[i,j].set_title(title[idx])
        axs[i,j].legend(frameon=False, loc='upper right')

        
        idx=idx+1

        if idx==8:
            break

fig.text(0.5, 0.015, 'Integrated Yield', ha='center', fontsize=16)
fig.text(0.015, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=16)

fig.tight_layout()
plt.show()


