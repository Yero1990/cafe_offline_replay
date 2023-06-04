import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy

    
"""
---------------------------
Author: C. Yero
Date: June 04, 2023
email: cyero@jlab.org,
       cyero002@gmail.com
---------------------------

-------------
description:
-------------
this script does the following:
1) reads the data yield over the systematically-varied cuts for each (target, kinematics)
2) histograms the individual (target,kin) set, single ratios or double ratios 
3) fits the histograms with a gaussian-distribution to determine mean,standard_deviation
   where the standard_deviation / mean is a direct measure of the systematic spread
"""


def printHelp():
    print('''\n
    ---------------------------------------------------------------------------------------------------------------
    Usage                                                              Description                       
    ---------------------------------------------------------------------------------------------------------------

    ipython fit_systematics.py <tgt> <kin>                             only calculate systematics for
    e.g., ipython fit_systematics.py Be9 MF                            single (target,kin) configuration
                             
    ipython fit_systematics.py single <tgt1> <kin1> <tgt2> <kin2>      calculate systematics for a singles ratio: 
    e.g., ipython fit_systematics.py single Be9 SRC Be9 MF             R = (target1,kin1) / (target2,kin2)
    
    ipython fit_systematics.py double <tgt1>  <tgt2>                   calculate systematics for a double ratio: 
    e.g., ipython fit_systematics.py double Be9 C12                    R = (target1_SRC/MF) / (target2_SRC/MF)

    tgt: Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
    kiN: MF, SRC
    ---------------------------------------------------------------------------------------------------------------

    ''')



# define gaussian fit function 
def gaus(X,C,X_mean,sigma):
    return C*np.exp(-(X-X_mean)**2/(2*sigma**2))



# define functions to get normalization quantities from .csv file
def get_var(tgt, kin, var, npass):
    csv_file='../../../summary_files/%s/cafe_prod_%s_%s_report_summary.csv' % (npass,tgt,kin)

    # read .csv file
    df           = pd.read_csv(csv_file, comment='#')
    charge       = df['charge'] # [mC]
    hms_trk_eff  = unumpy.uarray(df['hTrkEff'],         df['hTrkEff_err'])
    shms_trk_eff = unumpy.uarray(df['pTrkEff'],         df['pTrkEff_err'])
    total_LT     = unumpy.uarray(df['tLT'],             df['tLT_err_Bi'])
    mult_trk_eff = np.array(df['multi_track_eff'])  

    # calculate average hms/shms track efficiency and live time
    # (these may be useful later on, for normalizing yield)
    avg_hms_trk_eff  =  hms_trk_eff.mean()
    avg_shms_trk_eff =  shms_trk_eff.mean()
    avg_total_LT     =  total_LT.mean()
    avg_mult_trk_eff =  mult_trk_eff.mean()
    # sum over all total charge
    total_charge = charge.sum()

    if(var=='charge'):
        return total_charge
    if(var=='hms_trk_eff'):
        return  unumpy.nominal_values(avg_hms_trk_eff)
    if(var=='shms_trk_eff'):
        return  unumpy.nominal_values(avg_shms_trk_eff)
    if(var=='total_LT'):
        return  unumpy.nominal_values(avg_total_LT)
    if(var=='multi_trk_eff'):
        return  avg_mult_trk_eff


    
# Usage: ipython fit_systematics.py  Be9 MF
if (len(sys.argv) == 3):
    
       
    target = sys.argv[1]    # Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
    kin    = sys.argv[2]    # MF, SRC

    
    
    fname='../output/cafe_systematics_%s_%s.csv' % (target, kin)
    
    df_data = pd.read_csv(fname, comment='#')
    df_data.to_numpy()

    

    if kin=="SRC":
        fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True)
        fig.set_size_inches(10,8, forward=True)
        plt.subplots_adjust(bottom=0.75)
        plt.subplots_adjust(left=0.56)
        fig.text(0.5, 0.001, 'Integrated Yield', ha='center', fontsize=14)
        fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=14)

    elif kin=="MF":
        fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True)
        fig.set_size_inches(12,8, forward=True)
        plt.subplots_adjust(bottom=0.8)
        plt.subplots_adjust(left=0.3)
        plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
        fig.text(0.5, 0.01, 'Integrated Yield', ha='center', fontsize=14)
        fig.text(0.001, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=14)

    nbins = 20

    # define systematic cut arrays to loop over
    if kin=="SRC":
        syst_name=['syst_total_real', 'syst_dPm_min_real', 'syst_dPm_max_real', 'syst_dQ2_real', 'syst_dXbj_real', 'syst_dthrq_real', 'syst_dHcoll_real', 'syst_dScoll_real']
        clr=['silver', 'lime', 'forestgreen', 'violet', 'gold', 'crimson', 'lime', 'plum']
        title=[r'(total)', r'($\delta Pm_{min}$)', r'($\delta Pm_{max}$)', r'($\delta Q^{2}$)', r'($\delta X^{2}_{bj}$)', r'($\delta\theta_{rq}$)', r'($\delta HMS_{coll}$)', r'($\delta SHMS_{coll}$)']

    elif kin=="MF":
        syst_name=['syst_total_real','syst_dPm_max_real','syst_dQ2_real','syst_dEm_real','syst_dHcoll_real','syst_dScoll_real']
        clr=['silver', 'tomato', 'violet', 'dodgerblue', 'darkorange', 'plum']
        title=[r'(total)', r'($\delta Pm_{max}$)',r'($\delta Q^{2}$)', r'($\delta Em_{max}$)', r'($\delta HMS_{coll}$)', r'($\delta SHMS_{coll}$)']

    
    # loop oved subplot indices (ith-row, jth-col)
    if kin=="SRC":
        row=[0,1,2]
        col=[0,1,2]
    elif kin=="MF":
        row=[0,1]
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
            
            print('(idx, row, col): %i %i %i'%(idx,i,j))
            print(syst_name[idx])


        
            
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

            rel_err = sig_fit / mu_fit
             
            print('pcov[1,1]:', np.sqrt(pcov[1,1]))
            print('pcov[2,2]:', np.sqrt(pcov[2,2]))
            # create evenly-spaced points for plotting fit function
            xmin=df_data[syst_name[idx]].min()
            xmax=df_data[syst_name[idx]].max()
            x_fit = np.linspace(xmin, xmax, 500) 

            bin_w = xedge[1] - xedge[0]
            if(bin_w < 1.):
                mu_fit_err = 0
                sig_fit_err = 0
            # plot fit function and write fit parameters as label for legend
            axs[i,j].plot(x_fit,gaus(x_fit,*popt), color='r', label=r'$\mu:{0:.0f}\pm{1:.0f}$''\n''$\sigma:{2:.0f}\pm{3:.0f}$''\n''$\sigma$ / $\mu:{4:.3f}$'.format(mu_fit, mu_fit_err, sig_fit, sig_fit_err, rel_err))
            axs[i,j].set_title(title[idx])
            axs[i,j].legend(frameon=False, loc='upper right')
            axs[i,j].tick_params(labelbottom=True)

            
            idx=idx+1

            if(kin=="SRC" and idx==8):
                break
            if(kin=="MF" and idx==6):
                break  
        
    # common title
    fig.suptitle('%s %s Integrated Yield'%(target,kin), fontsize=20)
    fig.tight_layout()
    plt.show()
    #plt.savefig('%s_%s_systematics.pdf'%(target,kin))


# ------------------------------------
# CALCULATE SINGLE RATIO SYSTEMATICS
# ------------------------------------

# Usage: ipython fit_systematics.py single target1 kin1 target2 kin2
# to calculate the ratio:  R = target1_kin1 / target2_kin2
elif (len(sys.argv)==6 and sys.argv[1] == "single"):

    target1 = sys.argv[2]    # Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
    kin1    = sys.argv[3]    # MF, SRC

    target2 = sys.argv[4]    # Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
    kin2    = sys.argv[5]    # MF, SRC
    
    fname1='../output/cafe_systematics_%s_%s.csv' % (target1, kin1)
    fname2='../output/cafe_systematics_%s_%s.csv' % (target2, kin2)

    print('\n reading files: \n --> %s \n --> %s' % (fname1, fname2))
    print('\n calculating single ratio systematics: %s %s / %s %s' % (target1, kin1, target2, kin2))
    
    df_data1 = pd.read_csv(fname1, comment='#')
    df_data1.to_numpy()

    df_data2 = pd.read_csv(fname2, comment='#')
    df_data2.to_numpy()

    # calculate the ratios on an entry-by-entry basis
    # where each entry is a specific combination of cuts selected from
    # cuts file of N cuts. The cuts file has been generated
    # from randomly-sampled gaussian distributions  about the central value
    # for each of the kinematical cuts
    
    R_single  = df_data1['syst_total_real']  /  df_data2['syst_total_real'] 


    # set up the canvas
    fig, axs = plt.subplots()
    fig.set_size_inches(10,8, forward=True)

    nbins = 20

    yhist, xedge, patches = axs.hist(R_single, nbins, density=False, histtype='stepfilled', facecolor='silver', alpha=0.75)
    xhist = (xedge[:-1] + xedge[1:])/2  # bin center
    
    #---------------------------------------
    # Gaussian least-square fitting process
    #---------------------------------------
    
    # get data info (peak, mean, sigma) to be used as input fit param
    data_peak = yhist.max()
    data_mu   = R_single.mean()
    data_sig  = R_single.std()

  
    # perform fit (provided fit function-> gaus, data-> xhist,yhist, and fit param p0)
    popt,pcov = curve_fit(gaus,xhist,yhist ,p0=[data_peak, data_mu, data_sig], maxfev=5000) 

    # get fit parameters
    mu_fit,  mu_fit_err  = popt[1], np.sqrt(pcov[1,1])
    sig_fit, sig_fit_err = popt[2], np.sqrt(pcov[2,2])

    rel_err = sig_fit / mu_fit
    
    # create evenly-spaced points for plotting fit function
    xmin=R_single.min()
    xmax=R_single.max()
    x_fit = np.linspace(xmin, xmax, 500) 
    
    # plot fit function and write fit parameters as label for legend
    axs.plot(x_fit,gaus(x_fit,*popt), color='r', label=r'$\mu:{0:.3f}\pm{1:.3f}$''\n''$\sigma:{2:.3f}\pm{3:.3f}$''\n''$\sigma$ / $\mu:{4:.3f}$'.format(mu_fit, mu_fit_err, sig_fit, sig_fit_err, rel_err))
    axs.legend(frameon=False, loc='upper right')
    
            
            
    #fig.tight_layout()
    if(target1 != target2):
        plt.title('Single Ratio %s %s / %s %s (systematics)' % (target1,kin1,target2,kin2), fontsize=16)
    elif(target1==target2):
        plt.title('Single Ratio %s %s / %s (systematics)' % (target1,kin1,kin2), fontsize=16)        
    plt.xlabel('Single Ratio Counts', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
        
    plt.show()
    #plt.savefig('%s_%s_single_ratio_systematics.pdf'%(target,kin))




# ------------------------------------
# CALCULATE DOUBLE RATIO SYSTEMATICS
# ------------------------------------

# Usage: ipython fit_systematics.py double target1 target2 
# to calculate the ratio:  R = (target1_SRC/target1_MF) / )target2_SRC/target2_MF)
elif (len(sys.argv)==4 and sys.argv[1] == "double"):

    target1 = sys.argv[2]    # Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
    target2 = sys.argv[3]    # Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
    
    fname1_SRC='../output/cafe_systematics_%s_SRC.csv' % (target1)
    fname1_MF='../output/cafe_systematics_%s_MF.csv' % (target1)

    fname2_SRC='../output/cafe_systematics_%s_SRC.csv' % (target2)
    fname2_MF='../output/cafe_systematics_%s_MF.csv' % (target2)


    print('\n reading files: \n --> %s \n --> %s  \n --> %s \n --> %s' % (fname1_SRC, fname1_MF, fname2_SRC, fname2_MF))
    print('\n calculating double ratio systematics: %s / %s ' % (target1, target2))

    df_data1_SRC = pd.read_csv(fname1_SRC, comment='#')
    df_data1_SRC.to_numpy()
    
    df_data1_MF = pd.read_csv(fname1_MF, comment='#')
    df_data1_MF.to_numpy()

    df_data2_SRC = pd.read_csv(fname2_SRC, comment='#')
    df_data2_SRC.to_numpy()
    
    df_data2_MF = pd.read_csv(fname2_MF, comment='#')
    df_data2_MF.to_numpy()  

    # calculate the ratios on an entry-by-entry basis
    # where each entry is a specific combination of cuts selected from
    # cuts file of N cuts. The cuts file has been generated
    # from randomly-sampled gaussian distributions  about the central value
    # for each of the kinematical cuts
    
    R_double  = (df_data1_SRC['syst_total_real'] /  df_data1_MF['syst_total_real']) /  (df_data2_SRC['syst_total_real'] /  df_data2_MF['syst_total_real'])


    # set up the canvas
    fig, axs = plt.subplots()
    fig.set_size_inches(10,8, forward=True)

    nbins = 20

    yhist, xedge, patches = axs.hist(R_double, nbins, density=False, histtype='stepfilled', facecolor='silver', alpha=0.75)
    xhist = (xedge[:-1] + xedge[1:])/2  # bin center
    
    #---------------------------------------
    # Gaussian least-square fitting process
    #---------------------------------------
    
    # get data info (peak, mean, sigma) to be used as input fit param
    data_peak = yhist.max()
    data_mu   = R_double.mean()
    data_sig  = R_double.std()

  
    # perform fit (provided fit function-> gaus, data-> xhist,yhist, and fit param p0)
    popt,pcov = curve_fit(gaus,xhist,yhist ,p0=[data_peak, data_mu, data_sig], maxfev=5000) 

    # get fit parameters
    mu_fit,  mu_fit_err  = popt[1], np.sqrt(pcov[1,1])
    sig_fit, sig_fit_err = popt[2], np.sqrt(pcov[2,2])

    rel_err = sig_fit / mu_fit
    
    # create evenly-spaced points for plotting fit function
    xmin=R_double.min()
    xmax=R_double.max()
    x_fit = np.linspace(xmin, xmax, 500) 
    
    # plot fit function and write fit parameters as label for legend
    axs.plot(x_fit,gaus(x_fit,*popt), color='r', label=r'$\mu:{0:.3f}\pm{1:.3f}$''\n''$\sigma:{2:.3f}\pm{3:.3f}$''\n''$\sigma$ / $\mu:{4:.3f}$'.format(mu_fit, mu_fit_err, sig_fit, sig_fit_err, rel_err))
    axs.legend(frameon=False, loc='upper right')
    
            
            
    #fig.tight_layout()
   
    plt.title('Double Ratio %s_SRC/MF / %s_SRC/MF (systematics)' % (target1,target2), fontsize=16)
    plt.xlabel('Double Ratio Counts', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
        
    plt.show()


else:
    printHelp()
    sys.exit()    
