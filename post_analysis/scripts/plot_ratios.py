import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy


if(len(sys.argv)!=2):
    print('''\n
    -----------------------------------------
    Usage: ipython plot_ratios.py <passN>
    -----------------------------------------
    ''')
    sys.exit()


npass = sys.argv[1]       # "pass3"


ifname= 'cafe_ratios_%s_testing.csv' % (npass)

# ----------------
# read input file
# ----------------
df = pd.read_csv(ifname, comment='#')

targ = df['target']


# A_SRC / A_MF
singleR_per_proton                = np.array(df['singleR_per_proton'])
singleR_per_proton_stat_err       = np.array(df['singleR_per_proton_stat_err'])
singleR_per_proton_norm_syst_err  = np.array(df['singleR_per_proton_norm_syst_err'])
singleR_per_proton_RC_syst_err    = np.array(df['singleR_per_proton_RC_syst_err'])
singleR_per_proton_syst_err       = np.array(df['singleR_per_proton_syst_err'])
singleR_per_proton_tot_err        = np.array(df['singleR_per_proton_tot_err'])

# A_MF / C12_MF
singleR_A_c12_mf                 = np.array(df['singleR_A_c12_mf'])
singleR_A_c12_mf_stat_err        = np.array(df['singleR_A_c12_mf_stat_err'])
singleR_A_c12_mf_norm_syst_err   = np.array(df['singleR_A_c12_mf_norm_syst_err'])
singleR_A_c12_mf_RC_syst_err     = np.array(df['singleR_A_c12_mf_RC_syst_err'])
singleR_A_c12_mf_cut_syst_err    = np.array(df['singleR_A_c12_mf_cut_syst_err'])
singleR_A_c12_mf_syst_err        = np.array(df['singleR_A_c12_mf_syst_err'])
singleR_A_c12_mf_tot_err         = np.array(df['singleR_A_c12_mf_tot_err'])

# A_SRC / C12_SRC
singleR_A_c12_src                = np.array(df['singleR_A_c12_src'])
singleR_A_c12_src_stat_err       = np.array(df['singleR_A_c12_src_stat_err'])
singleR_A_c12_src_norm_syst_err  = np.array(df['singleR_A_c12_src_norm_syst_err'])
singleR_A_c12_src_RC_syst_err    = np.array(df['singleR_A_c12_src_RC_syst_err'])
singleR_A_c12_src_cut_syst_err   = np.array(df['singleR_A_c12_src_cut_syst_err'])
singleR_A_c12_src_syst_err       = np.array(df['singleR_A_c12_src_syst_err'])
singleR_A_c12_src_tot_err        = np.array(df['singleR_A_c12_src_tot_err'])

# A_SRC/MF / C12_SRC/MF
doubleR               = np.array(df['doubleR'])
doubleR_stat_err      = np.array(df['doubleR_stat_err'])
doubleR_norm_syst_err = np.array(df['doubleR_norm_syst_err'])
doubleR_RC_syst_err   = np.array(df['doubleR_RC_syst_err'])
doubleR_cut_syst_err  = np.array(df['doubleR_cut_syst_err'])
doubleR_syst_err      = np.array(df['doubleR_syst_err'])
doubleR_tot_err       = np.array(df['doubleR_tot_err'])

# -----------
# Read Models
# -----------

# Justin Estees / Andrew Denniston Model
doubleR_Jmodel               = np.array(df['doubleR_Jmodel'], dtype=float)

singleR_A_c12_mf_av18        = np.array(df['singleR_A_c12_mf_av18'],  dtype=float)
singleR_A_c12_src_av18       = np.array(df['singleR_A_c12_src_av18'], dtype=float)
doubleR_av18                 = np.array(df['doubleR_av18'],           dtype=float)

singleR_A_c12_mf_osu         = np.array(df['singleR_A_c12_mf_osu'],   dtype=float)
singleR_A_c12_src_osu        = np.array(df['singleR_A_c12_src_osu'],  dtype=float)
doubleR_osu                  = np.array(df['doubleR_osu'],            dtype=float)

# mask NaN values for plotting
doubleR_Jmodel               = np.ma.masked_invalid(doubleR_Jmodel)

singleR_A_c12_mf_av18        = np.ma.masked_invalid(singleR_A_c12_mf_av18)
singleR_A_c12_src_av18       = np.ma.masked_invalid(singleR_A_c12_src_av18)
doubleR_av18                 = np.ma.masked_invalid(doubleR_av18)

singleR_A_c12_mf_osu        = np.ma.masked_invalid(singleR_A_c12_mf_osu)
singleR_A_c12_src_osu       = np.ma.masked_invalid(singleR_A_c12_src_osu)
doubleR_osu                 = np.ma.masked_invalid(doubleR_osu)



A = df['A'] 
NoZ = df['NoZ'] 
NmZoA = df['NmZoA'] 
y_rel = np.zeros([len(doubleR_Jmodel)])  # y-axis (at zero) for plotting relative errors

compare_flag = False  # compare to previous pass ?
compare_simc_flag = False  #compare to simc ?
error_breakdown = False     # flag to display breakdown of absolute uncertainties in the ratios
rel_err_breakdown = False   # flag to display breakdown of relative uncertainties in the ratios

if(compare_simc_flag):

    # read SIMC single_ratio (C12, Fe56, Au197)
    df_simc = pd.read_csv('cafe_simc_summary.csv', comment='#')
    A_simc = df_simc['A']
    singleR_per_proton_simc      = np.array(df_simc['singleR_A_to_C12_MF'])
    singleR_per_proton_err_simc  = np.array(df_simc['singleR_A_to_C12_MF_err'])




if(compare_flag):
    fname= 'cafe_ratios_pass2.csv' 

    # read input file
    df = pd.read_csv(fname, comment='#')
    
    targ_2 = df['target']
    singleR_per_proton_2      = np.array(df['singleR_per_proton'])
    singleR_per_proton_err_2  = np.array(df['singleR_per_proton_err'])

    singleR_A_c12_mf_2        = np.array(df['singleR_A_c12_mf'])
    singleR_A_c12_mf_err_2    = np.array(df['singleR_A_c12_mf_err'])
    
    singleR_A_c12_src_2        = np.array(df['singleR_A_c12_src'])
    singleR_A_c12_src_err_2    = np.array(df['singleR_A_c12_src_err'])

    doubleR_2      = np.array(df['doubleR'])
    doubleR_err_2  = np.array(df['doubleR_err'])
    A_2 = df['A'] 
    NoZ_2 = df['NoZ'] 
    NmZoA_2 = df['NmZoA'] 


#--------------------------
# PLOT Double Ratio vs. A
#--------------------------
print(' 100*doubleR_stat_err/doubleR ======> ',  100*doubleR_stat_err/doubleR)
fig1= plt.figure()

if(error_breakdown):
    # break-down of syst. contributions
    plt.errorbar(A, doubleR, doubleR_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(A, doubleR, doubleR_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(A, doubleR, doubleR_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(A, doubleR, doubleR_cut_syst_err,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='cut sensitivity (syst)')
    plt.errorbar(A, doubleR, doubleR_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')

elif(rel_err_breakdown):
    plt.plot(A, 100*doubleR_stat_err/doubleR,      marker='_', markersize=10, mfc='r', mec='r',  markeredgewidth=2, linestyle='None', label='statistical')
    plt.plot(A, 100*doubleR_norm_syst_err/doubleR, marker='_', markersize=10, mfc='b', mec='b',  markeredgewidth=2, linestyle='None', label='normalization (syst)')
    plt.plot(A, 100*doubleR_RC_syst_err/doubleR,   marker='_', markersize=10, mfc='g', mec='g',  markeredgewidth=2, linestyle='None', label='radiative corr (syst)')
    plt.plot(A, 100*doubleR_cut_syst_err/doubleR,  marker='_', markersize=10, mfc='m', mec='m',  markeredgewidth=2, linestyle='None', label='cut sensitivity (syst)')
    plt.plot(A, 100*doubleR_tot_err/doubleR,       marker='_', markersize=10, mfc='k', mec='k',  markeredgewidth=2, linestyle='None', label='total error')
   

else:
    plt.errorbar(A, doubleR, doubleR_stat_err, marker='o', markersize=7, alpha=.9, mfc='r', ecolor='r', mec='r', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical', zorder=1)
    plt.errorbar(A, doubleR, doubleR_tot_err, marker='o', markersize=7, alpha=.9, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', label='total', zorder=2)

    # PLOT MODELS
    plt.plot(A, doubleR_Jmodel, marker='s', markersize=7, alpha=.5, mfc='darkorange', mec='darkorange', linestyle='None', label='Justin/Andrew Model', zorder=3)
    #plt.plot(A, doubleR_av18,   marker='*', markersize=15, alpha=.5, mfc='g', mec='g', linestyle='None', label='AV18', zorder=3)
    #plt.plot(A, doubleR_osu,    marker='P', markersize=15, alpha=.5, mfc='b', mec='b', linestyle='None', label='OSU', zorder=3)


    
if(compare_flag ):
    plt.errorbar(A_2, doubleR_2, doubleR_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

if(rel_err_breakdown):
    plt.title('CaFe Double Ratio Relative Error vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel(r'Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe Double Ratio vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel('A / C12', fontsize=16)

    # add target names to plot
    for i, tgt in enumerate(targ):
        print('i, tgt -> ',i, tgt)
        print('(x,y) ->  ',A[i], doubleR[i] )

        # standard (x,y) coordinates
        x = A[i] + 2
        y = doubleR[i]

        if tgt=="Fe54":
            x = A[i] 
            y = doubleR[i] + 0.07
        elif tgt=="Ca48":
            x = A[i] - 13
        elif tgt=="Ca40":
            x = A[i] - 12
        elif tgt=="Au197":
            x = A[i] - 70
        elif tgt=="B11":
            x = A[i] - 2
            y = doubleR[i] + 0.05
        elif tgt=="Be9" or tgt=="B10" or tgt=="C12":
            x = A[i] + 1
        else:
            x = A[i] + 2

        plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_doubleR_vs_A.pdf')


'''
#--------------------------
# PLOT Single Ratio vs. A
#--------------------------

#--------------------------
# Single Ratio A_SRC / A_MF 
#--------------------------
fig1_s= plt.figure()

if(error_breakdown):
    plt.errorbar(A, singleR_per_proton, singleR_per_proton_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(A, singleR_per_proton, singleR_per_proton_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(A, singleR_per_proton, singleR_per_proton_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(A, singleR_per_proton, singleR_per_proton_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')

elif(rel_err_breakdown):
    plt.errorbar(A, y_rel, 100*singleR_per_proton_stat_err/singleR_per_proton,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(A, y_rel, 100*singleR_per_proton_norm_syst_err/singleR_per_proton, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(A, y_rel, 100*singleR_per_proton_RC_syst_err/singleR_per_proton,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(A, y_rel, 100*singleR_per_proton_tot_err/singleR_per_proton,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')

else:
    plt.errorbar(A, singleR_per_proton, singleR_per_proton_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(A, singleR_per_proton, singleR_per_proton_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')



if(compare_flag ):
    plt.errorbar(A_2, singleR_per_proton_2, singleR_per_proton_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')

if(rel_err_breakdown):
    plt.title('CaFe SRC/MF Single Ratio (per proton) Rel. Error vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel(r'Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe SRC/MF Single Ratio (per proton) vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel('SRC / MF', fontsize=16)

    # add target names to plot
    for i, tgt in enumerate(targ):
        print('i, tgt -> ',i, tgt)
        print('(x,y) ->  ',A[i], singleR_per_proton[i] )

        # standard (x,y) coordinates
        x = A[i] + 2
        y = singleR_per_proton[i]

        if tgt=="Fe54":
            x = A[i] + 15
        elif tgt=="Ca48":
            x = A[i] + 13
        elif tgt=="Ca40":
            x = A[i] + 12
        elif tgt=="Au197":
            x = A[i] - 70
        elif tgt=="B11":
            x = A[i] 
            y = singleR_per_proton[i] +0.0005
        elif tgt=="Be9" or tgt=="B10" or tgt=="C12":
            x = A[i] + 1
        else:
            x = A[i] + 2

        plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_singleR_vs_A.pdf')

# --------------------------------------------
'''
#-----------------------------
# Single Ratio A_MF / C12_MF 
#-----------------------------
fig1_a2c_mf= plt.figure()

if(error_breakdown):
    # break-down of syst. contributions
    plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_cut_syst_err,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='cut sensitivity (syst)')
    plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
    
elif(rel_err_breakdown):
    plt.plot(A, 100*singleR_A_c12_mf_stat_err/singleR_A_c12_mf,      marker='_', markersize=10, mfc='r', mec='r', markeredgewidth=2, linestyle='None', label='statistical')
    plt.plot(A, 100*singleR_A_c12_mf_norm_syst_err/singleR_A_c12_mf, marker='_', markersize=10, mfc='b', mec='b', markeredgewidth=2, linestyle='None', label='normalization (syst)')
    plt.plot(A, 100*singleR_A_c12_mf_RC_syst_err/singleR_A_c12_mf,   marker='_', markersize=10, mfc='g', mec='g', markeredgewidth=2, linestyle='None', label='radiative corr (syst)')
    plt.plot(A, 100*singleR_A_c12_mf_cut_syst_err/singleR_A_c12_mf,  marker='_', markersize=10, mfc='m', mec='m', markeredgewidth=2, linestyle='None', label='cut sensitivity (syst)')
    plt.plot(A, 100*singleR_A_c12_mf_tot_err/singleR_A_c12_mf,       marker='_', markersize=10, mfc='k', mec='k', markeredgewidth=2, linestyle='None', label='total error')
    
else:
    plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_stat_err,      marker='o', markersize=7, alpha=.9, mfc='r', ecolor='r', mec='r', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical', zorder=1)
    plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_tot_err,       marker='o', markersize=7, alpha=.9, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error', zorder=1)

    # PLOT MODELS
    plt.plot(A, singleR_A_c12_mf_av18,   marker='*', markersize=15, alpha=.5, mfc='g', mec='g', linestyle='None', label='AV18', zorder=2)
    plt.plot(A, singleR_A_c12_mf_osu,    marker='P', markersize=15, alpha=.5, mfc='b', mec='b', linestyle='None', label='OSU', zorder=2)


if(compare_flag ):
    plt.errorbar(A_2, singleR_A_c12_mf_2, singleR_A_c12_mf_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
if(compare_simc_flag ):
    plt.errorbar(A_simc, singleR_per_proton_simc, singleR_per_proton_err_simc, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='SIMC')
        
plt.xscale('log')

plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

if(rel_err_breakdown):
    plt.title('CaFe MF Single Ratio (per proton) Rel. Error vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel(r'A / C12  Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe MF Single Ratio (per proton) vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel('A / C12', fontsize=16) 

    # add target names to plot
    for i, tgt in enumerate(targ):
        print('i, tgt -> ',i, tgt)
        print('(x,y) ->  ',A[i], singleR_A_c12_mf[i] )
        
        # standard (x,y) coordinates
        x = A[i] + 2
        y = singleR_A_c12_mf[i] 
   
        if tgt=="Au197":
            x = A[i] - 60
        elif tgt=="Be9" or tgt=="B10" or tgt=="B11" :
            x = A[i] + 0.8
            y = y + 0.05
        elif tgt=="C12":
            x = A[i] + 0.8
            y = y - 0.02
        elif tgt=="Ca40" or tgt=="Ca48" or tgt=="Fe54" :
            x = A[i] - 2
            y = y - 0.06
        plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_MFsingleR_vs_A.pdf')
# --------------------------------------------


#-----------------------------
# Single Ratio A_SRC / C12_SRC 
#-----------------------------
fig1_a2c= plt.figure()

if(error_breakdown):
    # break-down of syst. contributions
    plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='statistical')
    plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='normalization (syst)')
    plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='radiative corr (syst)')
    plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_cut_syst_err,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='cut sensitivity (syst)')
    plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='total error')

elif(rel_err_breakdown):    
    plt.plot(A, 100*singleR_A_c12_src_stat_err/singleR_A_c12_src,      marker='_', markersize=10, mfc='r', mec='r', markeredgewidth=2, linestyle='None',  label='statistical')
    plt.plot(A, 100*singleR_A_c12_src_norm_syst_err/singleR_A_c12_src, marker='_', markersize=10, mfc='b', mec='b', markeredgewidth=2, linestyle='None',  label='normalization (syst)')
    plt.plot(A, 100*singleR_A_c12_src_RC_syst_err/singleR_A_c12_src,   marker='_', markersize=10, mfc='g', mec='g', markeredgewidth=2, linestyle='None',  label='radiative corr (syst)')
    plt.plot(A, 100*singleR_A_c12_src_cut_syst_err/singleR_A_c12_src,  marker='_', markersize=10, mfc='m', mec='m', markeredgewidth=2, linestyle='None',  label='cut sensitivity (syst)')
    plt.plot(A, 100*singleR_A_c12_src_tot_err/singleR_A_c12_src,       marker='_', markersize=10, mfc='k', mec='k', markeredgewidth=2, linestyle='None',  label='total error')

else:
    plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='statistical', zorder=1)
    plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='total error', zorder=1)

    # PLOT MODELS
    plt.plot(A, singleR_A_c12_src_av18,   marker='*', markersize=15, alpha=.5, mfc='g', mec='g', linestyle='None', label='AV18', zorder=2)
    plt.plot(A, singleR_A_c12_src_osu,    marker='P', markersize=15, alpha=.5, mfc='b', mec='b', linestyle='None', label='OSU', zorder=2)

    
if(compare_flag ):
    plt.errorbar(A_2, singleR_A_c12_src_2, singleR_A_c12_src_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')

plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

if(rel_err_breakdown):
    plt.title('CaFe SRC Single Ratio (per proton) Rel. Error vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel(r'A / C12 Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe SRC Single Ratio (per proton) vs. A', fontsize=18)
    plt.xlabel('A', fontsize=16)
    plt.ylabel('A / C12', fontsize=16)

    # add target names to plot
    for i, tgt in enumerate(targ):
        print('i, tgt -> ',i, tgt)
        print('(x,y) ->  ',A[i], singleR_A_c12_src[i] )
        
        # standard (x,y) coordinates
        x = A[i] + 2
        y = singleR_A_c12_src[i]
        
        if tgt=="Fe54":
            x = A[i] + 6
        elif tgt=="Ca48":
            x = A[i]
            y = singleR_A_c12_src[i] - 0.09
        elif tgt=="Ca40":
            x = A[i] - 2
            y = singleR_A_c12_src[i] - 0.09
        elif tgt=="Au197":
            x = A[i] - 60
        elif tgt=="Be9":
            x = A[i] - 1
            y = singleR_A_c12_src[i] +0.09
        elif tgt=="B10":
            x = A[i] - 1
            y = singleR_A_c12_src[i] +0.08
        elif tgt=="B11":
            x = A[i] - 1
            y = singleR_A_c12_src[i] +0.08         
        elif tgt=="C12":
            x = A[i] + 1

        else:
            x = A[i] + 2

        plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_SRCsingleR_vs_A.pdf')
# --------------------------------------------

'''
#--------------------------
# PLOT Single Ratio vs. N/Z
#--------------------------

# A_SRC / A_MF
fig2_s = plt.figure()

if(error_breakdown):
    plt.errorbar(NoZ, singleR_per_proton, singleR_per_proton_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(NoZ, singleR_per_proton, singleR_per_proton_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(NoZ, singleR_per_proton, singleR_per_proton_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(NoZ, singleR_per_proton, singleR_per_proton_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')

elif(rel_err_breakdown):
    plt.errorbar(NoZ, y_rel, 100*singleR_per_proton_stat_err/singleR_per_proton,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(NoZ, y_rel, 100*singleR_per_proton_norm_syst_err/singleR_per_proton, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(NoZ, y_rel, 100*singleR_per_proton_RC_syst_err/singleR_per_proton,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(NoZ, y_rel, 100*singleR_per_proton_tot_err/singleR_per_proton,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')

else:
    plt.errorbar(NoZ, singleR_per_proton, singleR_per_proton_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(NoZ, singleR_per_proton, singleR_per_proton_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')



if(compare_flag ):
    plt.errorbar(NoZ_2, singleR_per_proton_2, singleR_per_proton_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

if(rel_err_breakdown):
    plt.title('CaFe SRC/MF Single Ratio (per proton) Rel. Error vs. N/Z', fontsize=18)
    plt.xlabel('N/Z', fontsize=16)
    plt.ylabel(r'Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe SRC/MF Single Ratio (per proton) vs. N/Z', fontsize=18)
    plt.xlabel('N/Z', fontsize=16)
    plt.ylabel('SRC / MF', fontsize=16)

    for i, tgt in enumerate(targ):
    
        x = NoZ[i] + 0.01
        y = singleR_per_proton[i]

        if tgt=="Au197":
            x = NoZ[i] - 0.07
 
        plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_singleR_vs_NoZ.pdf')


# A_SRC / C12_SRC
fig2_s_src = plt.figure()
if(error_breakdown):
    # break-down of syst. contributions
    plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='statistical')
    plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='normalization (syst)')
    plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='radiative corr (syst)')
    plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_cut_syst_err,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='cut sensitivity (syst)')
    plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='total error')

elif(rel_err_breakdown):    
    plt.errorbar(NoZ, y_rel, 100*singleR_A_c12_src_stat_err/singleR_A_c12_src,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='statistical')
    plt.errorbar(NoZ, y_rel, 100*singleR_A_c12_src_norm_syst_err/singleR_A_c12_src, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='normalization (syst)')
    plt.errorbar(NoZ, y_rel, 100*singleR_A_c12_src_RC_syst_err/singleR_A_c12_src,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='radiative corr (syst)')
    plt.errorbar(NoZ, y_rel, 100*singleR_A_c12_src_cut_syst_err/singleR_A_c12_src,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='cut sensitivity (syst)')
    plt.errorbar(NoZ, y_rel, 100*singleR_A_c12_src_tot_err/singleR_A_c12_src,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='total error')

else:
    plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='statistical')
    plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None',  label='total error')



if(compare_flag ):
    plt.errorbar(NoZ_2, singleR_A_c12_src_2, singleR_A_c12_src_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

if(rel_err_breakdown):
    plt.title('CaFe SRC Single Ratio (per proton) Rel. Error vs. N/Z', fontsize=18)
    plt.xlabel('N/Z', fontsize=16)
    plt.ylabel(r'A / C12 Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe SRC Single Ratio (per proton) vs. N/Z', fontsize=18)
    plt.xlabel('N/Z', fontsize=16)
    plt.ylabel('A / C12', fontsize=16)

    for i, tgt in enumerate(targ):
    
        x = NoZ[i] + 0.01
        y = singleR_A_c12_src[i]
 
        plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_SRCsingleR_vs_NoZ.pdf')
'''

#--------------------------
# PLOT Double Ratio vs. N/Z
#--------------------------
fig2 = plt.figure()
if(error_breakdown):
    # break-down of syst. contributions
    plt.errorbar(NoZ, doubleR, doubleR_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(NoZ, doubleR, doubleR_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(NoZ, doubleR, doubleR_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(NoZ, doubleR, doubleR_cut_syst_err,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='cut sensitivity (syst)')
    plt.errorbar(NoZ, doubleR, doubleR_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')

elif(rel_err_breakdown):
    plt.plot(NoZ, 100*doubleR_stat_err/doubleR,      marker='_', markersize=10, mfc='r', mec='r', markeredgewidth=2, linestyle='None', label='statistical')
    plt.plot(NoZ, 100*doubleR_norm_syst_err/doubleR, marker='_', markersize=10, mfc='b', mec='b', markeredgewidth=2, linestyle='None', label='normalization (syst)')
    plt.plot(NoZ, 100*doubleR_RC_syst_err/doubleR,   marker='_', markersize=10, mfc='g', mec='g', markeredgewidth=2, linestyle='None', label='radiative corr (syst)')
    plt.plot(NoZ, 100*doubleR_cut_syst_err/doubleR,  marker='_', markersize=10, mfc='m', mec='m', markeredgewidth=2, linestyle='None', label='cut sensitivity (syst)')
    plt.plot(NoZ, 100*doubleR_tot_err/doubleR,       marker='_', markersize=10, mfc='k', mec='k', markeredgewidth=2, linestyle='None', label='total error')
   

else:
    plt.errorbar(NoZ, doubleR, doubleR_stat_err, marker='o', markersize=7, alpha=.9, mfc='r', ecolor='r', mec='r', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical', zorder=1)
    plt.errorbar(NoZ, doubleR, doubleR_tot_err, marker='o', markersize=7, alpha=.9, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', label='total',zorder=2)

    # PLOT MODELS
    plt.plot(NoZ, doubleR_Jmodel, marker='s', markersize=7, alpha=.5, mfc='darkorange', mec='darkorange', linestyle='None', label='Justin/Andrew Model', zorder=3)
    #plt.plot(NoZ, doubleR_av18,   marker='*', markersize=15, alpha=.5, mfc='g', mec='g', linestyle='None', label='AV18', zorder=3)
    #plt.plot(NoZ, doubleR_osu,    marker='P', markersize=15, alpha=.5, mfc='b', mec='b', linestyle='None', label='OSU', zorder=3)


if(compare_flag ):
    plt.errorbar(NoZ_2, doubleR_2, doubleR_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

if(rel_err_breakdown):
    plt.title('CaFe Double Ratio Relative Error vs. N/Z', fontsize=18)
    plt.xlabel('N/Z', fontsize=16)
    plt.ylabel(r'Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe Double Ratio vs. N/Z', fontsize=18)
    plt.xlabel('N/Z', fontsize=16)
    plt.ylabel('A / C12', fontsize=16)
    
    for i, tgt in enumerate(targ):
        print('i, tgt -> ',i, tgt)
        print('(x,y) ->  ',NoZ[i], doubleR[i] )
        
        x = NoZ[i] + 0.01
        y = doubleR[i]
        
        if tgt=="Au197":
            x = NoZ[i] - 0.07
 
        plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_doubleR_vs_NoZ.pdf')

'''
#-------------------------------
# PLOT Double Ratio vs. (N-Z)/A
#-------------------------------
fig3 = plt.figure()

if(error_breakdown):
    print('OK')
    # break-down of syst. contributions
    plt.errorbar(NmZoA, doubleR, doubleR_stat_err,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(NmZoA, doubleR, doubleR_norm_syst_err, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(NmZoA, doubleR, doubleR_RC_syst_err,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(NmZoA, doubleR, doubleR_cut_syst_err,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='cut sensitivity (syst)')
    plt.errorbar(NmZoA, doubleR, doubleR_tot_err,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')

elif(rel_err_breakdown):
    plt.errorbar(NmZoA, y_rel, 100*doubleR_stat_err/doubleR,      marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(NmZoA, y_rel, 100*doubleR_norm_syst_err/doubleR, marker='o', markersize=7, mfc='b', ecolor='b', mec='b', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='normalization (syst)')
    plt.errorbar(NmZoA, y_rel, 100*doubleR_RC_syst_err/doubleR,   marker='o', markersize=7, mfc='g', ecolor='g', mec='g', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='radiative corr (syst)')
    plt.errorbar(NmZoA, y_rel, 100*doubleR_cut_syst_err/doubleR,  marker='o', markersize=7, mfc='m', ecolor='m', mec='m', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='cut sensitivity (syst)')
    plt.errorbar(NmZoA, y_rel, 100*doubleR_tot_err/doubleR,       marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
    plt.ylim(-10,10)

else:
    plt.errorbar(NmZoA, doubleR, doubleR_stat_err, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
    plt.errorbar(NmZoA, doubleR, doubleR_tot_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total')
    plt.plot(NmZoA, doubleR_Jmodel, marker='o', markersize=7, mfc='darkorange', mec='darkorange', linestyle='None', label='Justin/Andrew Model')



if(compare_flag ):
    plt.errorbar(NmZoA_2, doubleR_2, doubleR_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')


if(rel_err_breakdown):
    plt.title('CaFe Double Ratio Relative Error vs. (N-Z)/A', fontsize=18)
    plt.xlabel('(N-Z)/A', fontsize=16)
    plt.ylabel(r'Relative Error ($\%$)', fontsize=16)
else:
    plt.title('CaFe Double Ratio vs. A', fontsize=18)
    plt.xlabel('(N-Z)/A', fontsize=16)
    plt.ylabel('A / C12', fontsize=16)
    
    for i, tgt in enumerate(targ):
        print('i, tgt -> ',i, tgt)
        print('(x,y) ->  ',NmZoA[i], doubleR[i] )
        
        x = NmZoA[i] + 0.005
        y = doubleR[i]
        
        if tgt=="Au197":
            x = NmZoA[i] - 0.03
 
        plt.text(x, y, tgt)


plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_doubleR_vs_NmZoA.pdf')
'''

plt.show()

