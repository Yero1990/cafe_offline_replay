import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import sys




ifname= 'cafe_triplet.csv' 

df = pd.read_csv(ifname, comment='#')

# A_SRC / Ca48_SRC
singleR_A_ca48_src                = np.array(df['singleR_A_ca48_src'])
singleR_A_ca48_src_stat_err       = np.array(df['singleR_A_ca48_src_stat_err'])
singleR_A_ca48_src_norm_syst_err  = np.array(df['singleR_A_ca48_src_norm_syst_err'])
singleR_A_ca48_src_RC_syst_err    = np.array(df['singleR_A_ca48_src_RC_syst_err'])
singleR_A_ca48_src_cut_syst_err   = np.array(df['singleR_A_ca48_src_cut_syst_err'])
singleR_A_ca48_src_syst_err       = np.array(df['singleR_A_ca48_src_syst_err'])
singleR_A_ca48_src_tot_err        = np.array(df['singleR_A_ca48_src_tot_err'])



# A_MF / Ca48_MF
singleR_A_ca48_mf                = np.array(df['singleR_A_ca48_mf'])
singleR_A_ca48_mf_stat_err       = np.array(df['singleR_A_ca48_mf_stat_err'])
singleR_A_ca48_mf_norm_syst_err  = np.array(df['singleR_A_ca48_mf_norm_syst_err'])
singleR_A_ca48_mf_RC_syst_err    = np.array(df['singleR_A_ca48_mf_RC_syst_err'])
singleR_A_ca48_mf_cut_syst_err   = np.array(df['singleR_A_ca48_mf_cut_syst_err'])
singleR_A_ca48_mf_syst_err       = np.array(df['singleR_A_ca48_mf_syst_err'])
singleR_A_ca48_mf_tot_err        = np.array(df['singleR_A_ca48_mf_tot_err'])


# double ratio
doubleR                = np.array(df['doubleR'])
doubleR_stat_err       = np.array(df['doubleR_stat_err'])
doubleR_norm_syst_err  = np.array(df['doubleR_norm_syst_err'])
doubleR_RC_syst_err    = np.array(df['doubleR_RC_syst_err'])
doubleR_cut_syst_err   = np.array(df['doubleR_cut_syst_err'])
doubleR_syst_err       = np.array(df['doubleR_syst_err'])
doubleR_tot_err        = np.array(df['doubleR_tot_err'])







A = df['A'] 
NoZ = np.round(df['NoZ'], 1) 
NmZoA = df['NmZoA']

fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False)


# A / Ca48 vs. A (MF)
ax0.errorbar(A, singleR_A_ca48_mf, yerr=singleR_A_ca48_mf_stat_err, marker='o', markersize=7, mfc='r', mec='r', ecolor='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
ax0.errorbar(A, singleR_A_ca48_mf, yerr=singleR_A_ca48_mf_tot_err, marker='o', markersize=7, mfc='k', mec='k', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
ax0.set_title('CaFe Triplet: A / Ca48')
#ax0.set_xlabel('A', fontsize=15)
#ax0.set_xscale('log')
ax0.tick_params(axis='both', which='major', labelsize=15)


# A / Ca48 vs. A (SRC)
ax1.errorbar(A, singleR_A_ca48_src, yerr=singleR_A_ca48_src_stat_err, marker='o', markersize=7, mfc='r', mec='r', ecolor='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
ax1.errorbar(A, singleR_A_ca48_src, yerr=singleR_A_ca48_src_tot_err, marker='o', markersize=7, mfc='k', mec='k', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
#ax1.set_xlabel('A', fontsize=15)
#ax1.set_xscale('log')
ax1.tick_params(axis='both', which='major', labelsize=15)


# A / Ca48 vs. A (double)
ax2.errorbar(A, doubleR, yerr=doubleR_stat_err, marker='o', markersize=7, mfc='r', mec='r', ecolor='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
ax2.errorbar(A, doubleR, yerr=doubleR_tot_err, marker='o', markersize=7, mfc='k', mec='k', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
ax2.set_xlabel('A', fontsize=15)
#ax1.set_xscale('log')
ax2.tick_params(axis='both', which='major', labelsize=15)

plt.xticks([40,48,54])


plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

#----------

'''
fig2, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False)


# A / Ca48 vs. N/Z (MF)
ax0.errorbar(NoZ, singleR_A_ca48_mf, yerr=singleR_A_ca48_mf_stat_err, marker='o', markersize=7, mfc='r', mec='r', ecolor='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
ax0.errorbar(NoZ, singleR_A_ca48_mf, yerr=singleR_A_ca48_mf_tot_err, marker='o', markersize=7, mfc='k', mec='k', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
ax0.set_title('CaFe Triplet: A / Ca48')
#ax0.set_xlabel('A', fontsize=15)
#ax0.set_xscale('log')
ax0.tick_params(axis='both', which='major', labelsize=15)


# A / Ca48 vs. N/Z (SRC)
ax1.errorbar(NoZ, singleR_A_ca48_src, yerr=singleR_A_ca48_src_stat_err, marker='o', markersize=7, mfc='r', mec='r', ecolor='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
ax1.errorbar(NoZ, singleR_A_ca48_src, yerr=singleR_A_ca48_src_tot_err, marker='o', markersize=7, mfc='k', mec='k', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
#ax1.set_xlabel('A', fontsize=15)
#ax1.set_xscale('log')
ax1.tick_params(axis='both', which='major', labelsize=15)


# A / Ca48 vs. N/Z (double)
ax2.errorbar(NoZ, doubleR, yerr=doubleR_stat_err, marker='o', markersize=7, mfc='r', mec='r', ecolor='r', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='statistical')
ax2.errorbar(NoZ, doubleR, yerr=doubleR_tot_err, marker='o', markersize=7, mfc='k', mec='k', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None', label='total error')
ax2.set_xlabel('N/Z', fontsize=15)
#ax1.set_xscale('log')
ax2.tick_params(axis='both', which='major', labelsize=15)

plt.xticks(NoZ)


plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

'''



plt.show()
