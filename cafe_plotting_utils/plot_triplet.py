
#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# This file when run asks for no inputs from the user but should be in the same folder as the appropriately formated .csv files from the ratiohistos.C code
# this file produces cafe triplet plots.
# This file has axis adjusted for per-nucleus normalization
# Note that a lot of values are hard coded and should ideally be read in from a .xlsv file
#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=

import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import sys




ifname= 'cafe_triplet_to40.csv' 

df = pd.read_csv(ifname, comment='#')

# A_SRC / Ca40_SRC
singleR_A_Ca40_src                = np.array(df['singleR_A_ca40_src'])
singleR_A_Ca40_src_stat_err       = np.array(df['singleR_A_ca40_src_stat_err'])
singleR_A_Ca40_src_norm_syst_err  = np.array(df['singleR_A_ca40_src_norm_syst_err'])
singleR_A_Ca40_src_RC_syst_err    = np.array(df['singleR_A_ca40_src_RC_syst_err'])
singleR_A_Ca40_src_cut_syst_err   = np.array(df['singleR_A_ca40_src_cut_syst_err'])
singleR_A_Ca40_src_syst_err       = np.array(df['singleR_A_ca40_src_syst_err'])
singleR_A_Ca40_src_tot_err        = [0.0, 0.040, 0.063]#np.array(df['singleR_A_ca40_src_tot_err'])
#0, 0.037, 0.038
#[1.0, 1.096, 1.638]

# A_MF / Ca40_MF
singleR_A_Ca40_mf                = np.array(df['singleR_A_ca40_mf'])
singleR_A_Ca40_mf_stat_err       = np.array(df['singleR_A_ca40_mf_stat_err'])
singleR_A_Ca40_mf_norm_syst_err  = np.array(df['singleR_A_ca40_mf_norm_syst_err'])
singleR_A_Ca40_mf_RC_syst_err    = np.array(df['singleR_A_ca40_mf_RC_syst_err'])
singleR_A_Ca40_mf_cut_syst_err   = np.array(df['singleR_A_ca40_mf_cut_syst_err'])
singleR_A_Ca40_mf_syst_err       = np.array(df['singleR_A_ca40_mf_syst_err'])
singleR_A_Ca40_mf_tot_err        = np.array(df['singleR_A_ca40_mf_tot_err'])


# double ratio
doubleR                = np.array(df['doubleR'])
doubleR_stat_err       = np.array(df['doubleR_stat_err'])
doubleR_norm_syst_err  = np.array(df['doubleR_norm_syst_err'])
doubleR_RC_syst_err    = np.array(df['doubleR_RC_syst_err'])
doubleR_cut_syst_err   = np.array(df['doubleR_cut_syst_err'])
doubleR_syst_err       = np.array(df['doubleR_syst_err'])
doubleR_tot_err        = np.array(df['doubleR_tot_err'])


















ifname2 = 'cafe_triplet_to48.csv' 

df2 = pd.read_csv(ifname2, comment='#')

# A_SRC / Ca48_SRC
singleR_A_Ca48_src                = np.array(df2['singleR_A_ca48_src'])
singleR_A_Ca48_src_stat_err       = np.array(df2['singleR_A_ca48_src_stat_err'])
singleR_A_Ca48_src_norm_syst_err  = np.array(df2['singleR_A_ca48_src_norm_syst_err'])
singleR_A_Ca48_src_RC_syst_err    = np.array(df2['singleR_A_ca48_src_RC_syst_err'])
singleR_A_Ca48_src_cut_syst_err   = np.array(df2['singleR_A_ca48_src_cut_syst_err'])
singleR_A_Ca48_src_syst_err       = np.array(df2['singleR_A_ca48_src_syst_err'])
singleR_A_Ca48_src_tot_err        = np.array(df2['singleR_A_ca48_src_tot_err'])



# A_MF / Ca48_MF
singleR_A_Ca48_mf                = np.array(df2['singleR_A_ca48_mf'])
singleR_A_Ca48_mf_stat_err       = np.array(df2['singleR_A_ca48_mf_stat_err'])
singleR_A_Ca48_mf_norm_syst_err  = np.array(df2['singleR_A_ca48_mf_norm_syst_err'])
singleR_A_Ca48_mf_RC_syst_err    = np.array(df2['singleR_A_ca48_mf_RC_syst_err'])
singleR_A_Ca48_mf_cut_syst_err   = np.array(df2['singleR_A_ca48_mf_cut_syst_err'])
singleR_A_Ca48_mf_syst_err       = np.array(df2['singleR_A_ca48_mf_syst_err'])
singleR_A_Ca48_mf_tot_err        = np.array(df2['singleR_A_ca48_mf_tot_err'])




# -----------
# Read Models Per Proton
# -----------

# OSU (D. Furnshtal) / SRG
singleR_A_ca40_mf_osu = [1, 0.97864, 0.99084]#[0.983, 0.962, 0.974]
singleR_A_ca40_src_osu = [1.0, 1.16379, 1.01724]#[1.16, 1.35, 1.18]

singleR_A_ca48_mf_osu = [0.983, 0.962, 0.974]
singleR_A_ca48_src_osu = [0.85926, 1.0, 0.87407]

# Colle / Sn=0 Pair
singleR_A_ca40_src_colle = [1.0, 1.2, 1.127848] #really the double

singleR_A_ca48_src_colle = [0.8333, 1.0, 0.93987] #really the double

# Spatial
singleR_A_ca40_src_Jmodel = [1.0, 1.29074, 1.14738]#[1.16, 1.35, 1.18]

singleR_A_ca48_src_Jmodel = [0.77475, 1.0, 0.888932]







A = df['A'] 
NoZ = np.round(df['NoZ'], 1) 
NmZoA = df['NmZoA']























######################################################################
######################################################################
######################################################################

fig, ax7 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax7.errorbar(A, [1.0, 1.096, 1.638], yerr=[0.0, 0.040, 0.063], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax7.plot([54], [1.322412], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax7.plot([54], [1.491594], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#Spatial
ax7.plot([54], [1.4662024], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#l=0, n=0
ax7.plot([54], [1.43], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#l=0, L=0
ax7.plot([54], [1.4662024], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)
ax7.plot([54], [1.43], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)

ax7.plot([48], [1.16379], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
ax7.plot([48], [1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
ax7.plot([48], [1.2], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
ax7.plot([48], [1.0], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
ax7.plot([48], [1.2], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)
ax7.plot([48], [1.0], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)

#ax7.set_title('SRC Enhancement', fontsize=21)
ax7.set_xlabel('A', fontsize=21)
ax7.set_ylabel('Per Nucleus SRC Ratio to $^{40}$Ca', fontsize=21)#, weight='bold')

ax7.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax7.set_xlim(38, 58)
ax7.set_ylim(0.9, 1.8)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=14, loc='lower right')

fig.savefig('SRC_Ca48Ca40_Nucleus.png')


######################################################################
######################################################################
######################################################################


fig, ax10 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax10.errorbar(A, [1.0*(20/20), 1.096*(20/20), 1.638*(20/26)], yerr=[0.0, 0.040, 0.048], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax10.plot([54], [1.43*(20/26)], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, L=0
ax10.plot([54], [1.322412*(20/26)], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax10.plot([54], [1.491594*(20/26)], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#Spatial
ax10.plot([54], [1.4662024*(20/26)], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, n=0

ax10.plot([48], [1.0*(20/20)], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
ax10.plot([48], [1.16379*(20/20)], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
ax10.plot([48], [1.29074*(20/20)], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
ax10.plot([48], [1.2*(20/20)], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

ax10.set_title('SRC Enhancement', fontsize=21)
ax10.set_xlabel('A', fontsize=21)
ax10.set_ylabel('Per Proton SRC Ratio to $^{40}$Ca', fontsize=21)#, weight='bold')

ax10.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax10.set_xlim(38, 58)
ax10.set_ylim(0.9, 1.4)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2, 1.3])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=14, loc='lower right')

fig.savefig('SRC_Ca48Ca40_Proton.png')


######################################################################
######################################################################
######################################################################


fig, ax11 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax11.errorbar(A, [1.0, 1.096*(20/28), 1.638*(20/28)], yerr=[0.0, 0.029, 0.045], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax11.plot([54], [1.43*(20/28)], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, L=0
ax11.plot([54], [1.322412*(20/28)], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax11.plot([54], [1.491594*(20/28)], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#Spatial
ax11.plot([54], [1.4662024*(20/28)], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0,n=0

ax11.plot([48], [1.0*(20/28)], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
ax11.plot([48], [1.16379*(20/28)], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
ax11.plot([48], [1.29074*(20/28)], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
ax11.plot([48], [1.2*(20/28)], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

ax11.set_title('SRC Enhancement', fontsize=21)
ax11.set_xlabel('A', fontsize=21)
ax11.set_ylabel('Per Neutron SRC Ratio to $^{40}$Ca', fontsize=21)#, weight='bold')

ax11.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax11.set_xlim(38, 58)
ax11.set_ylim(0.65, 1.3)
plt.xticks(fontsize = 15)
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('SRC_Ca48Ca40_Neutron.png')









######################################################################
######################################################################
######################################################################

fig, ax12 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax12.errorbar(A, [1.0*(1/1.096), 1.096*(1/1.096), 1.638*(1/1.096)], yerr=[0.0346, 0.0, 0.063], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax12.plot([54], [1.322412/1.16379], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, L=0
ax12.plot([54], [1.491594/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax12.plot([54], [1.4662024/1.2], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#Spatial
ax12.plot([54], [1.43/1.0], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#l=0, n=0
ax12.plot([54], [1.4662024/1.2], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)
ax12.plot([54], [1.43/1.0], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)

#ax12.plot([48], [1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
#ax12.plot([48], [1.16379/1.16379], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
#ax12.plot([48], [1.29074/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
#ax12.plot([48], [1.2/1.2], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

#Check L=0,l=0 ca40 value
ax12.plot([40], [1.0/1.16379], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
ax12.plot([40], [1.0/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
ax12.plot([40], [1.0/1.20], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
ax12.plot([40], [1.0/1.0], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
ax12.plot([40], [1.0/1.20], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)
ax12.plot([40], [1.0/1.0], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)


#ax12.set_title('SRC Enhancement', fontsize=21)
ax12.set_xlabel('A', fontsize=21)
ax12.set_ylabel('Per Nucleus SRC Ratio to $^{48}$Ca', fontsize=21)#, weight='bold')

ax12.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax12.set_xlim(38, 58)
ax12.set_ylim(0.7, 1.7)
plt.xticks(fontsize = 15)
plt.yticks([0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('SRC_Ca40Ca48_Nucleus.png')


######################################################################
######################################################################
######################################################################

fig, ax123 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax123.errorbar([54], [1.638*(1/1.096)], yerr=[0.063], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax123.plot([54], [1.322412/1.16379], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax123.plot([54], [1.491594/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#Spatial
ax123.plot([54], [1.4662024/1.2], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#l=0, n=0
ax123.plot([54], [1.43/1.0], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#l=0, L=0

#ax12.set_title('SRC Enhancement', fontsize=21)
ax123.set_xlabel('54', fontsize=21)
ax123.set_ylabel('Per Nucleus SRC Ratio $^{54}$Fe to $^{48}$Ca', fontsize=21)#, weight='bold')

ax123.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1])
ax123.set_xlim(52, 56)
ax123.set_ylim(1.0, 1.6)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
plt.yticks(fontsize = 15)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('SRC_Fe54Ca48_Nucleus.png')


######################################################################
######################################################################
######################################################################

fig, ax124 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
#last to ca40
#ax124.errorbar(A, [1.0, 1.097, 1.637], yerr=[0.0, 0.028, 0.043], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')
#last to ca48
ax124.errorbar(A, [1.0, 1.097, 1.493], yerr=[0.0, 0.028, 0.039], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

#last to ca40
ax124.plot([54], [1.138], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
#1.138, 1.323
ax124.plot([54], [1.196], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#LCA
#1.196, 1.342
ax124.plot([54], [1.222], marker='^', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', zorder=3)#l=0, n=0
#1.222, 1.466
ax124.plot([54], [1.3], marker='v', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', zorder=3)#L=0, l=0
#1.3

#last to ca48
#.plot([54], [1.138], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
#ax124.plot([54], [1.196], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#LCA
#ax124.plot([54], [1.22], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#l=0, n=0
#ax124.plot([54], [1.22], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, n=0
#ax124.plot([54], [1.3], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', zorder=3)#L=0, l=0
#ax124.plot([54], [1.3], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#L=0, l=0

ax124.plot([48], [1.163], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum-SRG', zorder=3)
ax124.plot([48], [1.122], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum-LCA', zorder=3)
#ax124.plot([48], [1.120], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
ax124.plot([48], [1.200], marker='^', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
#ax124.plot([48], [1.0], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
ax124.plot([48], [1.0], marker='v', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)

#ax12.set_title('SRC Enhancement', fontsize=21)
ax124.set_xlabel('', fontsize=21)
ax124.set_ylabel('Cross Section Ratios', fontsize=23)#, weight='bold')

ax124.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40, 48, 54])
ax124.set_xlim(38, 58)
#ax124.set_ylim(0.9, 1.7)
ax124.set_ylim(0.9, 1.6)
plt.xticks(fontsize = 15)
#plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
plt.yticks(fontsize = 17)

plt.legend(frameon=True, fontsize=14, loc='lower right')

fig.savefig('SRC_Ca40Ca48Fe54_Nucleus.png')



#filled_marker_style = dict(marker='o', linestyle=':', markersize=15,
 #                          color='darkgrey',
 #                          markerfacecolor='tab:blue',
 #                          markerfacecoloralt='lightsteelblue',
 #                          markeredgecolor='brown')




















fig, ax13 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax13.errorbar(A, [1.0*(20/20)*(1/1.096), 1.096*(20/20)*(1/1.096), 1.638*(20/26)*(1/1.096)], yerr=[0.037, 0.0, 0.048], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax13.plot([54], [1.53*(20/26)/1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0,L=0
ax13.plot([54], [1.322412*(20/26)/1.16379], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax13.plot([54], [1.491594*(20/26)/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#Spatial
ax13.plot([54], [1.4662024*(20/26)/1.4662024], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, n=0

#ax13.plot([48], [1.0*(20/20)/1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
#ax13.plot([48], [1.16379*(20/20)/1.16379], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
#ax13.plot([48], [1.29074*(20/20)/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
#ax13.plot([48], [1.2*(20/20)/1.4662024], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

#Check L=0,l=0 ca40 value
ax13.plot([40], [1.0*(20/20)/1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
ax13.plot([40], [1.0*(20/20)/1.16379], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
ax13.plot([40], [1.0*(20/20)/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
ax13.plot([40], [1.0*(20/20)/1.20], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

ax13.set_title('SRC Enhancement', fontsize=21)
ax13.set_xlabel('A', fontsize=21)
ax13.set_ylabel('Per Proton SRC Ratio to $^{48}$Ca', fontsize=21)#, weight='bold')

ax13.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax13.set_xlim(38, 58)
ax13.set_ylim(0.6, 1.35)
plt.xticks(fontsize = 15)
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('SRC_Ca40Ca48_Proton.png')


######################################################################
######################################################################
######################################################################


fig, ax14 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax14.errorbar(A, [1.0*(1/1.096)*(26/20), 1.096*(1/1.096)*(26/26), 1.638*(1/1.096)*(26/28)], yerr=[0.047, 0.0, 0.057], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax14.plot([54], [1.53*(26/28)/1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, L=0
ax14.plot([54], [1.322412*(26/28)/1.16379], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax14.plot([54], [1.491594*(26/28)/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#Spatial
ax14.plot([54], [1.4662024*(26/28)/1.2], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, n=0

#ax14.plot([48], [1.0*(20/26)/1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
#ax14.plot([48], [1.16379*(20/26)/1.16379], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
#ax14.plot([48], [1.29074*(20/26)/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
#ax14.plot([48], [1.2*(20/26)/1.2], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

#Check L=0,l=0 ca40 value
ax14.plot([40], [1.0*(26/20)/1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0, pairs', zorder=3)
ax14.plot([40], [1.0*(26/20)/1.16379], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
ax14.plot([40], [1.0*(26/20)/1.29074], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
ax14.plot([40], [1.0*(26/20)/1.20], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

ax14.set_title('SRC Enhancement', fontsize=21)
ax14.set_xlabel('A', fontsize=21)
ax14.set_ylabel('Per Neutron SRC Ratio to $^{48}$Ca', fontsize=21)#, weight='bold')

ax14.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax14.set_xlim(38, 58)
ax14.set_ylim(0.9, 1.5)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('SRC_Ca40Ca48_Neutron.png')



































































######################################################################
######################################################################
######################################################################

fig, ax9 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax9.errorbar([20, 28, 28], [1.0, 1.096, 1.636], yerr=[0.0, 0.04, 0.063], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax9.plot([28], [1.43], marker='*', markersize=18, alpha=1.0, mfc='none', mec='b', linestyle='None', zorder=3)#l=0, L=0
ax9.plot([28], [1.322412], marker='P', markersize=18, alpha=1.0, mfc='none', mec='b', linestyle='None', label='None', zorder=3)#SRG
ax9.plot([28], [1.491594], marker='o', markersize=18, alpha=1.0, mfc='none', mec='b', linestyle='None', label='None', zorder=3)#spatial
ax9.plot([28], [1.4662024], marker='^', markersize=18, alpha=1.0, mfc='none', mec='b', linestyle='None', label='None', zorder=3)#l=0, n=0

ax9.plot([28], [1.0], marker='*', markersize=18, alpha=1.0, mfc='none', mec='red', linestyle='None', label='l=0, L=0, pairs', zorder=3)
ax9.plot([28], [1.16379], marker='P', markersize=18, alpha=1.0, mfc='none', mec='red', linestyle='None', label='SRG', zorder=3)#osu
ax9.plot([28], [1.29074], marker='o', markersize=18, alpha=1.0, mfc='none', mec='red', linestyle='None', label='Spatial', zorder=3)#osu
ax9.plot([28], [1.2], marker='^', markersize=18, alpha=1.0, mfc='none', mec='red', linestyle='None', label='l=0, n=0 pairs', zorder=3)#osu




ax9.set_title('SRC Enhancement', fontsize=21)
ax9.set_xlabel('N', fontsize=21)
ax9.set_ylabel('Per Nucleus SRC Ratio to $^{40}$Ca', fontsize=21)#, weight='bold')

ax9.tick_params(axis='both', which='major', labelsize=15)
ax9.set_xlim(16, 32)
ax9.set_ylim(0.9, 1.8)
plt.xticks([20, 28])
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

#fig.savefig('SRC_Ca48Ca40_N_c.png')









######################################################################
######################################################################
######################################################################

fig, ax2 = plt.subplots(figsize=(8,7))

# A / Ca40 vs. A (SRC)
ax2.errorbar(A, [1.0, 1.139, 1.188], yerr=[0, 0.021, 0.026], marker='s', markersize=10, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax2.plot([48], [1.188], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG
ax2.plot([48], [1.291], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#Spatial
ax2.plot([48], [1.20], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#l=0, n=0

ax2.plot([54], [1.026], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)
ax2.plot([54], [1.148], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
ax2.plot([54], [1.128], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)


ax2.set_title('Double CaFe Triplet: A / Ca40')
ax2.set_xlabel('A', fontsize=15)
ax2.set_ylabel(r'(SRC/MF)$_A$ / (SRC/MF)$_{\rm Ca40}$', fontsize=12)#, weight='bold')

ax2.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('Double_Ca48Ca40.png')





















































######################################################################
######################################################################
######################################################################

fig, ax20 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax20.errorbar(A, [1.0, 0.963, 1.378], yerr=[0.0, 0.032, 0.046], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax20.plot([48], [0.979], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG

ax20.plot([54], [1.288], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)

#ax20.set_title('SRC Enhancement', fontsize=21)
ax20.set_xlabel('A', fontsize=21)
ax20.set_ylabel('Per Nucleus MF Ratio to $^{40}$Ca', fontsize=21)#, weight='bold')

ax20.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax20.set_xlim(38, 58)
ax20.set_ylim(0.9, 1.5)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('MF_Ca48Ca40_Nucleus.png')


######################################################################
######################################################################
######################################################################


fig, ax21 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax21.errorbar(A, [1.0*(20/20), 0.963*(20/20), 1.378*(20/26)], yerr=[0.0, 0.032, 0.035], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

ax21.plot([48], [0.979], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG

ax21.plot([54], [0.991], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)

#ax21.set_title('SRC Enhancement', fontsize=21)
ax21.set_xlabel('A', fontsize=21)
ax21.set_ylabel('Per Proton MF Ratio to $^{40}$Ca', fontsize=21)#, weight='bold')

ax21.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax21.set_xlim(38, 58)
ax21.set_ylim(0.85, 1.15)
plt.xticks(fontsize = 15)
plt.yticks([0.9, 1.0, 1.1])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('MF_Ca48Ca40_Proton.png')


######################################################################
######################################################################
######################################################################


fig, ax22 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax22.errorbar(A, [1.0, 0.963*(20/28), 1.378*(20/28)], yerr=[0.0, 0.023, 0.033], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')



ax22.plot([54], [1.288*(20/28)], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', zorder=3)#SRG

ax22.plot([48], [0.979*(20/28)], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='SRG', zorder=3)



#ax22.set_title('SRC Enhancement', fontsize=21)
ax22.set_xlabel('A', fontsize=21)
ax22.set_ylabel('Per Neutron MF Ratio to $^{40}$Ca', fontsize=21)#, weight='bold')

ax22.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax22.set_xlim(38, 58)
ax22.set_ylim(0.6, 1.2)
plt.xticks(fontsize = 15)
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1])
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('MF_Ca48Ca40_Neutron.png')









######################################################################
######################################################################
######################################################################

fig, ax23 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax23.errorbar(A, [1.0*(1/0.963), 0.963*(1/0.963), 1.378*(1/0.963)], yerr=[0.033, 0.0, 0.046], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

#ax23.set_title('SRC Enhancement', fontsize=21)
ax23.set_xlabel('A', fontsize=21)
ax23.set_ylabel('Per Nucleus MF Ratio to $^{48}$Ca', fontsize=21)#, weight='bold')

ax23.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax23.set_xlim(38, 58)
ax23.set_ylim(0.9, 1.55)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4])
plt.yticks(fontsize = 15)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('MF_Ca40Ca48_Nucleus.png')


######################################################################
######################################################################
######################################################################


fig, ax24 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax24.errorbar(A, [1.0*(20/20)*(1/0.963), 0.963*(20/20)*(1/0.963), 1.378*(20/26)*(1/0.963)], yerr=[0.033, 0.0, 0.035], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

#ax24.set_title('SRC Enhancement', fontsize=21)
ax24.set_xlabel('A', fontsize=21)
ax24.set_ylabel('Per Proton MF Ratio to $^{48}$Ca', fontsize=21)#, weight='bold')

ax24.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax24.set_xlim(38, 58)
ax24.set_ylim(0.9, 1.3)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2])
plt.yticks(fontsize = 15)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('MF_Ca40Ca48_Proton.png')


######################################################################
######################################################################
######################################################################


fig, ax25 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax25.errorbar(A, [1.0*(1/0.963)*(26/20), 0.963*(1/0.963)*(26/26), 1.378*(1/0.963)*(26/28)], yerr=[0.048, 0.0, 0.048], marker='s', markersize=15, mfc='k', mec='None', ecolor='k', elinewidth=1.2, capsize=4, markeredgewidth=1.2, linestyle='None')

#ax25.set_title('SRC Enhancement', fontsize=21)
ax25.set_xlabel('A', fontsize=21)
ax25.set_ylabel('Per Neutron MF Ratio to $^{48}$Ca', fontsize=21)#, weight='bold')

ax25.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([40,48,54])
ax25.set_xlim(38, 58)
ax25.set_ylim(0.9, 1.5)
plt.xticks(fontsize = 15)
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4])
plt.yticks(fontsize = 15)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('MF_Ca40Ca48_Neutron.png')
