import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
import pandas as pd
from matplotlib import rc

rc('text', usetex=True)

#df_scaler = pd.read_csv("ca48_report_files/fall2021_bcm/ca48_bcmFall2021.csv")
df_scaler = pd.read_csv("ca48_scaler_report_files/ca48_bcmInfo.csv")

plt.rcParams["font.family"] = "Times New Roman"


run = df_scaler['run']
kin = df_scaler['kin']

# fall 2021 bcm calib
T2_scaler = df_scaler['T2_scaler']  # T2 (SHMS EL-REAL)
bcm1_charge = df_scaler['q_bcm1']   # charge [mC]
bcm2_charge = df_scaler['q_bcm2']
bcm4a_charge = df_scaler['q_bcm4a']
bcm4b_charge = df_scaler['q_bcm4b']
bcm4c_charge = df_scaler['q_bcm4c']
I_bcm4a_avg  = df_scaler['I_bcm4a_avg']

# summer 2022 bcm calib
T2_scaler_new = df_scaler['T2_scaler_new']  # T2 (SHMS EL-REAL)
bcm1_charge_new = df_scaler['q_bcm1_new']   # charge [mC]
bcm2_charge_new = df_scaler['q_bcm2_new']
bcm4a_charge_new = df_scaler['q_bcm4a_new']
bcm4b_charge_new = df_scaler['q_bcm4b_new']
bcm4c_charge_new = df_scaler['q_bcm4c_new']
I_bcm4a_avg_new  = df_scaler['I_avg_new']


# calculate T2 scaler counts per charge (fall 2021 bcm calib)
T2_per_mC_bcm1 = T2_scaler / bcm1_charge
T2_per_mC_bcm2 = T2_scaler / bcm2_charge
T2_per_mC_bcm4a = T2_scaler / bcm4a_charge
T2_per_mC_bcm4b = T2_scaler / bcm4b_charge
T2_per_mC_bcm4c = T2_scaler / bcm4c_charge

# calculate T2 scaler counts per charge (summer 2022 bcm calib)
T2_per_mC_bcm1_new = T2_scaler_new / bcm1_charge_new
T2_per_mC_bcm2_new = T2_scaler_new / bcm2_charge_new
T2_per_mC_bcm4a_new = T2_scaler_new / bcm4a_charge_new
T2_per_mC_bcm4b_new = T2_scaler_new / bcm4b_charge_new
T2_per_mC_bcm4c_new = T2_scaler_new / bcm4c_charge_new

# normalize T2 by 1st SRC run (fall 2021 bcm calib)
T2_per_mC_bcm4a_norm = T2_per_mC_bcm4a / T2_per_mC_bcm4a[2]

# normalize T2 by 1st SRC run (summer 2022 bcm calib)
T2_per_mC_bcm4a_norm_new = T2_per_mC_bcm4a_new / T2_per_mC_bcm4a_new[2]

'''
# fall 2021 bcm calibrations
fig0, (ax1, ax2) = plt.subplots(2)
ax1.set_title('Charge-normalized T2 Scaler Counts Relative to 1st SRC Run', fontsize=16, fontweight='bold')
ax1.plot(run[:2],   T2_per_mC_bcm4a_norm[:2], marker='o', color='blue', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 MF')
ax1.plot(run[2:23], T2_per_mC_bcm4a_norm[2:23], marker='o', color='red', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 SRC')
ax1.plot(run[23:],  T2_per_mC_bcm4a_norm[23:], marker='o', color='green', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 MF (round 2)')
ax1.set_ylim([0.9, 1.05])
ax1.set_xticklabels([])
ax1.grid()
ax1.legend(loc='upper left',fontsize=12)
ax1.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)

ax2.set_title('Beam Current', fontsize=16, fontweight='bold')
ax2.plot(run[:2],   I_bcm4a_avg[:2], marker='o', color='blue', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 MF')
ax2.plot(run[2:23], I_bcm4a_avg[2:23], marker='o', color='red', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 SRC')
ax2.plot(run[23:],  I_bcm4a_avg[23:], marker='o', color='green', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 MF (round 2)')
ax2.set_ylim([25, 65])
ax2.grid()
ax2.legend(loc='upper left', fontsize=12)
ax2.set_ylabel('Averge Beam Current [uA]', fontsize=16)
ax2.tick_params(axis='both', which='both', labelsize=15)
ax2.set_xlabel('Run Number', fontsize=16)
plt.show()


# summer 2022 bcm calibrations
fig0, (ax1, ax2) = plt.subplots(2)

ax1.set_title('Charge-normalized T2 Scaler Counts Relative to 1st SRC Run', fontsize=16, fontweight='bold')
ax1.plot(run[:2],   T2_per_mC_bcm4a_norm_new[:2], marker='^', color='blue', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 MF (bcm_calib22)')
ax1.plot(run[2:23], T2_per_mC_bcm4a_norm_new[2:23], marker='^', color='red', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 SRC (bcm_calib22)')
ax1.plot(run[23:],  T2_per_mC_bcm4a_norm_new[23:], marker='^', color='green', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 MF (round 2, bcm_calib22)')
ax1.set_ylim([0.9, 1.05])
ax1.set_xticklabels([])
ax1.grid()
ax1.legend(loc='upper left',fontsize=12)
ax1.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)

ax2.set_title('Beam Current', fontsize=16, fontweight='bold')
ax2.plot(run[:2],   I_bcm4a_avg_new[:2], marker='^', color='blue', markerfacecolor='white',   markersize=8, linestyle='None', label='Ca48 MF (bcm_calib22)')
ax2.plot(run[2:23], I_bcm4a_avg_new[2:23], marker='^', color='red', markerfacecolor='white',  markersize=8, linestyle='None', label='Ca48 SRC (bcm_calib22)')
ax2.plot(run[23:],  I_bcm4a_avg_new[23:], marker='^', color='green', markerfacecolor='white', markersize=8, linestyle='None', label='Ca48 MF (round 2, bcm_calib22)')
ax2.set_ylim([25, 65])
ax2.grid()
ax2.legend(loc='upper left', fontsize=12)
ax2.set_ylabel('Averge Beam Current [uA]', fontsize=16)
ax2.tick_params(axis='both', which='both', labelsize=15)
ax2.set_xlabel('Run Number', fontsize=16)
plt.show()
'''

'''
# bcm calib sanity check
fig0, (ax1, ax2) = plt.subplots(2)
ax1.set_title('Fall 2021 BCM Calibration Parameters', fontsize=16, fontweight='bold')
ax1.plot(run,  bcm1_charge, marker='^', color='blue', markerfacecolor='None', markersize=9, linestyle='None', label='BCM1')
ax1.plot(run,  bcm2_charge, marker='s', color='green', markerfacecolor='None', markersize=9, linestyle='None', label='BCM2')
ax1.plot(run,  bcm4a_charge, marker='o', color='red', markerfacecolor='None', markersize=9, linestyle='None', label='BCM4A')
ax1.plot(run,  bcm4b_charge, marker='D', color='magenta', markerfacecolor='None', markersize=9, linestyle='None', label='BCM4B')
ax1.plot(run,  bcm4c_charge, marker='+', color='violet', markerfacecolor='None', markersize=9, linestyle='None', label='BCM4C')
ax1.set_ylim([-50, 230])
ax1.set_xticklabels([])
ax1.grid()
ax1.legend(loc='upper left')
ax1.set_ylabel('BCM Charge [mC]', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)

ax2.set_title('Summer 2022 BCM Calibration Parameters', fontsize=16, fontweight='bold')
ax2.plot(run,  bcm1_charge_new, marker='^', color='blue', markerfacecolor='None', markersize=9, linestyle='None', label='BCM1')
ax2.plot(run,  bcm2_charge_new, marker='s', color='green', markerfacecolor='None', markersize=9, linestyle='None', label='BCM2')
ax2.plot(run,  bcm4a_charge_new, marker='o', color='red', markerfacecolor='None', markersize=9, linestyle='None', label='BCM4A')
ax2.plot(run,  bcm4b_charge_new, marker='D', color='magenta', markerfacecolor='None', markersize=9, linestyle='None', label='BCM4B')
ax2.plot(run,  bcm4c_charge_new, marker='+', color='violet', markerfacecolor='None', markersize=9, linestyle='None', label='BCM4C')
ax2.set_ylim([0, 230])
ax2.grid()
ax2.legend(loc='upper left', fontsize=12)
ax2.set_ylabel('BCM Charge [mC]', fontsize=16)
ax2.tick_params(axis='both', which='both', labelsize=15)
ax2.set_xlabel('Run Number', fontsize=16)
plt.show()
'''



#-----------------------------------
# Select Ca48 SRC runs with Ib>50 uA and do exponential fit
# (parametrize the percentage drop in T2 scalers versus charge)
#-----------------------------------
'''
run_src = list(run[2:8]) + list(run[12:23])
T2_scl_per_mC_src = list(T2_per_mC_bcm4a_norm_new[2:8]) + list(T2_per_mC_bcm4a_norm_new[12:23])
Q_src = list(bcm4a_charge_new[2:8]) +  list(bcm4a_charge_new[12:23])
Q_csum_src = np.cumsum(Q_src)

T2_scl_per_mC_src = np.array(T2_scl_per_mC_src)
Q_csum_src = np.array(Q_csum_src)
print('x=',Q_csum_src)
print('y=',T2_scl_per_mC_src)

# define exponential-fit function with parameters (a,b)
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

# fit the data
popt, pcov = curve_fit(func, Q_csum_src,  T2_scl_per_mC_src, p0=[1, 1e-05, 200])

print('a = ', popt[0])
print('b = ', popt[1])
print('c = ', popt[2])

print(pcov)

# error in parameters is squareroot of diagonal elements of covariant matrix (diagonal elements -> variances (sig^2))
# off diagonal element represents correlations errors
p_sigma = np.sqrt(np.diag(pcov))

print(p_sigma)

# plot data 
fig0, (ax1) = plt.subplots(1)
ax1.set_title('Charge-normalized T2 Scaler Counts Relative to 1st SRC Run', fontsize=16, fontweight='bold')
ax1.plot(Q_csum_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 SRC data ($I_{b} > 50 \mu$A')
ax1.set_ylim([0.96, 1.02])
ax1.grid()
ax1.legend(loc='upper left',fontsize=12)
ax1.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax1.set_xlabel('Cumulative Charge Q [mC]', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)
# plot fit curve 
ax1.plot(Q_csum_src, func(Q_csum_src, *popt), 'r-', label=r"fit: $A e^{-\alpha Q}$ + C" "\n" "A = %.3E $\pm$ %.3E" "\n" r"$\alpha$ = %.3E $\pm$ %.3E" "\n" "C = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1], popt[2], p_sigma[2]) )
ax1.legend(loc='upper left',fontsize=15)

plt.show()
'''


#-----------------------------------
# Select Ca48 SRC runs with Ib<55 uA
# and FIT the percentage drop
# in T2 scalers versus current) -part1 
# (using the 1st 7 runs < 55 uA)
#-----------------------------------

'''
# 1st 5 outlier runs
#run_src = list(run[8:13])
#T2_scl_per_mC_src = list(T2_per_mC_bcm4a_norm_new[8:13])
#I_bcm4a_avg_src = list(I_bcm4a_avg_new[8:13])
#Q_src = list(bcm4a_charge_new[8:13])

#  1st 5 outlier runs + next 2 runs
run_src = list(run[8:15])
T2_scl_per_mC_src = list(T2_per_mC_bcm4a_norm_new[8:15])
I_bcm4a_avg_src = list(I_bcm4a_avg_new[8:15])
Q_src = list(bcm4a_charge_new[8:15])

Q_csum_src = np.cumsum(Q_src)

T2_scl_per_mC_src = np.array(T2_scl_per_mC_src)
I_bcm4a_avg_src = np.array(I_bcm4a_avg_src)
Q_csum_src = np.array(Q_csum_src)

print('Q=',Q_csum_src)
print('T2_scaler=',T2_scl_per_mC_src)
print('Ib = ', I_bcm4a_avg_src)

# define line fit
def func(x, m, b):
    return m * x + b

# fit the data
popt, pcov = curve_fit(func, I_bcm4a_avg_src, T2_scl_per_mC_src, p0=[0.5, 0.1])
p_sigma = np.sqrt(np.diag(pcov))


fig0, (ax1, ax2) = plt.subplots(2)
ax1.set_title('Charge-normalized T2 Scaler Counts Relative to 1st SRC Run', fontsize=16, fontweight='bold')

ax1.plot(I_bcm4a_avg_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 SRC data ($I_{b} < 55 \mu$A)')

ax1.set_ylim([0.93, 1.02])
ax1.grid()
ax1.legend(loc='upper left',fontsize=12)
ax1.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax1.set_xlabel(r'Beam Current [$\mu$A]', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)


# plot fit curve 
ax1.plot(I_bcm4a_avg_src, func(I_bcm4a_avg_src, *popt), 'r-', label=r"fit: m$x$ + b" "\n" "m = %.3E $\pm$ %.3E" "\n" "b = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1]) )

ax1.legend(loc='upper left',fontsize=12)

#-----

ax2.plot(Q_csum_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 SRC data ($I_{b} < 50 \mu$A')

ax2.set_ylim([0.93, 1.02])
ax2.grid()
ax2.legend(loc='upper left',fontsize=12)
ax2.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax2.set_xlabel(r'Beam Charge [mC]', fontsize=16)
ax2.tick_params(axis='both', which='both', labelsize=15)

plt.show()
'''



'''
#-----------------------------------
# Select Ca48 SRC runs with Ib<50 uA
# and FIT the percentage drop
# in T2 scalers versus current) -part1 
# (using the next two runs at I <~50 and 55 uA)
# runs 17048, 17049
#-----------------------------------

run_src = list(run[13:16])
T2_scl_per_mC_src = list(T2_per_mC_bcm4a_norm_new[13:16])
I_bcm4a_avg_src = list(I_bcm4a_avg_new[13:16])
Q_src = list(bcm4a_charge_new[13:16])
Q_csum_src = np.cumsum(Q_src)

T2_scl_per_mC_src = np.array(T2_scl_per_mC_src)
I_bcm4a_avg_src = np.array(I_bcm4a_avg_src)
Q_csum_src = np.array(Q_csum_src)

print('Q=',Q_csum_src)
print('T2_scaler=',T2_scl_per_mC_src)
print('Ib = ', I_bcm4a_avg_src)

# define line fit
def func(x, m, b):
    return m * x + b

# fit the data
popt, pcov = curve_fit(func, I_bcm4a_avg_src, T2_scl_per_mC_src, p0=[0.5, 1])
p_sigma = np.sqrt(np.diag(pcov))

print('pcov=',pcov)


fig0, (ax1, ax2) = plt.subplots(2)
ax1.set_title('Charge-normalized T2 Scaler Counts Relative to 1st SRC Run', fontsize=16, fontweight='bold')

ax1.plot(I_bcm4a_avg_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 SRC data ($I_{b} ~ 50, 55) \mu$A')

ax1.set_ylim([0.93, 1.02])
ax1.grid()
ax1.legend(loc='upper left',fontsize=12)
ax1.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax1.set_xlabel(r'Beam Current [$\mu$A]', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)


# plot fit curve 
ax1.plot(I_bcm4a_avg_src, func(I_bcm4a_avg_src, *popt), 'r-', label=r"fit: m$x$ + b" "\n" "m = %.3E $\pm$ %.3E" "\n" "b = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1]) )

ax1.legend(loc='upper left',fontsize=12)


#-----

ax2.plot(Q_csum_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 SRC data ($I_{b} ~ 50, 55) \mu$A')

ax2.set_ylim([0.93, 1.02])
ax2.grid()
ax2.legend(loc='upper left',fontsize=12)
ax2.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax2.set_xlabel(r'Beam Charge [mC]', fontsize=16)
ax2.tick_params(axis='both', which='both', labelsize=15)

plt.show()
'''


'''
#-----------------------------------
# Select last 3 Ca48 MF runs
# and FIT the percentage drop
# in T2 scalers versus current)
# RUNS 17093, 17904, 17096
#-----------------------------------

run_src = list(run[-3:])
T2_scl_per_mC_src = list(T2_per_mC_bcm4a_norm_new[-3:])
I_bcm4a_avg_src = list(I_bcm4a_avg_new[-3:])
Q_src = list(bcm4a_charge_new[-3:])
Q_csum_src = np.cumsum(Q_src)

T2_scl_per_mC_src = np.array(T2_scl_per_mC_src)
I_bcm4a_avg_src = np.array(I_bcm4a_avg_src)
Q_csum_src = np.array(Q_csum_src)

print('Q=',Q_csum_src)
print('T2_scaler=',T2_scl_per_mC_src)
print('Ib = ', I_bcm4a_avg_src)

# define line fit
def func(x, m, b):
    return m * x + b

# fit the data
popt, pcov = curve_fit(func, I_bcm4a_avg_src, T2_scl_per_mC_src, p0=[0.5, 0.1])
p_sigma = np.sqrt(np.diag(pcov))


fig0, (ax1, ax2) = plt.subplots(2)
ax1.set_title('Charge-normalized T2 Scaler Counts Relative to 1st SRC Run', fontsize=16, fontweight='bold')

ax1.plot(I_bcm4a_avg_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 MF data (last 3)')

ax1.set_ylim([0.93, 1.02])
ax1.grid()
ax1.legend(loc='upper left',fontsize=12)
ax1.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax1.set_xlabel(r'Beam Current [$\mu$A]', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)


# plot fit curve 
ax1.plot(I_bcm4a_avg_src, func(I_bcm4a_avg_src, *popt), 'r-', label=r"fit: m$x$ + b" "\n" "m = %.3E $\pm$ %.3E" "\n" "b = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1]) )

ax1.legend(loc='upper left',fontsize=12)


#-----

ax2.plot(Q_csum_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 MF data (last 3)')

ax2.set_ylim([0.93, 1.02])
ax2.grid()
ax2.legend(loc='upper left',fontsize=12)
ax2.set_ylabel('Relative T2 scalers/mC', fontsize=16)
ax2.set_xlabel(r'Beam Charge [mC]', fontsize=16)
ax2.tick_params(axis='both', which='both', labelsize=15)

plt.show()
'''



# With the known slope fromnthe fits (in relative T2/mC per uA), one can correct the
# runs relative to 55 uA bem current.

# slope was taken to be the average of the three lineaer fits:  m = 1.636e-3, b= 0.897
# now that the current dependecy is parametrized, we know by how much to correct for beam current
# T2_corrected = T2 + (1.636e-3 * 55 - 1.636e-3 * 35),    T2 (@ 55 uA): 1.636e-3 * 55,  T2 (@ Ib): 1.636e-3 * Ib
# the difference between these two, is the correction factor to be added for T2 (@ Ib)
# T2_corrected = T2 (@ Ib) + ( T2 (@ 55 uA) - T2 (@ Ib) )


#-----------------------------------
# Select Ca48 SRC runs with Ib<55 uA
# and FIT the percentage drop
# in T2 scalers versus current) -part1 
# (using the next two runs at I <~50 and 55 uA)
# runs 17048, 17049
#-----------------------------------

run_src = list(run)
T2_scl_per_mC_src = list(T2_per_mC_bcm4a_norm_new)
I_bcm4a_avg_src = list(I_bcm4a_avg_new)
Q_src = list(bcm4a_charge_new)
Q_csum_src = np.cumsum(Q_src)

T2_scl_per_mC_src = np.array(T2_scl_per_mC_src)
I_bcm4a_avg_src = np.array(I_bcm4a_avg_src)
Q_csum_src = np.array(Q_csum_src)

# Apply a correction factor to up-scale data to 55 uA beam current
m_slope1 = 1.712e-3  # slope of 1st 7 outlier runs
m_slope1_err = 1.120e-4

m_slope2 = 1.573e-3  # slope of 3 last Ca48 MF runs
m_slope2_err = 2.572e-6

# calculate weighted average
w1 = 1./m_slope1_err**2
w2 = 1./m_slope2_err**2

m_slope_avg = (m_slope1*w1 + m_slope2*w2)/(w1+w2)
m_slope_avg_err = 1./np.sqrt(w1 + w2)

# calculate un-weighted average
m_slope_avg_uw = (m_slope1 + m_slope2)/2.
m_slope_avg_uw_err = np.sqrt(0.5 * (m_slope1_err**2 + m_slope2_err))

#print(m_slope_avg )
#print(m_slope_avg_err )

print('unweighted avg slope = ', m_slope_avg_uw )
print('unweighted avg slope err =', m_slope_avg_uw_err )



# calculating corrections using separate slopes
T2_scl_per_mC_src_corr1 = T2_scl_per_mC_src[:-3] + (m_slope1*55. - m_slope1*I_bcm4a_avg_src[:-3])  # select all but last 3 runs (to apply correction 1)
T2_scl_per_mC_src_corr2 = T2_scl_per_mC_src[-3:] + (m_slope2*55. - m_slope2*I_bcm4a_avg_src[-3:])  # select only the last 3 runs (to apply correction 2)

T2_scl_per_mC_src_corr1_err = np.abs(55-I_bcm4a_avg_src[:-3]) * m_slope1_err
T2_scl_per_mC_src_corr2_err = np.abs(55-I_bcm4a_avg_src[-3:]) * m_slope2_err

T2_scl_corr = list(T2_scl_per_mC_src_corr1) + list(T2_scl_per_mC_src_corr2)
T2_scl_corr_err = list(T2_scl_per_mC_src_corr1_err) + list(T2_scl_per_mC_src_corr2_err)

print('T2_scl_per_mC_src_corr1 =',T2_scl_per_mC_src_corr1 )
print('T2_scl_per_mC_src_corr2 =',T2_scl_per_mC_src_corr2 )
print('T2_scl_corr =',T2_scl_corr )

# calculating corrections using the average of the slopes
T2_scl_per_mC_src_corr_avg = T2_scl_per_mC_src + (m_slope_avg*55. - m_slope_avg*I_bcm4a_avg_src)  
T2_scl_per_mC_src_corr_avg_err = np.abs(55-I_bcm4a_avg_src) * m_slope_avg_err

# calculating corrections using the average of the slopes (unweighted)
T2_scl_per_mC_src_corr_uw_avg = T2_scl_per_mC_src + (m_slope_avg_uw*55. - m_slope_avg_uw*I_bcm4a_avg_src)  
T2_scl_per_mC_src_corr_uw_avg_err = np.abs(55-I_bcm4a_avg_src) * m_slope_avg_uw_err
  

#print('Q=',Q_csum_src)
print('T2_scaler=',T2_scl_per_mC_src)
#print('T2_scaler_corr=',T2_scl_per_mC_src_corr)

print('Ib = ', I_bcm4a_avg_src)
Ib_corr  = [55] * len(I_bcm4a_avg_src)

# define line fit
def func(x, a, b, c):
    return a * np.exp(-b*x) + c

# fit the data
#popt, pcov = curve_fit(func, Q_csum_src, T2_scl_per_mC_src, p0=[1, 1e-05, 200])
#popt, pcov = curve_fit(func, Q_csum_src, T2_scl_corr, p0=[1, 1e-05, 200], sigma=T2_scl_corr_err)
popt, pcov = curve_fit(func, Q_csum_src[:-3], T2_scl_corr[:-3], p0=[1, 1e-05, 200], sigma=T2_scl_corr_err[:-3])  # exclude 3 last points in fit


#popt, pcov = curve_fit(func, Q_csum_src[:-3], T2_scl_per_mC_src_corr_avg[:-3], p0=[1, 1e-05, 200], sigma=T2_scl_per_mC_src_corr_avg_err[:-3])
#popt, pcov = curve_fit(func, Q_csum_src, T2_scl_per_mC_src_corr_uw_avg, p0=[1, 1e-05, 200], sigma=T2_scl_per_mC_src_corr_uw_avg_err)

p_sigma = np.sqrt(np.diag(pcov))


fig0, (ax1) = plt.subplots(1)
ax1.set_title('Charge-normalized T2 Scaler Counts Relative to 1st SRC Run', fontsize=16, fontweight='bold')

# beam current vs. run
#ax1.plot(run, I_bcm4a_avg_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Uncorrected')
#ax1.plot(run[:-3], Ib_corr[:-3], marker='^', color='red', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA')
#ax1.plot(run[-3:], Ib_corr[-3:], marker='^', color='green', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA' '\n' '(only last 3 runs)')

# T2 scaler/mC vs. run
#ax1.plot(run, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Uncorrected')
#ax1.plot(run[:-3], T2_scl_per_mC_src_corr1, marker='^', color='red', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA')
#ax1.plot(run[-3:], T2_scl_per_mC_src_corr2, marker='^', color='green', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA' '\n' '(only last 3 runs)')

#T2 scaler/mC vs. cumulative charge
ax1.plot(Q_csum_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Uncorrected')
ax1.errorbar(Q_csum_src[:-3], T2_scl_per_mC_src_corr1, T2_scl_per_mC_src_corr1_err, marker='^', color='red', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA')
ax1.errorbar(Q_csum_src[-3:], T2_scl_per_mC_src_corr2, T2_scl_per_mC_src_corr2_err, marker='^', color='green', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA' '\n' '(only last 3 runs)')

#ax1.errorbar(Q_csum_src, T2_scl_per_mC_src_corr_avg, T2_scl_per_mC_src_corr_avg_err, marker='^', color='red', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA' '\n' '(weighted avg. of slopes)')
#ax1.errorbar(Q_csum_src, T2_scl_per_mC_src_corr_uw_avg, T2_scl_per_mC_src_corr_uw_avg_err, marker='^', color='blue', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 Corrected to 55 uA' '\n' '(un-weighted avg. of slopes)')

# plot fit curve 
#ax1.plot(I_bcm4a_avg_src, func(I_bcm4a_avg_src, *popt), 'r-', label=r"fit: m$x$ + b" "\n" "m = %.3E $\pm$ %.3E" "\n" "b = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1]) )
#ax1.plot(Q_csum_src, func(Q_csum_src, *popt), 'r-', label=r"fit: $A e^{-\alpha Q}$ + C" "\n" "A = %.3E $\pm$ %.3E" "\n" r"$\alpha$ = %.3E $\pm$ %.3E" "\n" "C = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1], popt[2], p_sigma[2]) )
ax1.plot(Q_csum_src[:-3], func(Q_csum_src[:-3], *popt), 'r-', label=r"fit: $A e^{-\alpha Q}$ + C" "\n" "A = %.3E $\pm$ %.3E" "\n" r"$\alpha$ = %.3E $\pm$ %.3E" "\n" "C = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1], popt[2], p_sigma[2]) )

ax1.set_ylim([0.93, 1.02])
#ax1.set_ylim([25, 65])

ax1.grid()
ax1.set_ylabel('Relative T2 scalers/mC', fontsize=16)
#ax1.set_ylabel('Average Beam Current [uA]', fontsize=16)


ax1.set_xlabel(r'Cumulative Charge [mC]', fontsize=16)
#ax1.set_xlabel(r'Beam Current [$\mu$A]', fontsize=16)
#ax1.set_xlabel(r'Run Number', fontsize=16)

ax1.tick_params(axis='both', which='both', labelsize=15)



ax1.legend(loc='upper right',fontsize=12)


#-----

#ax2.plot(Q_csum_src, T2_scl_per_mC_src, marker='^', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 SRC data ($I_{b} ~ 50, 55) \mu$A')

#ax2.set_ylim([0.93, 1.02])
#ax2.grid()
#ax2.legend(loc='upper left',fontsize=12)
#ax2.set_ylabel('Relative T2 scalers/mC', fontsize=16)
#ax2.set_xlabel(r'Beam Charge [mC]', fontsize=16)
#ax2.tick_params(axis='both', which='both', labelsize=15)


'''
# PLOT T1 (SHMS 3/4) /T2 (SHMS EL-REAL) ratio vs. beam current for 17036 (1st Ca48 SRC run) and last 3 Ca-48 MF runs (17093, 17094, 17096)
ratio = [591833944. / 401471731., 221892312 / 160130254., 57575009./41111681., 495627331./335721439] #T1 SCALERS /T2 SCALERS
run = [17036, 17093, 17094, 17096]
Ib = [52.359, 28.606, 33.127, 59.726]
plt.plot(Ib, ratio, marker='o', color='r',  markerfacecolor='white', markersize=8, linestyle='None', label=r'Ca48 T1/T2 scalers' )
plt.xlabel('Beam Current [uA]')
plt.ylabel('SHMS T1 (3/4) / T2 scalers (EL-REAL)')
'''

plt.show()
