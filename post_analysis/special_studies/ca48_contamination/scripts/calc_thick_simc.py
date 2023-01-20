import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
import pandas as pd
from matplotlib import rc
from matplotlib import pyplot

rc('text', usetex=True)
plt.rcParams["font.family"] = "Times New Roman"

# using fit results: evaluate data counts / mC and extract thickness
def get_Hthick(N, Q):
    # y = mx + b -->  Ydata = m*hthick + b
    m = popt[0]  # slope
    m_err = p_sigma[0]
    b = popt[1]  # y-intercept
    b_err = p_sigma[1]

    Ydata = N / Q
    hthick = (Ydata - b)/m

    print('Data: ', N,'Counts', Q,'[mC]', Ydata,'[counts/mC]')
    print('H-thick: ', hthick, '[g/cm^2]')
    return [Ydata, hthick]  # [data counts/mC,  H-thickness [g/cm2] ]

#df_scaler = pd.read_csv("ca48_report_files/fall2021_bcm/ca48_bcmFall2021.csv")
df = pd.read_csv("simc_tgthick_variation.csv")

x     = df['thickness']/1000.  #g/cm^2
y     = df['N']          # ep elastic counts (with |Em| < 0.02 && Pm <0.03 cuts
y_err = df['N_err']      # error in counts

# define line fit
def func(x, m, b):
    return m * x + b

# fit the data
popt, pcov = curve_fit(func, x[:8], y[:8], sigma=y_err[:8], p0=[30, 0.1])
p_sigma = np.sqrt(np.diag(pcov))

fig0, (ax1) = plt.subplots(1)
ax1.set_title(r'$ep$ elastic counts vs. H-thickness ', fontsize=16, fontweight='bold')
ax1.errorbar(x[:8], y[:8], y_err[:8], marker='o', color='black', markerfacecolor='white', markersize=8, linestyle='None', label=r'SIMC')

ax1.set_xlim([0., 0.009])
ax1.set_ylim([0., 300])

ax1.grid()
ax1.legend(loc='upper left',fontsize=12)
ax1.set_ylabel(r'Counts / mC', fontsize=16)
ax1.set_xlabel(r'H-thickness [g/cm$^{2}$]', fontsize=16)
ax1.tick_params(axis='both', which='both', labelsize=15)


# plot fit curve 
ax1.plot(x[:8], func(x[:8], *popt), 'r-', label=r"fit: m$x$ + b" "\n" "m = %.3E $\pm$ %.3E" "\n" "b = %.3E $\pm$ %.3E" % (popt[0], p_sigma[0], popt[1], p_sigma[1]) )

# extract H-thickness from data (for each Ca48 MF run number)
Ydata_16978, hthick_16978 = get_Hthick(1156, 5.485)  #input: measured ep elastic counts, total charge [mC]
Ydata_16979, hthick_16979 = get_Hthick(12780, 76.223)
Ydata_17093, hthick_17093 = get_Hthick(1756, 49.725)
Ydata_17094, hthick_17094 = get_Hthick(377.8, 12.675)
Ydata_17096, hthick_17096 = get_Hthick(2914, 98.943)


plt.vlines(hthick_16978, 0, Ydata_16978, linestyle="dashed", color='k')
plt.hlines(Ydata_16978, 0, hthick_16978, linestyle="dashed", color='k', label=r'(data) 16978: %.3E [g/cm$^{2}$]'%(hthick_16978))

plt.vlines(hthick_16979, 0, Ydata_16979, linestyle="dashed", color='b')
plt.hlines(Ydata_16979, 0, hthick_16979, linestyle="dashed", color='b', label=r'(data) 16979: %.3E [g/cm$^{2}$]'%(hthick_16979))

plt.vlines(hthick_17093, 0, Ydata_17093, linestyle="dashed", color='g')
plt.hlines(Ydata_17093, 0, hthick_17093, linestyle="dashed", color='g', label=r'(data) 17093: %.3E [g/cm$^{2}$]'%(hthick_17093))

plt.vlines(hthick_17094, 0, Ydata_17094, linestyle="dashed", color='magenta')
plt.hlines(Ydata_17094, 0, hthick_17094, linestyle="dashed", color='magenta', label=r'(data) 17094: %.3E [g/cm$^{2}$]'%(hthick_17094))

plt.vlines(hthick_17096, 0, Ydata_17096, linestyle="dashed", color='r')
plt.hlines(Ydata_17096, 0, hthick_17096, linestyle="dashed", color='r', label=r'(data) 17096: %.3E [g/cm$^{2}$]'%(hthick_17096))

ax1.legend(loc='upper left',fontsize=15, bbox_to_anchor=(0.9, 0.9))

plt.show()




