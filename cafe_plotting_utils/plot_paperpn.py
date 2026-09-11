#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# This file when run asks for no inputs from the user but should be in the same folder as the appropriately formated .csv files from the ratiohistos.C code
# this file produces various plots used in my thesis and the nature paper
# Do note that many values are hard coded and really should be read from a .csv or .xlsv file
# This file has axis adjusted for per-nucleus normalization
#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=


import numpy as np
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
from matplotlib import colormaps
#from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#from matplotlib import pyplot as plt
#@import numpy as np


#6 plots,\
#light nucleu 
    #per proton
    #per Nucleus
#all nuclei
    ##per nuclues log log
    #per proton vs n/z data only
    #per proton vs A data only
    #per proton va A with theory
#Spatial (JAM) (jam)




ig, ax11 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax11.errorbar(
[0.385, 0.405, 0.435, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.6925],
[0.997, 1.046, 1.066, 1.101, 1.140, 1.137, 1.105, 1.116, 1.112, 1.306, 1.042, 1.140], yerr=
[0.050, 0.053, 0.046, 0.050, 0.054, 0.060, 0.061, 0.069, 0.075, 0.100, 0.084, 0.146], 
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', ecolor='r', linestyle='None', label='Ca48/Ca40', zorder=3)

#ax11.errorbar(
#[0.385, 0.405, 0.435, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.6925],
#[1.23811, 1.30667, 1.24314, 1.24473, 1.21468, 1.29108, 1.25105, 1.20387, 1.339, 1.43987, 1.1647, 1.19066], yerr=
#[0.0674707, 0.0708893, 0.0583572, 0.0622255, 0.0648018, 0.0748347, 0.0761684, 0.0831433, 0.0969575, 0.121598, 0.102991, 0.169913], 
#marker='P', markersize=15, alpha=1.0, mfc='g', mec='None', ecolor='g', linestyle='None', label='Fe54/Ca40', zorder=3)

ax11.errorbar(
[0.385, 0.405, 0.435, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.6925],
[1.614, 1.625, 1.516, 1.469, 1.385, 1.477, 1.472, 1.402, 1.565, 1.434, 1.452, 1.358], yerr=
[0.090, 0.089, 0.072, 0.074, 0.074, 0.085, 0.090, 0.097, 0.114, 0.117, 0.131, 0.194], 
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', ecolor='b', linestyle='None', label='Fe54/Ca48', zorder=3)

plt.hlines(y=[1.090, 1.488], xmin=0.3, xmax=0.8, colors=['r', 'b'], linestyles=['--', '--'])


ax11.set_title('', fontsize=21)
ax11.set_xlabel('P$_{miss}$ [GeV/c]', fontsize=21)
ax11.set_ylabel('Cross Section Ratios', fontsize=21)#, weight='bold')

ax11.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
ax11.set_xlim(0.3, 0.8)
ax11.set_ylim(0.0, 1.8)
plt.xticks(fontsize = 15)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
plt.yticks(fontsize = 15)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('PerNucleusRatio.png')
plt.close()













fig, ax12 = plt.subplots(figsize=(10,7))

# A / Ca40 vs. A (SRC)
ax12.errorbar(
[0.385, 0.405, 0.435, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.7],#0.6925
[0.997, 1.046, 1.066, 1.101, 1.140, 1.137, 1.105, 1.116, 1.112, 1.306, 1.042, 1.140], yerr=
[0.050, 0.053, 0.046, 0.050, 0.054, 0.060, 0.061, 0.069, 0.075, 0.100, 0.084, 0.146], 
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', ecolor='r', linestyle='None', label='Ca48/Ca40', zorder=3)

#ax12.errorbar(
#[0.385, 0.405, 0.435, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.7],#0.6925
#[1.23811, 1.30667, 1.24314, 1.24473, 1.21468, 1.29108, 1.25105, 1.20387, 1.339, 1.43987, 1.1647, 1.19066], yerr=
#[0.0674707, 0.0708893, 0.0583572, 0.0622255, 0.0648018, 0.0748347, 0.0761684, 0.0831433, 0.0969575, 0.121598, 0.102991, 0.169913], 
#marker='P', markersize=15, alpha=1.0, mfc='g', mec='None', ecolor='g', linestyle='None', label='Fe54/Ca40', zorder=3)

ax12.errorbar(
[0.385, 0.405, 0.435, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.7],#0.6925
[1.242, 1.250, 1.166, 1.130, 1.066, 1.136, 1.133, 1.079, 1.204, 1.103, 1.117, 1.045], yerr=
[0.069, 0.069, 0.055, 0.057, 0.057, 0.066, 0.069, 0.075, 0.088, 0.090, 0.101, 0.150], 
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', ecolor='b', linestyle='None', label='Fe54/Ca48', zorder=3)

plt.hlines(y=[1.090, 1.144], xmin=0.3, xmax=0.8, colors=['r', 'b'], linestyles=['--', '--'])

ax12.set_title('', fontsize=21)
ax12.set_xlabel('P$_{miss}$ [GeV/c]', fontsize=21)
ax12.set_ylabel('Cross Section Ratios', fontsize=21)#, weight='bold')

ax12.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.3, 0.4, 0.5, 0.6, 0.7])
ax12.set_xlim(0.3, 0.8)
ax12.set_ylim(0.0, 1.5)
plt.xticks(fontsize = 15)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
plt.yticks(fontsize = 15)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('PerProtonRatio.png')
plt.close()




























fig, ax13 = plt.subplots(figsize=(10,7))

ax13.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022], yerr=
#[0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', ecolor='r', linestyle='None', label='Ca48/Ca40', zorder=3)
#marker='s', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Ca48/Ca40', zorder=3)

ax13.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517], yerr=
#[0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', ecolor='b', linestyle='None', label='Fe54/Ca48', zorder=3)
#marker='o', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Fe54/Ca48', zorder=3)

#plt.hlines(y=[1.093, 1.474], xmin=0.3, xmax=0.8, colors=['r', 'b'], linestyles=['--', '--'])
plt.hlines(y=[1.093, 1.474], xmin=0.39125, xmax=0.68375, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])





x2=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y2=np.array([1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517])
error2 = np.array([0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141])
plt.plot(x2, y2, linestyle='none')
plt.fill_between(x2, y2-error2, y2+error2, alpha=0.2, edgecolor='b', facecolor='b')

x1=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y1=np.array([1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022])
error1 = np.array([0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085])
plt.plot(x1, y1, linestyle='none')
plt.fill_between(x1, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')


ax13.set_title('', fontsize=21)
ax13.set_xlabel('Missing Momentum [GeV/c]', fontsize=23)
ax13.set_ylabel('Cross Section Ratio', fontsize=23)#, weight='bold')

ax13.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax13.set_xlim(0.375, 0.725)
ax13.set_ylim(0.0, 1.8)
plt.xticks(fontsize = 20)
#plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75])
#plt.yticks([0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8])
plt.yticks([0.0, 0.5, 1.0, 1.5])

plt.yticks(fontsize = 20)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('PerNucleusRatio1.png')
plt.close()
















fig, ax14 = plt.subplots(figsize=(10,7))

ax14.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022], yerr=
#[0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', ecolor='r', linestyle='None', label='Ca48/Ca40', zorder=3)
#marker='s', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Ca48/Ca40', zorder=3)

ax14.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517], yerr=
#[0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', ecolor='b', linestyle='None', label='Fe54/Ca48', zorder=3)
#marker='o', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Fe54/Ca48', zorder=3)

#plt.hlines(y=[1.093, 1.474], xmin=0.3, xmax=0.8, colors=['r', 'b'], linestyles=['--', '--'])
plt.hlines(y=[1.093, 1.474], xmin=0.375, xmax=0.725, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])





x2=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y2=np.array([1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517])
error2 = np.array([0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141])
plt.plot(x2, y2, linestyle='none')
plt.fill_between(x2, y2-error2, y2+error2, alpha=0.2, edgecolor='b', facecolor='b')

x1=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y1=np.array([1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022])
error1 = np.array([0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085])
plt.plot(x1, y1, linestyle='none')
plt.fill_between(x1, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')


ax14.set_title('', fontsize=21)
ax14.set_xlabel('Missing Momentum [GeV/c]', fontsize=23)
ax14.set_ylabel('Cross Section Ratio', fontsize=23)#, weight='bold')

ax14.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax14.set_xlim(0.375, 0.725)
ax14.set_ylim(0.0, 1.8)
plt.xticks(fontsize = 20)
#plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75])
#plt.yticks([0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8])
plt.yticks([0.0, 0.5, 1.0, 1.5])

plt.yticks(fontsize = 20)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('PerNucleusRatio2.png')
plt.close()

















fig, ax15 = plt.subplots(figsize=(10,7))

ax15.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022], yerr=
#[0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', ecolor='r', linestyle='None', label='Ca48/Ca40', zorder=3)
#marker='s', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Ca48/Ca40', zorder=3)

ax15.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517], yerr=
#[0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', ecolor='b', linestyle='None', label='Fe54/Ca48', zorder=3)
#marker='o', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Fe54/Ca48', zorder=3)

#plt.hlines(y=[1.093, 1.474], xmin=0.3, xmax=0.8, colors=['r', 'b'], linestyles=['--', '--'])
plt.hlines(y=[1.093, 1.474], xmin=0.375, xmax=0.725, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[3,3])





x2=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y2=np.array([1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517])
error2 = np.array([0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141])
plt.plot(x2, y2, linestyle='none')
plt.fill_between(x2, y2-error2, y2+error2, alpha=0.2, edgecolor='b', facecolor='b')

x1=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y1=np.array([1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022])
error1 = np.array([0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085])
plt.plot(x1, y1, linestyle='none')
plt.fill_between(x1, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')


ax15.set_title('', fontsize=21)
ax15.set_xlabel('Missing Momentum [GeV/c]', fontsize=21)
ax15.set_ylabel('Cross Section Ratio', fontsize=21)#, weight='bold')

ax15.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax15.set_xlim(0.375, 0.725)
ax15.set_ylim(0.0, 1.8)
plt.xticks(fontsize = 18)
#plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75])
#plt.yticks([0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8])
plt.yticks([0.0, 0.5, 1.0, 1.5])

plt.yticks(fontsize = 18)

#plt.legend(frameon=False, fontsize=16, loc='lower right')

#fig.savefig('PerNucleusRatio3.png')
plt.close()





#0.00875, 0.03375, 0.06625, 0.09875, 0.13125, 0.16375, 0.19625, 0.22875, 0.26125, 0.29375, 0.32625, 0.35875, 0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.58625, 0.61875, 0.65125, 0.68375, 0.71625, 0.74875, 0.78125, 0.81375
#0, 0, 0, 0.025229, 0.050458, 0.0754922, 0.100937, 0.301323, 5.39065, 21.3759, 29.6345, 26.3195, 24.3768, 22.0904, 17.7611, 15.1022, 12.9531, 9.94315, 8.74402, 7.0931, 5.93054, 4.51153, 3.5039, 3.26146, 1.63015, 1.9062
#0, 0, 0, 0.025229, 0.0356792, 0.043586, 0.0504685, 0.0869892, 0.367662, 0.73195, 0.861669, 0.81192, 0.781548, 0.743883, 0.667082, 0.615055, 0.569713, 0.49907, 0.468092, 0.42167, 0.38526, 0.336293, 0.296155, 0.286062, 0.202207, 0.218668









fig, ax23 = plt.subplots(figsize=(10,7))

ax23.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022], yerr=
#[0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', ecolor='r', linestyle='None', label='Ca48/Ca40', zorder=3, capsize=0)
#marker='s', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Ca48/Ca40', zorder=3)

ax23.errorbar(
[0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375],
[1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517], yerr=
#[0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', ecolor='b', linestyle='None', label='Fe54/Ca48', zorder=3, capsize=0)
#marker='o', markersize=15, alpha=1.0, mfc='black', mec='None', ecolor='black', linestyle='None', label='Fe54/Ca48', zorder=3)

#plt.hlines(y=[1.093, 1.474], xmin=0.3, xmax=0.8, colors=['r', 'b'], linestyles=['--', '--'])
plt.hlines(y=[1.093, 1.474], xmin=0.39125, xmax=0.68375, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])

#momentum
#plt.hlines(y=[1.322412/1.16379, 1.16379], xmin=0.39125, xmax=0.68375, colors=['b', 'b'], linestyles=['--', '--'], linewidth=[5,5])

#Spatial (JAM)
#plt.hlines(y=[1.513/1.319, 1.319], xmin=0.39125, xmax=0.68375, colors=['r', 'r'], linestyles=['--', '--'], linewidth=[5,5])

#l=0, n=0
#plt.hlines(y=[1.2, 1.4662024/1.2], xmin=0.39125, xmax=0.68375, colors=['g', 'g'], linestyles=['--', '--'], linewidth=[5,5])

#L=0, n=0
#plt.hlines(y=[1.0, 1.43], xmin=0.39125, xmax=0.68375, colors=['g', 'g'], linestyles=['--', '--'], linewidth=[5,5])




x2=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y2=np.array([1.607, 1.511, 1.508, 1.352, 1.458, 1.463, 1.418, 1.574, 1.285, 1.517])
error2 = np.array([0.069, 0.069, 0.071, 0.069, 0.081, 0.087, 0.096, 0.110, 0.105, 0.141])
plt.plot(x2, y2, linestyle='none')
plt.fill_between(x2, y2-error2, y2+error2, alpha=0.2, edgecolor='b', facecolor='b')

x1=np.array([0.39125, 0.42375, 0.45625, 0.48875, 0.52125, 0.55375, 0.588625, 0.61875, 0.65125, 0.68375])
y1=np.array([1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022])
error1 = np.array([0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085])
plt.plot(x1, y1, linestyle='none')
plt.fill_between(x1, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')




#y5=y4/y3

#plt.plot(x4, y5, color='blue', linestyle='solid', linewidth = 3, marker='none', markerfacecolor='none', markersize=12)


plt.text(0.7, 1.474, r'$\dfrac{^{54} \mathrm{Fe} }{ ^{48} \mathrm{Ca} }$', fontsize=25, color='b')
plt.text(0.7, 1.093, r'$\dfrac{^{48} \mathrm{Ca} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color='r')




ax23.set_title('', fontsize=21)
ax23.set_xlabel('Missing Momentum (Pmiss) [GeV/c]', fontsize=25)
ax23.set_ylabel('Cross Section Ratio', fontsize=25)#, weight='bold')

ax23.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax23.set_xlim(0.38, 0.75)
ax23.set_ylim(0.0, 1.8)
plt.xticks(fontsize = 22)
#plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75])
#plt.yticks([0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8])
plt.yticks([0.0, 0.5, 1.0, 1.5])

plt.yticks(fontsize = 22)

#plt.legend(frameon=False, fontsize=16, loc='lower right')
fig.set_tight_layout(True)

fig.savefig('PerNucleusRatio6.png')
plt.close()























fig123, main_ax = plt.subplots(figsize=(10,7))

categories = ['$^{9}$Be', '$^{10}$B', '$^{11}$B', '$^{12}$C']
values = [9, 10, 11, 12]
categories2 = ['$^{40}$Ca', '$^{48}$Ca', '$^{54}$Fe']
values2 = [40, 48, 54]

inset_ax = main_ax.inset_axes(
   [0.05, 0.55, 0.38, 0.38],  # [x, y, width, height] w.r.t. axes
    xlim=[8.5, 12.5], ylim=[0.45, 1.1], # sets viewport &amp; tells relation to main axes
    yticks=[0.5, 0.7, 0.9, 1.1]#, xticks=[9, 10, 11, 12]
    #xticks(fontsize = 22)#, weight='bold')
    #plt.yticks(fontsize = 22)
)
#main_ax.indicate_inset_zoom(inset_ax, edgecolor="black")
inset_ax.set_xticks(values)
inset_ax.set_xticklabels(categories)
inset_ax.set_xticks(values, categories, fontsize = 15)#, weight='bold')
inset_ax.set_yticks([0.5, 0.7, 0.9, 1.1], [0.5, 0.7, 0.9, 1.1], fontsize = 15)#, weight='bold')
#inset_ax.set_yticks(fontsize = 22)#mark_inset(main_ax, inset_ax, loc1=3, loc2=4, fc="none", ec="0.5")
inset_ax.spines['bottom'].set_color('r')
inset_ax.spines['top'].set_color('r')
inset_ax.spines['right'].set_color('r')
inset_ax.spines['left'].set_color('r')
inset_ax.spines['bottom'].set_linewidth(1)
inset_ax.spines['top'].set_linewidth(1)
inset_ax.spines['right'].set_linewidth(1)
inset_ax.spines['left'].set_linewidth(1)
plt.hlines(y=0.44, xmin=8.0, xmax=13.0, color='r', linestyles='-', linewidth=1)
plt.vlines(x=13.0, ymin=0.44, ymax=1.12, color='r', linestyles='-', linewidth=1)
plt.hlines(y=1.12, xmin=8.0, xmax=13.0, color='r', linestyles='-', linewidth=1)
plt.vlines(x=8.0, ymin=0.44, ymax=1.12, color='r', linestyles='-', linewidth=1)

#inset_ax2 = main_ax.inset_axes(
#   [0.61, 0.07, 0.38, 0.38],  # [x, y, width, height] w.r.t. axes
#    xlim=[38, 56], ylim=[2, 8.5], # sets viewport &amp; tells relation to main axes
#    yticks=[2.0, 4.0, 6.0, 8.0]#, xticks=[40, 48, 54]
#)
#main_ax.indicate_inset_zoom(inset_ax2, edgecolor="black")
#inset_ax2.set_xticks(values2)
#inset_ax2.set_xticklabels(categories2)
#inset_ax2.set_xticks(values2, categories2, fontsize = 15)#, weight='bold')
#inset_ax2.set_yticks([2, 4, 6, 8], [2, 4, 6, 8], fontsize = 15)#, weight='bold')
#mark_inset(main_ax, inset_ax2, loc1=1, loc2=3, fc="none", ec="0.5")
#inset_ax2.spines['bottom'].set_color('b')
#inset_ax2.spines['top'].set_color('b')
#inset_ax2.spines['right'].set_color('b')
#inset_ax2.spines['left'].set_color('b')
#inset_ax2.spines['bottom'].set_linewidth(1)
#inset_ax2.spines['top'].set_linewidth(1)
#inset_ax2.spines['right'].set_linewidth(1)
#inset_ax2.spines['left'].set_linewidth(1)
#plt.hlines(y=2.5, xmin=36.0, xmax=60.0, color='b', linestyles='-', linewidth=1)
#plt.vlines(x=60.0, ymin=2.5, ymax=9.0, color='b', linestyles='-', linewidth=1)
#plt.hlines(y=9.0, xmin=36.0, xmax=60.0, color='b', linestyles='-', linewidth=1)
#plt.vlines(x=36.0, ymin=2.5, ymax=9.0, color='b', linestyles='-', linewidth=1)

inset_ax.errorbar(
[9, 10, 11, 12],
[0.53, 0.77, 0.83, 1.0],
[0.03, 0.03, 0.03, 0.000],
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label=''
)

#inset_ax2.errorbar(
#[40, 48, 54],
##[3.09, 3.39, 5.06],
#[0.17, 0.20, 0.28],
#marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label=''
#)


#plt.errorbar([9, 10, 11, 12, 40, 48, 54], [0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167], [0.023, 0.029, 0.030, 0.0, 0.038, 0.047, 0.050], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')
plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197], 
[0.53, 0.77, 0.83, 1.0, 3.09, 3.39, 5.06, 16.62], 
[0.03, 0.03, 0.03, 0.00, 0.17, 0.20, 0.28, 1.29], 
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='Data')

# PLOT MODELS (per  nucleus to C vs A)  * NEED TO convert to N/Z for the PRL TOP 1
plt.plot([9, 10, 11, 40], [0.623, 0.747, 0.843, 2.847], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
inset_ax.plot([9, 10, 11], [0.623, 0.747, 0.843], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.683, 0.811, 0.866, 3.860, 4.490, 5.108, 20.856], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
inset_ax.plot([9, 10, 11], [0.683, 0.811, 0.866], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='None', zorder=1)

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.674, 0.797, 0.866, 3.68, 4.131, 4.938, 17.801], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)
inset_ax.plot([9, 10, 11], [0.674, 0.797, 0.866], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)

plt.plot([9, 40, 48, 54, 197], [0.631, 5.267, 6.320, 7.722, 35.392], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
inset_ax.plot([9], [0.631], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.667, 0.667, 0.833, 3.333, 3.333, 4.333, 14.483], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
inset_ax.plot([9, 10, 11], [0.667, 0.667, 0.833], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.plot([9, 11, 40, 48, 54, 197], [0.607, 0.850, 4.100, 4.600, 5.373, 20.013], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
inset_ax.plot([9, 11], [0.607, 0.850], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.xscale('log')
plt.yscale('log')

plt.xticks(fontsize = 22)#, weight='bold')
plt.yticks(fontsize = 22)#, weight='bold')
plt.ylim(0.4, 40.0)#40
plt.xlim(6.5, 300.0)#na

plt.title('', fontsize=25)
plt.xlabel('A', fontsize=25)#, weight='bold')
plt.ylabel(r'$\frac{\sigma_{\mathrm{A}}}{\sigma_{\mathrm{C}}}$', fontsize=40, rotation='horizontal')#, weight='bold')
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])
plt.tight_layout()

plt.legend(frameon=False, fontsize=14, loc='lower right')

plt.savefig('SRCsingle_vs_A_all_c12_pn_inset.png')
plt.close()










fig17 = plt.figure(figsize=(10,7))

plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197], 
[0.534, 0.768, 0.834, 1.0, 3.09, 3.387, 5.055, 16.617], 
[0.026, 0.027, 0.029, 0.000, 0.165, 0.196, 0.277, 1.290], 
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PLOT MODELS
plt.plot([9, 10, 11, 40], [0.623, 0.747, 0.843, 2.847], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([9, 10, 11, 40, 48, 54, 197], [0.683, 0.811, 0.866, 3.860, 4.490, 5.108, 20.856], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([9, 10, 11, 40, 48, 54, 197], [0.52, 0.717, 0.85, 3.867, 5.1, 5.85, 22.778], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)
plt.plot([9, 40, 48, 54, 197], [0.631, 5.267, 6.320, 7.722, 35.392], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([9, 40, 48, 54, 197], [0.631, 5.267, 6.320, 7.722, 35.392], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
plt.plot([9, 10, 11, 40, 48, 54, 197], [0.667, 0.667, 0.833, 3.333, 3.333, 4.767, 14.483], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([9, 10, 11, 40, 48, 54, 197], [0.667, 0.667, 0.833, 3.333, 3.333, 4.767, 14.483], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.xscale('log')
plt.yscale('log')

plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.ylim(0.4, 40.0)

plt.title('', fontsize=18)
plt.xlabel('A', fontsize=21)#, weight='bold')
plt.ylabel(r'Per Nucleus Cross Section Ratios', fontsize=21)#, weight='bold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.legend(frameon=False, fontsize=14, loc='upper left')
plt.savefig('SRCsingle_vs_A_all_c12_pn.png')
plt.close()




fig16, main_ax = plt.subplots(figsize=(10,7))
categories = ['$^{9}_{4}$Be', '$^{10}_{5}$B', '$^{11}_{5}$B', '$^{12}_{6}$C']
values = [9, 10, 11, 12]

plt.errorbar(
[9, 10, 11, 12], 
[0.534, 0.768, 0.834, 1.0],
[0.026, 0.027, 0.029, 0.000],       marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PLOT MODELS
plt.plot([9, 10, 11], [0.623, 0.747, 0.843],   marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([9, 10, 11], [0.683, 0.811, 0.866],   marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([9, 10, 11], [0.52, 0.717, 0.85], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)
plt.plot([9], [0.631],   marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([9], [0.631],   marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
plt.plot([9, 10, 11], [0.667, 0.667, 0.833],   marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([9, 10, 11], [0.667, 0.667, 0.833],   marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

main_ax.set_xticks(values)
main_ax.set_xticklabels(categories)

plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.ylim(0.5, 1.025)#1.6
plt.xlim(8.75, 12.25)#1.6

#plt.title('', fontsize=18)
plt.xlabel('A', fontsize=21)#, weight='bold')
plt.ylabel(r'Per Nucleus Cross Section Ratios', fontsize=21)#, weight='bold')

plt.xticks(fontsize = 15)
#plt.xticks([9, 10, 11, 12])
plt.yticks(fontsize = 15)
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    
plt.legend(frameon=False, fontsize=14, loc='lower right')
fig16.savefig('SRCsingle_vs_A_light_c12_pn.png')
plt.close()















fig18= plt.figure(figsize=(10,7))
plt.errorbar([9, 10, 11, 12, 40, 48, 54, 197], 
[0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167, 1.262],
[0.038, 0.033,0.035, 0.000, 0.049, 0.059, 0.064, 0.098],
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PLOT MODELS
plt.plot([9, 10, 11, 40], [0.934, 0.896, 1.011, 0.854], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([9, 10, 11, 40, 48, 54, 197], [1.025, 0.974, 1.039, 1.158, 1.147, 1.179, 1.584], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([9, 10, 11, 40, 48, 54, 197], [0.78, 0.86, 1.02, 1.16, 1.53, 1.35, 1.73], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)
plt.plot([9, 40, 48, 54, 197], [0.947, 1.58, 1.896, 1.782, 2.688], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([9, 40, 48, 54, 197], [0.947, 1.58, 1.896, 1.782, 2.688], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
plt.plot([9, 10, 11, 40, 48, 54, 197], [1.0, 0.8, 1.0, 1.0, 1.0, 1.1, 1.1], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([9, 10, 11, 40, 48, 54, 197], [1.0, 0.8, 1.0, 1.0, 1.0, 1.1, 1.1], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.xscale('log')

plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.ylim(0.6, 2.0)
#plt.ylim(0.6, 3.0)

plt.title('', fontsize=18)
plt.xlabel('A', fontsize=21)#, weight='bold')
plt.ylabel(r'Per Proton Cross Section Ratios', fontsize=21)#, weight='bold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])

plt.legend(frameon=False, fontsize=14, loc='upper left')
fig18.savefig('SRCsingle_vs_A_all_c12_pp.png')
plt.close()






fig182= plt.figure(figsize=(10,7))
plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197], 
[0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167, 1.262], 
[0.038, 0.033,0.035, 0.000, 0.049, 0.059, 0.064, 0.098], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.xscale('log')

plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.ylim(0.7, 1.4)
#plt.ylim(0.6, 3.0)

plt.title('', fontsize=18)
plt.xlabel('A', fontsize=21)#, weight='bold')
plt.ylabel(r'Per Proton Cross Section Ratios', fontsize=21)#, weight='bold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])

fig182.savefig('SRCsingle_vs_A_all_c12_pp_nt.png')
plt.close()









fig18= plt.figure(figsize=(10,7))
plt.errorbar([1.25, 1.0, 1.2, 1.0, 1.0, 1.4, 1.077, 1.494], 
[0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167, 1.262], 
[0.038, 0.033,0.035, 0.000, 0.049, 0.059, 0.064, 0.098], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PLOT MODELS
plt.plot([1.25, 1.0, 1.2, 1.0], [0.934, 0.896, 1.011, 0.854], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([1.25, 1.0, 1.2, 1.0, 1.4, 1.077, 1.494], [1.025, 0.974, 1.039, 1.158, 1.147, 1.179, 1.584], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([1.25, 1.0, 1.2, 1.0, 1.4, 1.077, 1.494 ], [0.78, 0.86, 1.02, 1.16, 1.53, 1.35, 1.73], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)
plt.plot([1.25, 1.0, 1.4, 1.077, 1.494 ], [0.947, 1.58, 1.896, 1.782, 2.688], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([1.25, 1.0, 1.4, 1.077, 1.494 ], [0.947, 1.58, 1.896, 1.782, 2.688], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
plt.plot([1.25, 1.0, 1.2, 1.0, 1.4, 1.077, 1.494 ], [1.0, 0.8, 1.0, 1.0, 1.0, 1.1, 1.1], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([1.25, 1.0, 1.2, 1.0, 1.4, 1.077, 1.494 ], [1.0, 0.8, 1.0, 1.0, 1.0, 1.1, 1.1], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)


plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.xlim(0.9, 1.6)
#plt.ylim(0.6, 2.0)
plt.ylim(0.6, 3.0)

plt.title('', fontsize=18)
plt.xlabel('N/Z', fontsize=21)#, weight='bold')
plt.ylabel(r'Per Proton Cross Section Ratios', fontsize=21)#, weight='bold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.xticks([0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0])

plt.legend(frameon=False, fontsize=14, loc='upper left')
fig18.savefig('SRCsingle_vs_NZ_all_c12_pp.png')
plt.close()









fig181= plt.figure(figsize=(10,7))
plt.errorbar(
[1.25, 1.0, 1.2, 1.0, 1.0, 1.4, 1.077, 1.494 ], 
[0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167, 1.262], 
[0.038, 0.033,0.035, 0.000, 0.049, 0.059, 0.064, 0.098], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.text(1.27, 0.801, r'$^{9} \mathrm{Be}$', fontsize=15, color='black')
plt.text(1.02, 0.881, r'$^{10} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.22, 1.0, r'$^{11} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.02, 1.0, r'$^{12} \mathrm{C}$', fontsize=15, color='black')
plt.text(1.02, 0.927, r'$^{40} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(1.42, 1.016, r'$^{48} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(1.097, 1.167, r'$^{54} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(1.514, 1.262, r'$^{197} \mathrm{Au}$', fontsize=15, color='black')


plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.xlim(0.9, 1.6)
#plt.ylim(0.6, 2.0)
plt.ylim(0.7, 1.4)

plt.title('', fontsize=18)
plt.xlabel('N/Z', fontsize=21)#, weight='bold')
plt.ylabel(r'Per Proton Cross Section Ratios', fontsize=21)#, weight='bold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.xticks([0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])

fig181.savefig('SRCsingle_vs_NZ_all_c12_pp_nt.png')
plt.close()





fig19, main_ax = plt.subplots(figsize=(10,7))
categories = ['$^{9}_{4}$Be', '$^{10}_{5}$B', '$^{11}_{5}$B', '$^{12}_{6}$C']
values = [9, 10, 11, 12]

plt.errorbar(
[9, 10, 11, 12], 
[0.801, 0.921, 1.0, 1.0], 
[0.038, 0.033,0.035, 0.000],       marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PLOT MODELS
plt.plot([9, 10, 11], [0.934, 0.896, 1.011], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([9, 10, 11], [1.025, 0.974, 1.039], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([9, 10, 11], [0.78, 0.86, 1.02], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)
plt.plot([9], [0.947], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([9], [0.947], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
plt.plot([9, 10, 11], [1.0, 0.8, 1.0], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([9, 10, 11], [1.0, 0.8, 1.0], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)


main_ax.set_xticks(values)
main_ax.set_xticklabels(categories)


plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.ylim(0.7, 1.1)#1.6
plt.xlim(8.75, 12.25)#1.6

#plt.title('', fontsize=18)
plt.xlabel('A', fontsize=21)#, weight='bold')
plt.ylabel(r'Per Proton Cross Section Ratios', fontsize=21)#, weight='bold')

plt.xticks(fontsize = 15)
#plt.xticks([9, 10, 11, 12])
plt.yticks(fontsize = 15)
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1])
    
plt.legend(frameon=False, fontsize=14, loc='lower right')
fig19.savefig('SRCsingle_vs_A_light_c12_pp.png')
plt.close()













fig18= plt.figure(figsize=(10,7))
#0.006, 0.006, 0.006
plt.errorbar([12, 54, 197], [0.619, 0.577, 0.451], [0.0037, 0.0035, 0.0027], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')
#0.003, 0.003, 0.004
#plt.errorbar([12, 54, 197], [0.741, 0.734, 0.603], [0.0023, 0.0022, 0.0018], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

#plt.plot([9, 10, 11, 12, 40, 48, 54, 197], [1.0074, 1.0037, 1.0037, 1.0000, 0.9484, 0.9484, 0.926, 0.731], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
#plt.plot([9, 10, 11, 12, 40, 48, 54, 197], [1.0043, 1.0022, 1.0022, 1.0000, 0.970, 0.970, 0.957, 0.84], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18

plt.xscale('log')

plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.ylim(0.0, 1.0)
#plt.ylim(0.6, 3.0)

plt.title('', fontsize=18)
plt.xlabel('A', fontsize=21)#, weight='bold')
plt.ylabel(r'Radiative Correction Factor', fontsize=21)#, weight='bold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.yticks([0.6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

fig18.savefig('RadCorr_MF.png')
plt.close()
































fig, ax232 = plt.subplots(figsize=(10,7))


ax232.errorbar(
[9, 11, 12, 40, 48, 54, 197],
[1.14, 1.02, 1.00, 1.33, 1.36, 1.06, 1.20], 
yerr=[0.06, 0.04, 0.00, 0.07, 0.08, 0.06, 0.09], 
marker='o', markersize=15, alpha=1.0, mfc='pink', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
#ax232.errorbar(
#[9, 12, 40, 48, 54, 197],
#[1.18, 1.00, 1.71, 1.87, 1.53, 2.13], 
#yerr=[
#0.06, 0.00, 0.09, 0.11, 0.08, 0.17
#],  
#marker='^', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[1.18],
marker='o', markersize=15, alpha=1.0, mfc='pink', mec='black', linestyle='None', label='Spatial Overlap', zorder=3)


ax232.errorbar(
[9, 10, 11, 12, 40],
[1.17, 0.97, 1.01, 1, 0.92], 
yerr=[
0.06, 0.03, 0.04, 0.00, 0.05
], 
marker='s', markersize=15, alpha=1.0, mfc='blue', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[1.17], 
marker='s', markersize=15, alpha=1.0, mfc='blue', mec='black', linestyle='None', label='Momentum (Av18)', zorder=3)


#ax124.plot([48], [1.0], marker='v', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
ax232.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[1.28, 1.06, 1.04, 1.00, 1.25, 1.33, 1.01, 1.26], 
yerr=[
0.06, 0.04, 0.04, 0.00, 0.07, 0.08, 0.06, 0.10], 
marker='D', markersize=15, alpha=1.0, mfc='green', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
#ax232.errorbar(
#[9, 10, 11, 12, 40, 48, 54, 197],
#[1.28, 1.06, 1.04, 1.00, 1.25, 1.33, 1.01, 1.26], 
#yerr=[
#0.06, 0.04, 0.04, 0.00, 0.07, 0.08, 0.06, 0.10
#], 
#marker='D', markersize=15, alpha=0.2, mfc='black', mec='None', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[1.28], 
marker='D', markersize=15, alpha=1.0, mfc='green', mec='black', linestyle='None', label='Momentum (SRG)', zorder=3)



#ax124.plot([48], [1.0], marker='v', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
ax232.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[1.26, 1.04, 1.04, 1.00, 1.19, 1.22, 0.98, 1.07], 
yerr=[0.06, 0.04, 0.04, 0.00, 0.06, 0.07, 0.05, 0.08], 
marker='P', markersize=15, alpha=1.0, mfc='red', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[0.974],
marker='P', markersize=15, alpha=1.0, mfc='red', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)




#ax232.errorbar(
#[9, 10, 11, 12, 40, 48, 54, 197],
#[0.974, 0.934, 1.020, 1.000, 1.251, 1.506, 1.157, 1.371], 
#yerr=[
#0.047, 0.033, 0.036, 0.000, 0.067, 0.087, 0.063, 0.106
#], 
#marker='o', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
#plt.plot(
#[-9],
#[0.974],
#marker='o', markersize=15, alpha=1.0, mfc='None', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)

#ax232.errorbar(
#[9, 10, 11, 12, 40, 48, 54, 197],
#[1.0, 1.0, 1.0, 1.000, 1.0, 1.104, 1.025, 1.0], 
#yerr=[
#0.047, 0.033, 0.036, 0.000, 0.067, 0.087, 0.063, 0.106
#], 
#marker='o', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
#plt.plot(
#[-9],
#[0.974],
#marker='o', markersize=15, alpha=1.0, mfc='None', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)


#ax124.plot([48], [1.0], marker='v', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
ax232.errorbar(
[9, 12, 40, 48, 54, 197],
[1.18, 1.00, 1.71, 1.87, 1.53, 2.13], 
yerr=[
0.06, 0.00, 0.09, 0.11, 0.08, 0.17
], 
marker='^', markersize=15, alpha=1.0, mfc='orange', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
#ax232.errorbar(
#[9, 12, 40, 48, 54, 197],
#[1.18, 1.00, 1.71, 1.87, 1.53, 2.13], 
#yerr=[
#0.06, 0.00, 0.09, 0.11, 0.08, 0.17
#],  
#marker='^', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[1.18],
marker='^', markersize=15, alpha=1.0, mfc='orange', mec='black', linestyle='None', label='l=0, n=0', zorder=3)


#ax124.plot([48], [1.0], marker='v', markersize=15, alpha=1.0, mfc='lightsteelblue', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
ax232.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[1.25, 0.87, 1.00, 1.00, 1.08, 0.98, 0.94, 0.87], 
yerr=[
0.06, 0.03, 0.04, 0.00, 0.06, 0.06, 0.05, 0.07
], 
marker='v', markersize=15, alpha=1.0, mfc='black', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
#ax232.errorbar(
#[9, 10, 11, 12, 40, 48, 54, 197],
#[1.25, 0.87, 1.00, 1.00, 1.08, 0.98, 0.94, 0.87], 
#yerr=[
#0.06, 0.03, 0.04, 0.00, 0.06, 0.06, 0.05, 0.07
#],  
#marker='v', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9], 
[1.25],
marker='v', markersize=15, alpha=1.0, mfc='black', mec='black', linestyle='None', label='l=0, L=0', zorder=3)



"""
x1=np.array([9, 10, 11, 12, 40, 48, 54, 197])
y1=np.array([1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022])
error1 = np.array([
0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085
])
plt.plot(x1, y1, linestyle='none')
plt.fill_between(x1, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')
"""
plt.hlines(y=[1.0, 1.0], xmin=0.0, xmax=300, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])

ax232.set_title('', fontsize=21)
ax232.set_xlabel('A', fontsize=25)
ax232.set_ylabel('Model Ratio to Data', fontsize=25)#, weight='bold')

ax232.tick_params(axis='both', which='major', labelsize=15)

plt.xscale('log')
plt.xlim(7, 220.0)
ax232.set_ylim(0.75, 2.4)
plt.xticks(fontsize = 22)
plt.yticks([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4])
plt.yticks(fontsize = 22)

plt.legend(frameon=False, fontsize=20, loc='upper left')
fig.set_tight_layout(True)

fig.savefig('ModelToDataRatio.png')
plt.close()

















fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 7))


plt.subplot(1,2,1)
plt.errorbar(
[1.25, 1.0, 1.2, 1.0, 1.0, 1.4, 1.077, 1.494 ], 
[0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167, 1.262], 
[0.038, 0.033,0.035, 0.000, 0.049, 0.059, 0.064, 0.098], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.errorbar(
[1.0, 1.077, 1.154, 1.537], 
[1.0, 1.15, 1.36, 1.50], 
[0.0, 0.09, 0.08, 0.1], marker='o', markersize=15, mfc='None', ecolor='red', mec='red', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.text(1.27, 0.801, r'$^{9} \mathrm{Be}$', fontsize=15, color='black')
plt.text(1.02, 0.881, r'$^{10} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.22, 1.0, r'$^{11} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.02, 1.0, r'$^{12} \mathrm{C}$', fontsize=15, color='black')
plt.text(1.02, 0.927, r'$^{40} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(1.42, 1.016, r'$^{48} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(1.097, 1.167, r'$^{54} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(1.502, 1.266, r'$^{197} \mathrm{Au}$', fontsize=15, color='black')

plt.text(1.097, 1.107, r'$^{27} \mathrm{Al}$', fontsize=15, color='black')
plt.text(1.18, 1.35, r'$^{56} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(1.42, 1.48, r'$^{208} \mathrm{Pb}$', fontsize=15, color='black')

plt.xticks(fontsize = 18)#, weight='bold')
plt.yticks(fontsize = 18)#, weight='bold')
plt.xlim(0.975, 1.6)
#pltlim(0.6, 2.0)
plt.ylim(0.7, 1.55)

plt.title('', fontsize=18)
plt.xlabel('N/Z', fontsize=25)#, weight='bold')
#plt.ylabel(r'$\frac{\sigma_{\mathrm{A}}}{\sigma_{\mathrm{C}}}$', fontsize=40, rotation='horizontal')#, weight='bold')
#plt.ylabel(r'$(\frac{\sigma_{\mathrm{A}}}{Z}) / (\frac{\sigma_{\mathrm{C}}}{6})$', fontsize=25)#, weight='bold')
plt.ylabel(r'SRC Per Proton Ratio to C', fontsize=25)#, weight='bold')
#plt.ylabel(r'Per Proton Ratio to C', fontsize=25)#, weight='bold')
#plt.ylabel(r'$\frac{\frac{\sigma_{\mathrm{A}}}{Z}}{\frac{\sigma_{\mathrm{C}}}{6}}$', fontsize=30, rotation='horizontal')#, weight='bold')
#axes.yaxis.label.set_position((0.1, 0.5))
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5])


plt.subplot(1,2,2)
plt.errorbar([9, 10, 11, 12, 40, 48, 54, 197], 
[0.80, 0.92, 1.00, 1.00, 0.93, 1.02, 1.17, 1.26],
[0.04, 0.03, 0.04, 0.00, 0.05, 0.06, 0.06, 0.10],
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.errorbar(
[12, 27, 56, 208], 
[1.0, 1.15, 1.36, 1.50],  
[0.0, 0.09, 0.08, 0.1], 
 marker='o', markersize=15, mfc='None', ecolor='red', mec='red', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.xscale('log')

plt.xticks(fontsize = 22)#, weight='bold')
plt.yticks(fontsize = 22)#, weight='bold')
plt.ylim(0.7, 1.55)

plt.title('', fontsize=25)
plt.xlabel('A', fontsize=25)#, weight='bold')
#plt.ylabel(r'Per Proton Cross Section Ratios', fontsize=21)#, weight='bold')
#plt.xticks(fontsize = 15)
plt.yticks([])
#plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1..5])

plt.tight_layout()
#fig.subplots_adjust(left=0.1)
plt.subplots_adjust(wspace=0)

plt.savefig('sidebyside.png')
plt.close()






fig200= plt.figure(figsize=(10,7))
plt.errorbar(
[-0.003, -0.001, 0.001, 0.003, 0.005], 
[0.909, 0.946, 0.969, 0.954, 0.940], 
[0.009, 0.007, 0.006, 0.008, 0.009],
marker='s', markersize=0, mfc='k', ecolor='k', mec='k', elinewidth=3.0, capsize=5, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.hlines(y=0.909, xmin=-0.004, xmax=-0.002, color='black', linestyle='-', linewidth=1)
plt.hlines(y=0.946, xmin=-0.002, xmax=0.000, color='black', linestyle='-', linewidth=1)
plt.hlines(y=0.969, xmin=0.000, xmax=0.002, color='black', linestyle='-', linewidth=1)
plt.hlines(y=0.954, xmin=0.002, xmax=0.004, color='black', linestyle='-', linewidth=1)
plt.hlines(y=0.940, xmin=0.004, xmax=0.006, color='black', linestyle='-', linewidth=1)

plt.vlines(x=-0.004, ymin=0.000, ymax=0.909, color='black', linestyle='-', linewidth=1)
plt.vlines(x=-0.002, ymin=0.909, ymax=0.946, color='black', linestyle='-', linewidth=1)
plt.vlines(x=-0.00, ymin=0.946, ymax=0.969, color='black', linestyle='-', linewidth=1)
plt.vlines(x=0.002, ymin=0.954, ymax=0.969, color='black', linestyle='-', linewidth=1)
plt.vlines(x=0.004, ymin=0.940, ymax=0.954, color='black', linestyle='-', linewidth=1)
plt.vlines(x=0.006, ymin=0.000, ymax=0.940, color='black', linestyle='-', linewidth=1)

plt.hlines(y=0.952, xmin=-0.008, xmax=0.008, color='black', linestyle='--', linewidth=2)
plt.vlines(x=-0.002, ymin=0.000, ymax=0.946, color='r', linestyle='--', linewidth=2)

plt.xticks(fontsize = 15)#, weight='bold')
plt.yticks(fontsize = 15)#, weight='bold')
plt.xlim(-0.008, 0.008)
plt.ylim(0.0, 1.0)

plt.title('', fontsize=18)
plt.xlabel('SHMS Y\'$_{tar}$ (rad)', fontsize=21)#, weight='bold')
plt.ylabel(r'SHMS Y\'$_{tar}$ Ratio', fontsize=21)#, weight='bold')
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.xticks([-0.008, -0.006, -0.004, -0.002, 0.0, 0.002, 0.004, 0.006, 0.008])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

fig200.savefig('ShmsYtarRatioProtAbs.png')
plt.close()











































fig, axes2 = plt.subplots(nrows=2, ncols=1, figsize=(10, 14))

plt.subplot(2,1,1)

categories = ['$^{9}$Be', '$^{10}$B', '$^{11}$B', '$^{12}$C']
values = [9, 10, 11, 12]
categories2 = ['$^{40}$Ca', '$^{48}$Ca', '$^{54}$Fe']
values2 = [40, 48, 54]

inset_ax = main_ax.inset_axes(
   [0.05, 0.55, 0.38, 0.38],  # [x, y, width, height] w.r.t. axes
    xlim=[8.5, 12.5], ylim=[0.45, 1.1], # sets viewport &amp; tells relation to main axes
    yticks=[0.5, 0.7, 0.9, 1.1]#, xticks=[9, 10, 11, 12]
)
#main_ax.indicate_inset_zoom(inset_ax, edgecolor="black")
inset_ax.set_xticks(values)
inset_ax.set_xticklabels(categories)
#mark_inset(main_ax, inset_ax, loc1=3, loc2=4, fc="none", ec="0.5")
inset_ax.spines['bottom'].set_color('r')
inset_ax.spines['top'].set_color('r')
inset_ax.spines['right'].set_color('r')
inset_ax.spines['left'].set_color('r')
inset_ax.spines['bottom'].set_linewidth(1)
inset_ax.spines['top'].set_linewidth(1)
inset_ax.spines['right'].set_linewidth(1)
inset_ax.spines['left'].set_linewidth(1)
plt.hlines(y=0.44, xmin=8.0, xmax=13.0, color='r', linestyles='-', linewidth=1)
plt.vlines(x=13.0, ymin=0.44, ymax=1.12, color='r', linestyles='-', linewidth=1)
plt.hlines(y=1.12, xmin=8.0, xmax=13.0, color='r', linestyles='-', linewidth=1)
plt.vlines(x=8.0, ymin=0.44, ymax=1.12, color='r', linestyles='-', linewidth=1)

inset_ax2 = main_ax.inset_axes(
   [0.61, 0.07, 0.38, 0.38],  # [x, y, width, height] w.r.t. axes
    xlim=[38, 56], ylim=[2, 8.5], # sets viewport &amp; tells relation to main axes
    yticks=[2.0, 4.0, 6.0, 8.0]#, xticks=[40, 48, 54]
)
#main_ax.indicate_inset_zoom(inset_ax2, edgecolor="black")
inset_ax2.set_xticks(values2)
inset_ax2.set_xticklabels(categories2)
#mark_inset(main_ax, inset_ax2, loc1=1, loc2=3, fc="none", ec="0.5")
inset_ax2.spines['bottom'].set_color('b')
inset_ax2.spines['top'].set_color('b')
inset_ax2.spines['right'].set_color('b')
inset_ax2.spines['left'].set_color('b')
inset_ax2.spines['bottom'].set_linewidth(1)
inset_ax2.spines['top'].set_linewidth(1)
inset_ax2.spines['right'].set_linewidth(1)
inset_ax2.spines['left'].set_linewidth(1)
plt.hlines(y=2.5, xmin=36.0, xmax=60.0, color='b', linestyles='-', linewidth=1)
plt.vlines(x=60.0, ymin=2.5, ymax=9.0, color='b', linestyles='-', linewidth=1)
plt.hlines(y=9.0, xmin=36.0, xmax=60.0, color='b', linestyles='-', linewidth=1)
plt.vlines(x=36.0, ymin=2.5, ymax=9.0, color='b', linestyles='-', linewidth=1)

inset_ax.errorbar(
[9, 10, 11, 12],
[0.53, 0.77, 0.83, 1.00],
[0.03, 0.03, 0.03, 0.00],
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label=''
)

inset_ax2.errorbar(
[40, 48, 54],
[3.09, 3.39, 5.06],
[0.17, 0.20, 0.28],
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label=''
)


#plt.errorbar([9, 10, 11, 12, 40, 48, 54], [0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167], [0.023, 0.029, 0.030, 0.0, 0.038, 0.047, 0.050], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')
plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197], 
[0.53, 0.77, 0.83, 1.00, 3.09, 3.39, 5.06, 16.62], 
[0.03, 0.03, 0.03, 0.00, 0.17, 0.20, 0.28, 1.29], 
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PLOT MODELS
plt.plot([9, 10, 11, 40], [0.623, 0.747, 0.843, 2.847], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
inset_ax.plot([9, 10, 11], [0.623, 0.747, 0.843], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
inset_ax2.plot([40], [2.847], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.683, 0.811, 0.866, 3.860, 4.490, 5.108, 20.856], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='None', zorder=1)
plt.plot([9, 10, 11, 40, 48, 54, 197], [0.683, 0.811, 0.866, 3.860, 4.490, 5.108, 20.856], marker='D', markersize=15, alpha=0.2, mfc='black', mec='None', linestyle='None', label='Momentum (SRG)', zorder=1)
inset_ax.plot([9, 10, 11], [0.683, 0.811, 0.866], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='None', zorder=1)
inset_ax.plot([9, 10, 11], [0.683, 0.811, 0.866], marker='D', markersize=15, alpha=0.2, mfc='black', mec='None', linestyle='None', label='Momentum (SRG)', zorder=1)
inset_ax2.plot([40, 48, 54], [3.860, 4.490, 5.108], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='None', zorder=1)
inset_ax2.plot([40, 48, 54], [3.860, 4.490, 5.108], marker='D', markersize=15, alpha=0.2, mfc='black', mec='None', linestyle='None', label='Momentum (SRG)', zorder=1)

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.52, 0.717, 0.85, 3.867, 5.1, 5.85, 22.778], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)
inset_ax.plot([9, 10, 11], [0.52, 0.717, 0.85], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)
inset_ax2.plot([40, 48, 54], [3.867, 5.1, 5.85], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)

plt.plot([9, 40, 48, 54, 197], [0.631, 5.267, 6.320, 7.722, 35.392], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([9, 40, 48, 54, 197], [0.631, 5.267, 6.320, 7.722, 35.392], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
#inset_ax.plot([9], [0.631], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
inset_ax.plot([9], [0.631], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
#inset_ax2.plot([40, 48, 54], [5.267, 6.320, 7.722], marker='^', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
inset_ax2.plot([40, 48, 54], [5.267, 6.320, 7.722], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.667, 0.667, 0.833, 3.333, 3.333, 4.767, 14.483], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([9, 10, 11, 40, 48, 54, 197], [0.667, 0.667, 0.833, 3.333, 3.333, 4.767, 14.483], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
inset_ax.plot([9, 10, 11], [0.667, 0.667, 0.833], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
inset_ax.plot([9, 10, 11], [0.667, 0.667, 0.833], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)
inset_ax2.plot([40, 48, 54], [3.333, 3.333, 4.767], marker='v', markersize=15, alpha=0.2, mfc='black', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
inset_ax2.plot([40, 48, 54], [3.333, 3.333, 4.767], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.xscale('log')
plt.yscale('log')

plt.xticks(fontsize = 22)#, weight='bold')
plt.yticks(fontsize = 22)#, weight='bold')
plt.ylim(0.4, 40.0)#40
plt.xlim(6.5, 300.0)#na

plt.title('', fontsize=25)
plt.xlabel('', fontsize=22)#, weight='bold')
plt.ylabel(r'$\frac{\sigma_{\mathrm{A}}}{\sigma_{\mathrm{C}}}$', fontsize=40, rotation='horizontal')#, weight='bold')
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)










plt.subplot(2,1,2)


plt.errorbar(
[9, 10, 11, 12, 40],
[1.166, 0.97, 1.01, 1, 0.92], 
yerr=[
0.056, 0.035, 0.035, 0.000, 0.049
], 
marker='s', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[1.166], 
marker='s', markersize=15, alpha=1.0, mfc='None', mec='black', linestyle='None', label='Momentum (Av18)', zorder=3)



plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[1.28, 1.057, 1.04, 1, 1.25, 1.33, 1.01, 1.26], 
yerr=[
0.061, 0.038, 0.036, 0.000, 0.067, 0.077, 0.055, 0.097
], 
marker='D', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[1.28, 1.057, 1.04, 1, 1.25, 1.33, 1.01, 1.26], 
yerr=[
0.061, 0.038, 0.036, 0.000, 0.067, 0.077, 0.055, 0.097
], 
marker='D', markersize=15, alpha=0.2, mfc='black', mec='None', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[1.28], 
marker='D', markersize=15, alpha=0.2, mfc='black', mec='None', linestyle='None', label='Momentum (SRG)', zorder=3)


plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[0.974, 0.934, 1.020, 1.000, 1.251, 1.506, 1.157, 1.371], 
yerr=[
0.047, 0.033, 0.036, 0.000, 0.067, 0.087, 0.063, 0.106
], 
marker='o', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[0.974],
marker='o', markersize=15, alpha=1.0, mfc='None', mec='black', linestyle='None', label='Spatial (JAM)', zorder=3)


#ax232.errorbar(
#[9, 12, 40, 48, 54, 197],
#[1.183, 1.000, 1.705, 1.866, 1.528, 2.130], 
#yerr=[
#0.057, 0.000, 0.091, 0.108, 0.084, 0.165
#], 
#marker='^', markersize=15, alpha=0.2, mfc='black', mec='None', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.errorbar(
[9, 12, 40, 48, 54, 197],
[1.183, 1.000, 1.705, 1.866, 1.528, 2.130], 
yerr=[
0.057, 0.000, 0.091, 0.108, 0.084, 0.165
], 
marker='^', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9],
[1.183],
marker='^', markersize=15, alpha=1.0, mfc='None', mec='black', linestyle='None', label='l=0, n=0', zorder=3)



plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[1.249, 0.869, 1.000, 1.000, 1.079, 0.984, 0.943, 0.872], 
yerr=[
0.060, 0.031, 0.035, 0.000, 0.058, 0.057, 0.052, 0.068
], 
marker='v', markersize=15, alpha=0.2, mfc='black', mec='None', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197],
[1.249, 0.869, 1.000, 1.000, 1.079, 0.984, 0.943, 0.872], 
yerr=[
0.060, 0.031, 0.035, 0.000, 0.058, 0.057, 0.052, 0.068
], 
marker='v', markersize=15, alpha=1.0, mfc='None', mec='black', ecolor = 'k', linestyle='None', label='', zorder=3, capsize=5)
plt.plot(
[-9], 
[1.249],
marker='v', markersize=15, alpha=0.2, mfc='black', mec='None', linestyle='None', label='l=0, L=0', zorder=3)
"""
x1=np.array([9, 10, 11, 12, 40, 48, 54, 197])
y1=np.array([1.025, 1.046, 1.108, 1.147, 1.124, 1.147, 1.049, 1.222, 1.258, 1.022])
error1 = np.array([
0.040, 0.042, 0.047, 0.052, 0.057, 0.062, 0.063, 0.081, 0.093, 0.085
])
plt.plot(x1, y1, linestyle='none')
plt.fill_between(x1, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')
"""
plt.hlines(y=[1.0, 1.0], xmin=0.0, xmax=300, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])

#plt.set_title('', fontsize=21)
plt.xlabel('A', fontsize=25)
plt.ylabel('Model Ratio to Data', fontsize=25)#, weight='bold')

plt.tick_params(axis='both', which='major', labelsize=15)

plt.xscale('log')
plt.xlim(7, 220.0)
ax232.set_ylim(0.75, 2.4)
plt.xticks(fontsize = 22)
plt.yticks([0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4])
plt.yticks(fontsize = 22)

plt.legend(frameon=False, fontsize=16, loc='upper left')













plt.tight_layout()
plt.savefig('testestestest.png')
plt.close()




















fig, axes = plt.subplots(figsize=(10, 7))

#A vs N/Z
plt.plot(
[9.0, 10.0, 11.0, 12.0, 40.0, 48.0, 54.0, 197.0],
[1.25, 1.0, 1.2, 1.0, 1.0, 1.4, 1.077, 1.494],  
marker='s', markersize=10, mfc='k', mec='k', linestyle='None', zorder=1, label='')

plt.plot(
[12.0, 27.0, 56.0, 208.0],
[1.0, 1.077, 1.154, 1.537],  
marker='o', markersize=15, mfc='None', mec='red', linestyle='None', zorder=1, label='')

plt.text(9.0, 1.27, r'$^{9} \mathrm{Be}$', fontsize=15, color='black')
plt.text(10.0, 1.02, r'$^{10} \mathrm{B}$', fontsize=15, color='black')
plt.text(11.0, 1.22, r'$^{11} \mathrm{B}$', fontsize=15, color='black')
plt.text(12.0, 1.02, r'$^{12} \mathrm{C}$', fontsize=15, color='black')
plt.text(40.0, 1.02, r'$^{40} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(48.0, 1.42, r'$^{48} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(54.0, 1.097, r'$^{54} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(197.0, 1.44, r'$^{197} \mathrm{Au}$', fontsize=15, color='black')

plt.text(27.0, 1.097, r'$^{27} \mathrm{Al}$', fontsize=15, color='black')
plt.text(56.0, 1.18, r'$^{56} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(208.0, 1.55, r'$^{208} \mathrm{Pb}$', fontsize=15, color='black')

plt.xticks(fontsize = 18)#, weight='bold')
plt.xscale('log')
plt.xlim(7.0, 280)
plt.xlabel('A', fontsize=25)#, weight='bold')

plt.yticks(fontsize = 18)#, weight='bold')
plt.ylim(0.7, 1.6)
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
plt.ylabel('N/Z', fontsize=25)#, weight='bold')

plt.title('', fontsize=18)

plt.savefig('AvNZ.png')
plt.close()
































fig, axes = plt.subplots(figsize=(10, 7))

plt.errorbar(
[1.25, 1.0, 1.2, 1.0, 1.0, 1.4, 1.077, 1.494 ], 
[0.801, 0.921, 1.0, 1.0, 0.927, 1.016, 1.167, 1.262], 
[0.038, 0.033,0.035, 0.000, 0.049, 0.059, 0.064, 0.098], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.errorbar(
[1.0, 1.077, 1.154, 1.537], 
[1.0, 1.15, 1.36, 1.50], 
[0.0, 0.09, 0.08, 0.1], marker='o', markersize=15, mfc='None', ecolor='red', mec='red', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.text(1.27, 0.801, r'$^{9} \mathrm{Be}$', fontsize=15, color='black')
plt.text(1.02, 0.881, r'$^{10} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.22, 1.0, r'$^{11} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.02, 1.0, r'$^{12} \mathrm{C}$', fontsize=15, color='black')
plt.text(1.02, 0.927, r'$^{40} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(1.42, 1.016, r'$^{48} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(1.097, 1.18, r'$^{54} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(1.502, 1.266, r'$^{197} \mathrm{Au}$', fontsize=15, color='black')

plt.text(1.097, 1.12, r'$^{27} \mathrm{Al}$', fontsize=15, color='black')
plt.text(1.18, 1.35, r'$^{56} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(1.42, 1.48, r'$^{208} \mathrm{Pb}$', fontsize=15, color='black')

plt.xticks(fontsize = 18)#, weight='bold')
plt.yticks(fontsize = 18)#, weight='bold')
plt.xlim(0.975, 1.6)
#pltlim(0.6, 2.0)
plt.ylim(0.7, 1.6)

plt.title('', fontsize=18)
plt.xlabel('N/Z', fontsize=25)#, weight='bold')
#plt.ylabel(r'$\frac{\sigma_{\mathrm{A}}}{\sigma_{\mathrm{C}}}$', fontsize=40, rotation='horizontal')#, weight='bold')
#plt.ylabel(r'$(\frac{\sigma_{\mathrm{A}}}{Z}) / (\frac{\sigma_{\mathrm{C}}}{6})$', fontsize=25)#, weight='bold')
plt.ylabel(r'SRC Per Proton Ratio to C', fontsize=25)#, weight='bold')
#plt.ylabel(r'Per Proton Ratio to C', fontsize=25)#, weight='bold')
#plt.ylabel(r'$\frac{\frac{\sigma_{\mathrm{A}}}{Z}}{\frac{\sigma_{\mathrm{C}}}{6}}$', fontsize=30, rotation='horizontal')#, weight='bold')
#axes.yaxis.label.set_position((0.1, 0.5))
#plt.yticks([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])

#plt.plot([0,1.6], [0,1.6], color='black', linestyle='--', linewidth=2, marker='o', alpha=0.5)


plt.tight_layout()

plt.savefig('sidebyside1.png')
plt.close()



# PRL Fig 2 (Top) per nucleus vs. N/Z
fig, axes = plt.subplots(figsize=(10, 7))

# per proton to C
#plt.errorbar(
#[1.25, 1.0, 1.2, 1.0, 1.0 ], 
#[0.801, 0.921, 1.0, 1.0, 0.927], 
#[0.038, 0.033,0.035, 0.000, 0.049], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PER NUCLEUS TO C
# Be9, B10, B11, C12
plt.errorbar(
[1.25, 1.0, 1.2, 1.0], 
[0.53, 0.77, 0.83, 1.00], 
[0.03, 0.03, 0.03, 0.00], marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')


# per nucleus to C (MEYTAL)
plt.errorbar(
[1.0], 
[1.0], 
[0.0], marker='o', markersize=15, mfc='None', ecolor='red', mec='red', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# PLOT MODELS (SRC Single Rati per  nucleus to C vs A or N/Z)  for the PRL TOP 1
# Be9 B10, B11, C12 --> A: 9, 10, 11, 12,  N/Z: 1.25, 1, 1.2, 1

# --- av18 ---
#plt.plot([9, 10, 11], [0.623, 0.747, 0.843], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([1.25, 1, 1.2], [0.623, 0.747, 0.843], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18

#plt.plot([9, 10, 11], [0.683, 0.811, 0.866], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([1.25, 1, 1.2], [0.683, 0.811, 0.866], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)

#plt.plot([9, 10, 11], [0.674, 0.797, 0.866], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)
plt.plot([1.25, 1, 1.2], [0.674, 0.797, 0.866], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)

#plt.plot([9], [0.631], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([1.25], [0.631], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

#plt.plot([9, 10, 11], [0.667, 0.667, 0.833], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([1.25, 1, 1.2], [0.667, 0.667, 0.833], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)

#plt.plot([9, 11], [0.607, 0.850], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
plt.plot([1.25, 1.2], [0.607, 0.850], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)


plt.text(1.27, 0.801, r'$^{9} \mathrm{Be}$', fontsize=15, color='black')
plt.text(1.02, 0.881, r'$^{10} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.22, 1.0, r'$^{11} \mathrm{B}$', fontsize=15, color='black')
plt.text(1.02, 1.0, r'$^{12} \mathrm{C}$', fontsize=15, color='black')


plt.xticks(fontsize = 18)#, weight='bold')
plt.yticks(fontsize = 18)#, weight='bold')

plt.xlim(0.975, 1.35)
plt.ylim(0.4, 1.2)

plt.title('', fontsize=18)
plt.xlabel('N/Z', fontsize=25)#, weight='bold')
plt.ylabel(r'SRC Per Nucleus Ratio to C', fontsize=25)#, weight='bold')

plt.xticks([1.0, 1.1, 1.2, 1.3])
plt.yticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])

#plt.plot([0,1.6], [0,1.6], color='black', linestyle='--', linewidth=2, marker='o', alpha=0.5)


plt.tight_layout()
plt.legend()

plt.savefig('prl_fig2_top.png')
plt.close()























fig, axes = plt.subplots(figsize=(10, 7))

plt.title('', fontsize=18)
plt.ylabel(r'SRC Per Proton Ratio to C', fontsize=25)#, weight='bold')

plt.errorbar([9, 10, 11, 12, 40, 48, 54, 197], 
[0.80, 0.92, 1.00, 1.00, 0.93, 1.02, 1.17, 1.26],
[0.04, 0.03, 0.04, 0.00, 0.05, 0.06, 0.06, 0.10],
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.errorbar(
[12, 27, 56, 208], 
[1.0, 1.15, 1.36, 1.50],  
[0.0, 0.09, 0.08, 0.1], 
 marker='o', markersize=15, mfc='None', ecolor='red', mec='red', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

plt.xscale('log')

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
plt.ylim(0.7, 1.6)

plt.title('', fontsize=25)
plt.xlabel('A', fontsize=25)
plt.yticks([])
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])

plt.tight_layout()

plt.savefig('sidebyside2.png')
plt.close()




#  Fig 2 PRL (bottom) Single Ratios
fig, axes = plt.subplots(figsize=(10, 7))

plt.title('', fontsize=18)
plt.ylabel(r'SRC Per Nucleus Ratio to C', fontsize=25)#, weight='bold')

# per proton to C
#plt.errorbar([9, 10, 11, 12, 40, 48, 54, 197], 
#[0.80, 0.92, 1.00, 1.00, 0.93, 1.02, 1.17, 1.26],
#[0.04, 0.03, 0.04, 0.00, 0.05, 0.06, 0.06, 0.10],
#marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# per nucleus to C
plt.errorbar([9, 10, 11, 12, 40, 48, 54, 197], 
[0.53, 0.77, 0.83, 1.00, 3.09, 3.39, 5.06, 16.62],
[0.03, 0.03, 0.03, 0.00, 0.17, 0.20, 0.28, 1.29],
marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')


# Meytal per nucleus to C ( conversion: Be: (per_proton X 4) /6, B10, B11: (per_proton X 5) /6, Ca40, 48: (per_proton X 20) /6, Fe54: (per_proton x26)/6, Au: (per_proton x 79/6)
plt.errorbar(
[12], 
[1],  
[0], 
 marker='o', markersize=15, mfc='None', ecolor='red', mec='red', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')


# PLOT MODELS (per  nucleus to C vs A)  
plt.plot([9, 10, 11, 40], [0.623, 0.747, 0.843, 2.847], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
inset_ax.plot([9, 10, 11], [0.623, 0.747, 0.843], marker='s', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.683, 0.811, 0.866, 3.860, 4.490, 5.108, 20.856], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
inset_ax.plot([9, 10, 11], [0.683, 0.811, 0.866], marker='D', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='None', zorder=1)

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.674, 0.797, 0.866, 3.68, 4.131, 4.938, 17.801], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)
inset_ax.plot([9, 10, 11], [0.674, 0.797, 0.866], marker='P', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)

plt.plot([9, 40, 48, 54, 197], [0.631, 5.267, 6.320, 7.722, 35.392], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
inset_ax.plot([9], [0.631], marker='^', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.plot([9, 10, 11, 40, 48, 54, 197], [0.667, 0.667, 0.833, 3.333, 3.333, 4.333, 14.483], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
inset_ax.plot([9, 10, 11], [0.667, 0.667, 0.833], marker='v', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)

plt.plot([9, 11, 40, 48, 54, 197], [0.607, 0.850, 4.100, 4.600, 5.373, 20.013], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='Spatial', zorder=3)
inset_ax.plot([9, 11], [0.607, 0.850], marker='o', markersize=15, alpha=1.0, mfc='none', mec='black', linestyle='None', label='', zorder=3)


plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
plt.xlim(8, 14)
plt.ylim(0.18, 1.22)

plt.title('', fontsize=25)
plt.xlabel('A', fontsize=25)
plt.yticks([0.2, 0.4, 0.8, 0.8, 1.0, 1.2])

plt.tight_layout()
plt.legend()
plt.savefig('prl_fig2_bottom.png')
plt.close()











fig, axes = plt.subplots(figsize=(10, 7))

#A vs N/Z
plt.plot(
[9.0, 10.0, 11.0, 12.0, 40.0, 48.0, 54.0, 197.0],
[1.25, 1.0, 1.2, 1.0, 1.0, 1.4, 1.077, 1.494],  
marker='s', markersize=10, mfc='k', mec='k', linestyle='None', zorder=1, label='')

plt.plot(
[12.0, 27.0, 56.0, 208.0],
[1.0, 1.077, 1.154, 1.537],  
marker='o', markersize=15, mfc='None', mec='red', linestyle='None', zorder=1, label='')

plt.text(9.0, 1.27, r'$^{9} \mathrm{Be}$', fontsize=15, color='black')
plt.text(10.0, 1.02, r'$^{10} \mathrm{B}$', fontsize=15, color='black')
plt.text(11.0, 1.22, r'$^{11} \mathrm{B}$', fontsize=15, color='black')
plt.text(12.0, 1.02, r'$^{12} \mathrm{C}$', fontsize=15, color='black')
plt.text(40.0, 1.02, r'$^{40} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(48.0, 1.42, r'$^{48} \mathrm{Ca}$', fontsize=15, color='black')
plt.text(54.0, 1.097, r'$^{54} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(197.0, 1.44, r'$^{197} \mathrm{Au}$', fontsize=15, color='black')

plt.text(27.0, 1.097, r'$^{27} \mathrm{Al}$', fontsize=15, color='black')
plt.text(56.0, 1.18, r'$^{56} \mathrm{Fe}$', fontsize=15, color='black')
plt.text(208.0, 1.55, r'$^{208} \mathrm{Pb}$', fontsize=15, color='black')

plt.xticks(fontsize = 18)#, weight='bold')
plt.xscale('log')
plt.xlim(7.0, 280)
plt.xlabel('A', fontsize=25)#, weight='bold')

plt.yticks(fontsize = 18)#, weight='bold')
plt.ylim(0.7, 1.6)
plt.yticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
plt.ylabel('N/Z', fontsize=25)#, weight='bold')

plt.title('', fontsize=18)

plt.savefig('AvNZ.png')
plt.close()
