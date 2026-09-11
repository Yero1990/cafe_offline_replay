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
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



'''
#---------PRL FIG 1, Per nucleus to C vs Pmiss--------------

ifname= 'SRCBe9Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y01 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10Pm.csv'
df = pd.read_csv(ifname, comment='#')
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
y04 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

for index, value in enumerate(y01):
    if value == 0:
        y01[index] = 0.1
for index, value in enumerate(e1):
    if value == 0:
        e1[index] = 0.1
for index, value in enumerate(y02):
    if value == 0:
        y02[index] = 0.1
for index, value in enumerate(e2):
    if value == 0:
        e2[index] = 0.1
for index, value in enumerate(y03):
    if value == 0:
        y03[index] = 0.1
for index, value in enumerate(e3):
    if value == 0:
        e3[index] = 0.1
for index, value in enumerate(y04):
    if value == 0:
        y04[index] = 0.1
for index, value in enumerate(e4):
    if value == 0:
        e4[index] = 0.1

for int in range(0,3):
    x = x[:-1]
    y01 = y01[:-1]
    e1 = e1[:-1]
    y02 = y02[:-1]
    e2 = e2[:-1]
    y03 = y03[:-1]
    e3 = e3[:-1]
    y04 = y04[:-1]
    e4 = e4[:-1]

for int in range(0,11):
    x = x[1:]
    y01 = y01[1:]
    e1 = e1[1:]
    y02 = y02[1:]
    e2 = e2[1:]
    y03 = y03[1:]
    e3 = e3[1:]
    y04 = y04[1:]
    e4 = e4[1:]

yBe9 = np.divide(y01, y04)
yB10 = np.divide(y02, y04)
yB11 = np.divide(y03, y04)
eBe9 = np.sqrt( ( (1/y04) * e1 )**2 + ( (y01/y04**2) * e4)**2)
eB10 = np.sqrt( ( (1/y04) * e2 )**2 + ( (y02/y04**2) * e4)**2)
eB11 = np.sqrt( ( (1/y04) * e3 )**2 + ( (y03/y04**2) * e4)**2)
nBe9 = 1 / (eBe9**2) * yBe9
nB10 = 1 / (eB10**2) * yB10
nB11 = 1 / (eB11**2) * yB11
dBe9 = 1 / (eBe9**2)
dB10 = 1 / (eB10**2)
dB11 = 1 / (eB11**2)
vBe9 = sum(nBe9) / sum(dBe9)
print(vBe9)
vB10 = sum(nB10) / sum(dB10)
vB11 = sum(nB11) / sum(dB11)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yBe9, marker='s', markersize=15, alpha=1.0, mfc='red', mec='None', linestyle='None', label=r'$^{9}$Be', zorder=3)
ax2091.plot(x, yB10, marker='o', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{10}$B', zorder=3)
ax2091.plot(x, yB11, marker='^', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{11}$B', zorder=3)

#plt.hlines(y=[vBe9, vB10, vB11], xmin=0.3, xmax=0.75, colors=['black', 'black', 'black'], linestyles=['--', '--', '--'], linewidth=[5,5,5])

#plt.text(0.7, 1.474, r'$\dfrac{^{9} \mathrm{Be} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')
#plt.text(0.7, 1.274, r'$\dfrac{^{10} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')
#plt.text(0.7, 1.074, r'$\dfrac{^{11} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor='green', facecolor='green')
plt.plot(x, yB11, linestyle='none')
plt.fill_between(x, yB11-eB11, yB11+eB11, alpha=0.5, edgecolor='blue', facecolor='blue')

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('SRC Per Nucleus Ratio to C', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax2091.set_xlim(0.36, 0.75)
#ax2091.set_ylim(0.35, 1.05)
ax2091.set_ylim(0., 1.05)

plt.xticks(fontsize = 22)
#plt.yticks([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0])

plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('prl_fig1_pn_to_C12_vs_pmiss_ratio.pdf')
plt.close()

#---------END: PRL FIG 1, Per nucleus to C vs Pmiss--------------
'''






#  Fig 2 PRL (bottom) Single Ratios
fig, axes = plt.subplots(figsize=(10, 7))

plt.title('', fontsize=18)
plt.ylabel(r'SRC Per Nucleus Ratio to C', fontsize=25)#, weight='bold')
#plt.ylabel(r'SRC Per Proton Ratio to C', fontsize=25)#, weight='bold')

# per proton to C
#plt.errorbar([9, 10, 11, 12, 40, 48, 54, 197], 
#[0.80, 0.92, 1.00, 1.00, 0.93, 1.02, 1.17, 1.26],
#[0.04, 0.03, 0.04, 0.00, 0.05, 0.06, 0.06, 0.10],
#marker='s', markersize=10, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

# per nucleus to C
plt.errorbar([9, 10, 11, 12, 40, 48, 54, 197], 
[0.53, 0.77, 0.83, 1.00, 3.09, 3.39, 5.06, 16.62],
[0.03, 0.03, 0.03, 0.00, 0.17, 0.20, 0.28, 1.29],
marker='s', markersize=15, mfc='k', ecolor='k', mec='k', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=4, label='')


# Meytal per nucleus to C ( conversion: Be: (per_proton X 4) /6, B10, B11: (per_proton X 5) /6, Ca40, 48: (per_proton X 20) /6, Fe54: (per_proton x26)/6, Au: (per_proton x 79/6)
#plt.errorbar(
#[12], 
#[1],  
#[0], 
# marker='o', markersize=15, mfc='None', ecolor='red', mec='red', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None', zorder=1, label='')

eps = 0.2  # offset to shift theory slightly from data for better visual


'''
# ----- PLOT MODELS (per  proton to C vs A) ----
plt.plot([9+eps, 10+eps, 11+eps], [0.935, 0.896, 1.012], marker='s', markersize=15, alpha=0.5, mfc='blue', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([9+eps, 10+eps, 11+eps], [1.025, 0.973, 1.0392], marker='D', markersize=15, alpha=0.5, mfc='blue', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([9+eps, 10+eps, 11+eps], [1.011, 0.956, 1.039], marker='P', markersize=15, alpha=0.5, mfc='blue', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)
#updated from Larry's spreadsheet 
plt.plot([9+eps, 10+eps, 11+eps], [0.899, 0.896, 0.9899], marker='^', markersize=15, alpha=0.5, mfc='green', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)
plt.plot([9+eps, 10+eps, 11+eps], [1.0, 0.8, 0.9995], marker='v', markersize=15, alpha=0.5, mfc='green', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)
plt.plot([9+eps, 11+eps], [0.9105, 1.02], marker='o', markersize=15, alpha=0.5, mfc='red', mec='black', linestyle='None', label='Spatial', zorder=3)
#------------------------------------------------
'''



# ----- PLOT MODELS (per  nucleus to C vs A) ----
plt.plot([9+eps, 10+eps, 11+eps], [0.623, 0.747, 0.843], marker='s', markersize=12, alpha=0.5, mfc='blue', mec='black', linestyle='None', label='Momentum (AV18)', zorder=1)#Av18
plt.plot([9+eps, 10+eps, 11+eps], [0.683, 0.811, 0.866], marker='D', markersize=12, alpha=0.5, mfc='blue', mec='black', linestyle='None', label='Momentum (SRG)', zorder=1)
plt.plot([9+eps, 10+eps, 11+eps], [0.674, 0.797, 0.866], marker='P', markersize=12, alpha=0.5, mfc='blue', mec='black', linestyle='None', label='Momentum (LCA)', zorder=3)

#updated from Larry's spreadsheet 
plt.plot([9+eps, 10+eps, 11+eps], [0.599, 0.747, 0.825], marker='^', markersize=12, alpha=0.5, mfc='green', mec='black', linestyle='None', label='l=0, n=0 pairs', zorder=3)

plt.plot([9+eps, 10+eps, 11+eps], [0.667, 0.667, 0.833], marker='v', markersize=12, alpha=0.5, mfc='green', mec='black', linestyle='None', label='l=0, L=0 pairs', zorder=3)

#Spatial Overlap Models

#--Spatial
#plt.plot([9+eps, 11+eps], [0.607, 0.850], marker='o', markersize=12, alpha=0.5, mfc='red', mec='black', linestyle='None', label='Spatial', zorder=3)

#--WS+CK
plt.plot([9+eps, 10+eps, 11+eps], [0.557, 0.725, 0.796], marker='o', markersize=12, alpha=0.5, mfc='red', mec='black', linestyle='None', label='WS+CK', zorder=3)

#Wiringa
plt.plot([9+eps, 10+eps, 11+eps], [0.56, 0.63, 0.81], marker='>', markersize=12, alpha=0.5, mfc='orange', mec='black', linestyle='None', label='Wiringa', zorder=3)


#-----------------------------------------------



# ---- per nucleus -----
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
plt.xlim(8.5, 12.5)
plt.ylim(0.4, 1.1)  


plt.title('', fontsize=25)
plt.xlabel('A', fontsize=25)

plt.yticks([0.5, 0.7, 0.9, 1.1])  

plt.xticks([9, 10, 11, 12])  

plt.tight_layout()
plt.legend(labelspacing=1.5, loc='lower right', fontsize=12)
plt.savefig('plb_fig2_pn.pdf')
plt.close()
#-------------------------


'''
# ----- per proton --------
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
plt.xlim(8.5, 12.5)
plt.ylim(0.4, 1.1)  


plt.title('', fontsize=25)
plt.xlabel('A', fontsize=25)

plt.yticks([0.5, 0.7, 0.9, 1.1])  

plt.xticks([9, 10, 11, 12])  

plt.tight_layout()
plt.legend(labelspacing=2.0, loc='lower right')
plt.savefig('prl_fig2_pp.pdf')
plt.close()
'''



#--------------------------









