
#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# This file when run asks for no inputs from the user but should be in the same folder as the appropriately formated .csv files from the ratiohistos.C code
# this file produces edelta, hdelta, epctime, pCalEtotTrkNorm, and W histograms
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










cmap = mpl.colormaps.get_cmap('Blues')
cmap2 = mpl.colormaps.get_cmap('YlOrBr')

rat0 = 0.3
rat01 = 0.45
rat1 = 0.6
rat2 = 0.8
rat3 = 1.0

rat4 = 0.4
rat5 = 0.6
rat6 = 0.8
rat7 = 1.0





ifname= 'SRCC12edelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197edelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40edelta.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48edelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54edelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax101.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=22.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.0,22.0], [0,0], [2200,2200], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

#plt.text(2.05, 100, r'$^{40}$Ca', fontsize=25, color='r')
#plt.text(2.05, 153, r'$^{48}$Ca', fontsize=25, color='b')
#plt.text(1.82, 190, r'$^{54}$Fe', fontsize=25, color='g')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'$\delta_{e}$ [%]', fontsize=25)
ax101.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

#sm = mpl.cm.ScalarMappable(cmap="YlOrBr", norm=mpl.colors.Normalize(30, 60))
#cbar = plt.colorbar(sm, ax=ax101, ticks=[], label="A")  
#cbar.set_label('A', size=22)

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([10.0, 12.0, 14.0, 16.0, 18.0])
ax101.set_xlim(10.0, 18.0)
ax101.set_ylim(0.0, 1000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_edelta_Cut.png')
plt.close()




ifname= 'SRCBe9edelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10edelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11edelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12edelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax102 = plt.subplots(figsize=(10,7))
ax102.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax102.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax102.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax102.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=22.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.0,22.0], [0,0], [2200,2200], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax102.set_title('', fontsize=21)
ax102.set_xlabel(r'$\delta_{e}$ [%]', fontsize=25)
ax102.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax102.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([10.0, 12.0, 14.0, 16.0, 18.0])
ax102.set_xlim(10.0, 18.0)
ax102.set_ylim(0.0, 60.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 20.0, 40.0, 60.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_edelta_Cut.png')
plt.close()













ifname= 'MFC12edelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197edelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40edelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48edelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54edelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax103 = plt.subplots(figsize=(10,7))
ax103.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax103.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax103.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax103.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax103.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=22.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.0,22.0], [0,0], [4000,4000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax103.set_title('', fontsize=21)
ax103.set_xlabel(r'$\delta_{e}$ [%]', fontsize=25)
ax103.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax103.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([6.0, 8.0, 10.0, 12.0, 14.0])
ax103.set_xlim(6.0, 14.0)
ax103.set_ylim(0.0, 3200.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_edelta_Cut.png')
plt.close()






ifname= 'MFBe9edelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10edelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11edelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12edelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax104 = plt.subplots(figsize=(10,7))
ax104.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax104.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax104.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax104.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=22.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.0,22.0], [0,0], [22000,22000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax104.set_title('', fontsize=21)
ax104.set_xlabel(r'$\delta_{e}$ [%]', fontsize=25)
ax104.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax104.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([6.0, 8.0, 10.0, 12.0, 14.0])
ax104.set_xlim(6.0, 14.0)
ax104.set_ylim(0.0, 6500.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2000.0, 4000.0, 6000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_edelta_Cut.png')
plt.close()















ifname= 'SRCC12hdelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40hdelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax105 = plt.subplots(figsize=(10,7))
ax105.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax105.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax105.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax105.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax105.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-10.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=10.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-10.0,10.0], [0,0], [1800,1800], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax105.set_title('', fontsize=21)
ax105.set_xlabel(r'$\delta_{h}$ [%]', fontsize=25)
ax105.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax105.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0])
ax105.set_xlim(-15.0, 15.0)
ax105.set_ylim(0.0, 450.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_hdelta_Cut.png')
plt.close()










ifname= 'SRCBe9hdelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax106 = plt.subplots(figsize=(10,7))
ax106.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax106.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax106.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax106.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-10.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=10.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-10.0,10.0], [0,0], [1800,1800], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax106.set_title('', fontsize=21)
ax106.set_xlabel(r'$\delta_{h}$ [%]', fontsize=25)
ax106.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax106.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0])
ax106.set_xlim(-15.0, 15.0)
ax106.set_ylim(0.0, 23.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 5.0, 10.0, 15.0, 20.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_hdelta_Cut.png')
plt.close()











ifname= 'MFC12hdelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40hdelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax107 = plt.subplots(figsize=(10,7))
ax107.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax107.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax107.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax107.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax107.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-10.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=10.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-10.0,10.0], [0,0], [1800,1800], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax107.set_title('', fontsize=21)
ax107.set_xlabel(r'$\delta_{h}$ [%]', fontsize=25)
ax107.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax107.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0])
ax107.set_xlim(-15.0, 20.0)
ax107.set_ylim(0.0, 650.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_hdelta_Cut.png')
plt.close()










ifname= 'MFBe9hdelta.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12hdelta.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax108 = plt.subplots(figsize=(10,7))
ax108.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax108.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax108.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax108.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-10.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=10.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-10.0,10.0], [0,0], [1800,1800], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax108.set_title('', fontsize=21)
ax108.set_xlabel(r'$\delta_{h}$ [%]', fontsize=25)
ax108.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax108.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0])
ax108.set_xlim(-15.0, 20.0)
ax108.set_ylim(0.0, 1300.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_hdelta_Cut.png')
plt.close()








ifname= 'SRCC12epctime.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197epctime.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40epctime.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48epctime.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54epctime.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax109 = plt.subplots(figsize=(10,7))
#ax109.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
#ax109.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax109.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
#ax109.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
#ax109.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-2.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=2.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-2.0,2.0], [0,0], [4500,4500], alpha=0.075, edgecolor='black', facecolor='black')

#plt.plot(x, y0, linestyle='none')
#plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
#plt.plot(x, y1, linestyle='none')
#plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
#plt.plot(x, y3, linestyle='none')
#plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
#plt.plot(x, y4, linestyle='none')
#plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax109.set_title('', fontsize=21)
ax109.set_xlabel(r'ep Coin Time [ns]', fontsize=25, labelpad=-2)
ax109.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax109.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-16.0, -12.0, -8.0, -4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
ax109.set_xlim(-16, 28.0)
ax109.set_ylim(0.0, 200.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 50.0, 100.0, 150.0, 200.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_epctime_Cut.png')
plt.close()







ifname= 'SRCBe9epctime.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10epctime.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11epctime.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12epctime.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax110 = plt.subplots(figsize=(10,7))
#ax110.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
#ax110.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
#ax110.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax110.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-2.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=2.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-2.0,2.0], [0,0], [4500,4500], alpha=0.075, edgecolor='black', facecolor='black')

#plt.plot(x, y1, linestyle='none')
#plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
#plt.plot(x, y2, linestyle='none')
#plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
#plt.plot(x, y3, linestyle='none')
#plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax110.set_title('', fontsize=21)
ax110.set_xlabel(r'ep Coin Time [ns]', fontsize=25, labelpad=-2)
ax110.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax110.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
ax110.set_xlim(-4, 28.0)
ax110.set_ylim(0.0, 60.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_epctime_Cut.png')
plt.close()









ifname= 'MFC12epctime.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197epctime.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40epctime.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48epctime.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54epctime.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax111 = plt.subplots(figsize=(10,7))
ax111.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax111.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax111.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax111.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax111.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-2.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=2.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-2.0,2.0], [0,0], [4500,4500], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax111.set_title('', fontsize=21)
ax111.set_xlabel(r'ep Coin Time [ns]', fontsize=25, labelpad=-2)
ax111.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax111.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
ax111.set_xlim(-4, 28.0)
ax111.set_ylim(0.0, 4200.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0, 4000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_epctime_Cut.png')
plt.close()







ifname= 'MFBe9epctime.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10epctime.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11epctime.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12epctime.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax112 = plt.subplots(figsize=(10,7))
ax112.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax112.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax112.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax112.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-2.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=2.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-2.0,2.0], [0,0], [45000,45000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax112.set_title('', fontsize=21)
ax112.set_xlabel(r'ep Coin Time [ns]', fontsize=25, labelpad=-2)
ax112.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax112.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
ax112.set_xlim(-4, 28.0)
ax112.set_ylim(0.0, 8500.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2000.0, 4000.0, 6000.0, 8000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_epctime_Cut.png')
plt.close()

































ifname= 'SRCC12pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax113 = plt.subplots(figsize=(10,7))
#ax113.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
#ax113.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
#ax113.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
#ax113.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax113.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.8, color='black', linestyle='--', linewidth=3)
plt.axvline(x=1.3, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.8,1.3], [0,0], [2500,2500], alpha=0.075, edgecolor='black', facecolor='black')

#plt.plot(x, y0, linestyle='none')
#plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
#plt.plot(x, y1, linestyle='none')
#plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
#plt.plot(x, y2, linestyle='none')
#plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
#plt.plot(x, y3, linestyle='none')
#plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax113.set_title('', fontsize=21)
ax113.set_xlabel(r'E$_{tot}$/P$_{Trk}$', fontsize=25, labelpad=-2)
ax113.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax113.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax113.set_xlim(0.0, 2.0)
ax113.set_ylim(0.0, 1000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_pCalEtotTrkNorm_Cut.png')
plt.close()









ifname= 'SRCBe9pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax114 = plt.subplots(figsize=(10,7))
#ax114.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
#ax114.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
#ax114.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax114.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.8, color='black', linestyle='--', linewidth=3)
plt.axvline(x=1.3, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.8,1.3], [0,0], [2500,2500], alpha=0.075, edgecolor='black', facecolor='black')

#plt.plot(x, y1, linestyle='none')
#plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
#plt.plot(x, y2, linestyle='none')
#plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
#plt.plot(x, y3, linestyle='none')
#plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax114.set_title('', fontsize=21)
ax114.set_xlabel(r'E$_{tot}$/P$_{Trk}$', fontsize=25, labelpad=-2)
ax114.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax114.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax114.set_xlim(0.0, 2.0)
ax114.set_ylim(0.0, 60.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_pCalEtotTrkNorm_Cut.png')
plt.close()







ifname= 'MFC12pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax115 = plt.subplots(figsize=(10,7))
ax115.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax115.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax115.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax115.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax115.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.8, color='black', linestyle='--', linewidth=3)
plt.axvline(x=1.3, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.8,1.3], [0,0], [25000,25000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax115.set_title('', fontsize=21)
ax115.set_xlabel(r'E$_{tot}$/P$_{Trk}$', fontsize=25, labelpad=-2)
ax115.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax115.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax115.set_xlim(0.0, 2.0)
ax115.set_ylim(0.0, 3500.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_pCalEtotTrkNorm_Cut.png')
plt.close()





ifname= 'MFBe9pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12pCalEtotTrkNorm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax116 = plt.subplots(figsize=(10,7))
ax116.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax116.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax116.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax116.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.8, color='black', linestyle='--', linewidth=3)
plt.axvline(x=1.3, color='black', linestyle='--', linewidth=3)
plt.fill_between([0.8,1.3], [0,0], [25000,25000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax116.set_title('', fontsize=21)
ax116.set_xlabel(r'E$_{tot}$/P$_{Trk}$', fontsize=25, labelpad=-2)
ax116.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax116.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax116.set_xlim(0.0, 2.0)
ax116.set_ylim(0.0, 7000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_pCalEtotTrkNorm_Cut.png')
plt.close()
























ifname= 'SRCC12W.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197W.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40W.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48W.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54W.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax117 = plt.subplots(figsize=(10,7))
ax117.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax117.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax117.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax117.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax117.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

"""
if cutvar:
    plt.axvline(x=0.375, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.35, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.4, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.7, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.6, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.8, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=0.375, color='black', linestyle='--', linewidth=3)
    plt.axvline(x=0.7, color='black', linestyle='--', linewidth=3)
"""
plt.fill_between([0.0,2.0], [0,0], [25000,25000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax117.set_title('', fontsize=21)
ax117.set_xlabel(r'Invariant Mass [GeV]', fontsize=25, labelpad=-2)
ax117.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax117.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax117.set_xlim(0.0, 2.0)
ax117.set_ylim(0.0, 100.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 20.0, 40.0, 60.0, 80.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_W_Cut.png')
plt.close()









ifname= 'SRCBe9W.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10W.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11W.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12W.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax118 = plt.subplots(figsize=(10,7))
ax118.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax118.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax118.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax118.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

"""
if cutvar:
    plt.axvline(x=0.35, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.375, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.4, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.6, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.7, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.8, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=0.375, color='black', linestyle='--', linewidth=3)
    plt.axvline(x=0.7, color='black', linestyle='--', linewidth=3)
"""

plt.fill_between([0.0,2.0], [0,0], [2500,2500], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax118.set_title('', fontsize=21)
ax118.set_xlabel(r'Invariant Mass [GeV]', fontsize=25, labelpad=-2)
ax118.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax118.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax118.set_xlim(0.0, 2.0)
ax118.set_ylim(0.0, 50.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 10.0, 20.0, 30.0, 40.0, 50.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_W_Cut.png')
plt.close()







ifname= 'MFC12W.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197W.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40W.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48W.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54W.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax119 = plt.subplots(figsize=(10,7))
ax119.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax119.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax119.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax119.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax119.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=0.27, color='black', linestyle='--', linewidth=3)

plt.fill_between([0.0,0.27], [0,0], [12000,12000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax119.set_title('', fontsize=21)
ax119.set_xlabel(r'Invariant Mass [GeV]', fontsize=25, labelpad=-2)
ax119.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax119.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax119.set_xlim(0.0, 2.2)
ax119.set_ylim(0.0, 7000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0, 4000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_W_Cut.png')
plt.close()





ifname= 'MFBe9W.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10W.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11W.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12W.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax120 = plt.subplots(figsize=(10,7))
ax120.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax120.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax120.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax120.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

"""
if cutvar:
    plt.axvline(x=0.27, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.25, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.29, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=0.27, color='black', linestyle='--', linewidth=3)
"""
plt.fill_between([0.0,2.0], [0,0], [14000,14000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax120.set_title('', fontsize=21)
ax120.set_xlabel(r'Invariant Mass [GeV]', fontsize=25, labelpad=-2)
ax120.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax120.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax120.set_xlim(0.0, 2.0)
ax120.set_ylim(0.0, 10000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_W_Cut.png')
plt.close()











"""
ifname= 'SRCC12epctimea.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40epctimea.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax121 = plt.subplots(figsize=(10,7))
ax121.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax121.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax121.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax121.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax121.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-2.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=2.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-2.0,2.0], [0,0], [4500,4500], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax121.set_title('', fontsize=21)
ax121.set_xlabel(r'ep Coin Time [ns]', fontsize=25, labelpad=-2)
ax121.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax121.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-16.0, -12.0, -8.0, -4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
ax121.set_xlim(-16, 28.0)
ax121.set_ylim(0.0, 1000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_epctime_Cut_alt.png')
plt.close()















ifname= 'MFC12epctimea.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40epctimea.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax123 = plt.subplots(figsize=(10,7))
ax123.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax123.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax123.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax123.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax123.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-2.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=2.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-2.0,2.0], [0,0], [4500,4500], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap(rat01), facecolor=cmap(rat01))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap(rat1), facecolor=cmap(rat1))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap(rat2), facecolor=cmap(rat2))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))

ax123.set_title('', fontsize=21)
ax123.set_xlabel(r'ep Coin Time [ns]', fontsize=25, labelpad=-2)
ax123.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax123.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
ax123.set_xlim(-4, 28.0)
ax123.set_ylim(0.0, 4200.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0, 4000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_epctime_Cut_alt.png')
plt.close()







ifname= 'MFBe9epctimea.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12epctimea.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax124 = plt.subplots(figsize=(10,7))
ax124.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax124.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax124.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax124.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.axvline(x=-2.0, color='black', linestyle='--', linewidth=3)
plt.axvline(x=2.0, color='black', linestyle='--', linewidth=3)
plt.fill_between([-2.0,2.0], [0,0], [45000,45000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax124.set_title('', fontsize=21)
ax124.set_xlabel(r'ep Coin Time [ns]', fontsize=25, labelpad=-2)
ax124.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax124.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0])
ax124.set_xlim(-4, 28.0)
ax124.set_ylim(0.0, 8500.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2000.0, 4000.0, 6000.0, 8000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_epctime_Cut_alt.png')
plt.close()
"""







