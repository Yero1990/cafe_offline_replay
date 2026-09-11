

#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# This file when run asks for no inputs from the user but should be in the same folder as the appropriately formated .csv files from the ratiohistos.C code
# this file produces Q2, thrq, xbj, pmiss, emiss and pmiss ratios.
# Do note bools on lines 43 and 44. False removes the 2 sigma dashed lines showing the cut variation. Simple turn it to true and run again to produce plots with the lines
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
import math
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

#cutvar = False
cutvar = True



'''
ifname= 'SRCC12Q2.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197Q2.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40Q2.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48Q2.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54Q2.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax101.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=1.8, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=1.7, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=1.9, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=1.8, color='black', linestyle='--', linewidth=3)

plt.fill_between([1.8,3.0], [0,0], [2200,2200], alpha=0.075, edgecolor='black', facecolor='black')

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
ax101.set_xlabel(r'Q$^{2}$ [GeV/c]$^{2}$', fontsize=25)
ax101.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

#sm = mpl.cm.ScalarMappable(cmap="YlOrBr", norm=mpl.colors.Normalize(30, 60))
#cbar = plt.colorbar(sm, ax=ax101, ticks=[], label="A")  
#cbar.set_label('A', size=22)

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax101.set_xlim(1.5, 3.0)
ax101.set_ylim(0.0, 850.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('SRC_Q2_CutVar.png')
else:
    fig.savefig('SRC_Q2_Cut.png')
plt.close()





ifname= 'SRCBe9Q2.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10Q2.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11Q2.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12Q2.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax201 = plt.subplots(figsize=(10,7))
ax201.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax201.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax201.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax201.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=1.8, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=1.7, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=1.9, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=1.8, color='black', linestyle='--', linewidth=3)

plt.fill_between([1.8,3.0], [0,0], [1200,1200], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax201.set_title('', fontsize=21)
ax201.set_xlabel(r'Q$^{2}$ [GeV/c]$^{2}$', fontsize=25)
ax201.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax201.set_xlim(1.5, 3.0)
ax201.set_ylim(0.0, 45.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 10.0, 20.0, 30.0, 40.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LSRC_Q2_CutVar.png')
else:
    fig.savefig('LSRC_Q2_Cut.png')
plt.close()














ifname= 'MFC12Q2.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197Q2.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40Q2.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48Q2.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54Q2.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax107 = plt.subplots(figsize=(10,7))
ax107.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax107.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax107.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax107.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax107.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=1.8, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=1.7, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=1.9, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=1.8, color='black', linestyle='--', linewidth=3)

plt.fill_between([1.8,3.0], [0,0], [45000,45000], alpha=0.075, edgecolor='black', facecolor='black')

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
ax107.set_xlabel(r'Q$^{2}$ [GeV/c]$^{2}$', fontsize=25)
ax107.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax107.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax107.set_xlim(1.0, 3.0)
ax107.set_ylim(0.0, 45000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 10000.0, 20000.0, 30000.0, 40000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('MF_Q2_CutVar.png')
else:
    fig.savefig('MF_Q2_Cut.png')
plt.close()







ifname= 'MFBe9Q2.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10Q2.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11Q2.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12Q2.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax207 = plt.subplots(figsize=(10,7))
ax207.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax207.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax207.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax207.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=1.8, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=1.9, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=1.7, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=1.8, color='black', linestyle='--', linewidth=3)

plt.fill_between([1.8,3.0], [0,0], [45000,45000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax207.set_title('', fontsize=21)
ax207.set_xlabel(r'Q$^{2}$ [GeV/c]$^{2}$', fontsize=25)
ax207.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax207.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax207.set_xlim(1.0, 3.0)
ax207.set_ylim(0.0, 12000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LMF_Q2_CutVar.png')
else:
    fig.savefig('LMF_Q2_Cut.png')
plt.close()















ifname= 'SRCC12thrq.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197thrq.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40thrq.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48thrq.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54thrq.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax102 = plt.subplots(figsize=(10,7))
ax102.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax102.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax102.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax102.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax102.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=40, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=36, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=44, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=40, color='black', linestyle='--', linewidth=3)

plt.fill_between([0.0,40.0], [0,0], [1800,1800], alpha=0.075, edgecolor='black', facecolor='black')

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

ax102.set_title('', fontsize=21)
ax102.set_xlabel(r'$\theta_{rq}$ [deg]', fontsize=25)
ax102.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax102.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0])
ax102.set_xlim(0.0, 70.0)
ax102.set_ylim(0.0, 600.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('SRC_thrq_CutVar.png')
else:
    fig.savefig('SRC_thrq_Cut.png')
plt.close()











ifname= 'SRCBe9thrq.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10thrq.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11thrq.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12thrq.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax202 = plt.subplots(figsize=(10,7))
ax202.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax202.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax202.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax202.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=40, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=36, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=44, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=40, color='black', linestyle='--', linewidth=3)

plt.fill_between([0.0,40.0], [0,0], [900,900], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax202.set_title('', fontsize=21)
ax202.set_xlabel(r'$\theta_{rq}$ [deg]', fontsize=25)
ax202.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax202.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0])
ax202.set_xlim(0.0, 70.0)
ax202.set_ylim(0.0, 35.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 10.0, 20.0, 30.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LSRC_thrq_CutVar.png')
else:
    fig.savefig('LSRC_thrq_Cut.png')
plt.close()














ifname= 'SRCC12xbj.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197xbj.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40xbj.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48xbj.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54xbj.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax103 = plt.subplots(figsize=(10,7))
ax103.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax103.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax103.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax103.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax103.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=1.2, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=1.1, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=1.3, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=1.2, color='black', linestyle='--', linewidth=3)

plt.fill_between([1.2,2.0], [0,0], [4500,4500], alpha=0.075, edgecolor='black', facecolor='black')

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
ax103.set_xlabel(r'x$_{bj}$', fontsize=25, labelpad=-2)
ax103.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax103.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.5, 1.0, 1.5, 2.0])
ax103.set_xlim(0.5, 2.0)
ax103.set_ylim(0.0, 1800.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('SRC_xbj_CutVar.png')
else:
    fig.savefig('SRC_xbj_Cut.png')
plt.close()








ifname= 'SRCBe9xbj.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10xbj.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11xbj.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12xbj.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax203 = plt.subplots(figsize=(10,7))
ax203.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax203.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax203.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax203.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=1.2, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=1.1, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=1.3, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=1.2, color='black', linestyle='--', linewidth=3)

plt.fill_between([1.2,2.0], [0,0], [1800,1800], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax203.set_title('', fontsize=21)
ax203.set_xlabel(r'x$_{bj}$', fontsize=25, labelpad=-2)
ax203.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax203.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.5, 1.0, 1.5, 2.0])
ax203.set_xlim(0.5, 2.0)
ax203.set_ylim(0.0, 80.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 20.0, 40.0, 60.0, 80.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LSRC_xbj_CutVar.png')
else:
    fig.savefig('LSRC_xbj_Cut.png')
plt.close()










ifname= 'SRCC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197Pm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax104 = plt.subplots(figsize=(10,7))
ax104.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax104.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax104.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax104.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax104.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

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

plt.fill_between([0.375,0.7], [0,0], [2500,2500], alpha=0.075, edgecolor='black', facecolor='black')

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

ax104.set_title('', fontsize=21)
ax104.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax104.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax104.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
ax104.set_xlim(0.2, 0.83)
ax104.set_ylim(0.0, 650.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('SRC_pmiss_CutVar.png')
else:
    fig.savefig('SRC_pmiss_Cut.png')
plt.close()









ifname= 'SRCBe9Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10Pm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11Pm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax104 = plt.subplots(figsize=(10,7))
ax104.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax104.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax104.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax104.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=0.35, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.375, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.4, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.6, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.7, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.8, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=0.375, color='black', linestyle='--', linewidth=3)
    plt.axvline(x=0.7, color='black', linestyle='--', linewidth=3)

plt.fill_between([0.375,0.7], [0,0], [250,250], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax104.set_title('', fontsize=21)
ax104.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax104.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax104.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
ax104.set_xlim(0.2, 0.83)
ax104.set_ylim(0.0, 80.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 20.0, 40.0, 60.0, 80.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LSRC_pmiss_CutVar.png')
else:
    fig.savefig('LSRC_pmiss_Cut.png')
plt.close()








ifname= 'MFC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197Pm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax108 = plt.subplots(figsize=(10,7))
ax108.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax108.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax108.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax108.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax108.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=0.27, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.25, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.29, color='black', linestyle=':', linewidth=3)
else:
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

ax108.set_title('', fontsize=21)
ax108.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax108.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax108.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.3])
ax108.set_xlim(0.0, 0.32)
ax108.set_ylim(0.0, 12000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2500.0, 5000.0, 7500.0, 10000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('MF_pmiss_CutVar.png')
else:
    fig.savefig('MF_pmiss_Cut.png')
plt.close()






ifname= 'MFBe9Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10Pm.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11Pm.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax208 = plt.subplots(figsize=(10,7))
ax208.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax208.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax208.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax208.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=0.27, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.25, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.29, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=0.27, color='black', linestyle='--', linewidth=3)

plt.fill_between([0.0,0.27], [0,0], [14000,14000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax208.set_title('', fontsize=21)
ax208.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax208.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax208.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.3])
ax208.set_xlim(0.0, 0.32)
ax208.set_ylim(0.0, 4000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0, 4000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LMF_pmiss_CutVar.png')
else:
    fig.savefig('LMF_pmiss_Cut.png')
plt.close()












ifname= 'MFC12Em.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197Em.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40Em.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48Em.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54Em.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax109 = plt.subplots(figsize=(10,7))
ax109.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax109.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax109.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax109.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax109.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=-0.02, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.09, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.085, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.095, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=-0.02, color='black', linestyle='--', linewidth=3)
    plt.axvline(x=0.09, color='black', linestyle='--', linewidth=3)

plt.fill_between([-0.02,0.09], [0,0], [20000,20000], alpha=0.075, edgecolor='black', facecolor='black')

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

ax109.set_title('', fontsize=21)
ax109.set_xlabel(r'Missing Energy [GeV]', fontsize=25, labelpad=-2)
ax109.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax109.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.1])
ax109.set_xlim(-0.03, 0.13)
ax109.set_ylim(0.0, 20000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 5000.0, 10000.0, 15000.0, 20000.0, ])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('MF_Em_CutVar.png')
else:
    fig.savefig('MF_Em_Cut.png')
plt.close()










ifname= 'MFBe9Em.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10Em.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11Em.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12Em.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax209 = plt.subplots(figsize=(10,7))
ax209.plot(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax209.plot(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax209.plot(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax209.plot(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=-0.02, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.09, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.085, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.095, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=-0.02, color='black', linestyle='--', linewidth=3)
    plt.axvline(x=0.09, color='black', linestyle='--', linewidth=3)

plt.fill_between([-0.02,0.09], [0,0], [30000,30000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))


ax209.set_title('', fontsize=21)
ax209.set_xlabel(r'Missing Energy [GeV]', fontsize=25, labelpad=-2)
ax209.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')


ax209.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.1])
ax209.set_xlim(-0.03, 0.13)
ax209.set_ylim(0.0, 7000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LMF_Em_CutVar.png')
else:
    fig.savefig('LMF_Em_Cut.png')
plt.close()














ifname= 'MFBe9Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y01 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10Pm.csv'
df = pd.read_csv(ifname, comment='#')
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12Pm.csv'
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

for int in range(0,4):
    x = x[:-1]
    y01 = y01[:-1]
    e1 = e1[:-1]
    y02 = y02[:-1]
    e2 = e2[:-1]
    y03 = y03[:-1]
    e3 = e3[:-1]
    y04 = y04[:-1]
    e4 = e4[:-1]

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
vB10 = sum(nB10) / sum(dB10)
vB11 = sum(nB11) / sum(dB11)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yBe9, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{9}$Be', zorder=3)
ax2091.plot(x, yB10, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{10}$B', zorder=3)
ax2091.plot(x, yB11, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{11}$B', zorder=3)

plt.hlines(y=[vBe9, vB10, vB11], xmin=0.0, xmax=0.3, colors=['black', 'black', 'black'], linestyles=['--', '--', '--'], linewidth=[5,5,5])

#plt.text(0.7, 1.474, r'$\dfrac{^{9} \mathrm{Be} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')
#plt.text(0.7, 1.274, r'$\dfrac{^{10} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')
#plt.text(0.7, 1.074, r'$\dfrac{^{11} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yB11, linestyle='none')
plt.fill_between(x, yB11-eB11, yB11+eB11, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.3])
ax2091.set_xlim(0.0, 0.3)
ax2091.set_ylim(0.0, 1.5)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_pmiss_ratio.png')
plt.close()



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
ax2091.plot(x, yBe9, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{9}$Be', zorder=3)
ax2091.plot(x, yB10, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{10}$B', zorder=3)
ax2091.plot(x, yB11, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{11}$B', zorder=3)

plt.hlines(y=[vBe9, vB10, vB11], xmin=0.3, xmax=0.75, colors=['black', 'black', 'black'], linestyles=['--', '--', '--'], linewidth=[5,5,5])

#plt.text(0.7, 1.474, r'$\dfrac{^{9} \mathrm{Be} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')
#plt.text(0.7, 1.274, r'$\dfrac{^{10} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')
#plt.text(0.7, 1.074, r'$\dfrac{^{11} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='b')

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yB11, linestyle='none')
plt.fill_between(x, yB11-eB11, yB11+eB11, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('SRC Per Nucleus Ratio to C', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax2091.set_xlim(0.35, 0.75)
ax2091.set_ylim(0.4, 1.05)
plt.xticks(fontsize = 22)
plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('prl_fig1_pn_to_C12_vs_pmiss_ratio.png')
plt.close()

#---------END: PRL FIG 1, Per nucleus to C vs Pmiss--------------



'''


ifname= 'MFCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y01 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

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

for int in range(0,4):
    x = x[:-1]
    y01 = y01[:-1]
    e1 = e1[:-1]
    y02 = y02[:-1]
    e2 = e2[:-1]
    y03 = y03[:-1]
    e3 = e3[:-1]
    y04 = y04[:-1]
    e4 = e4[:-1]

yBe9 = np.divide(y02, y01)
yB10 = np.divide(y03, y02)
eBe9 = np.sqrt(1/e2**2 + 1/e1**2)
eB10 = np.sqrt(1/e3**2 + 1/e2**2)
eBe9 = np.sqrt( ( (1/y01) * e2 )**2 + ( (y02/y01**2) * e1)**2)
eB10 = np.sqrt( ( (1/y02) * e2 )**2 + ( (y03/y02**2) * e3)**2)
nBe9 = 1 / (eBe9**2) * yBe9
nB10 = 1 / (eB10**2) * yB10
dBe9 = 1 / (eBe9**2)
dB10 = 1 / (eB10**2)
vBe9 = sum(nBe9) / sum(dBe9)
vB10 = sum(nB10) / sum(dB10)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yBe9, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax2091.plot(x, yB10, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.hlines(y=[vBe9, vB10], xmin=0.0, xmax=0.3, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])

#plt.text(0.7, 1.474, r'$\dfrac{^{48} \mathrm{Ca} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color='b')
#plt.text(0.7, 1.274, r'$\dfrac{^{54} \mathrm{Fe} }{ ^{48} \mathrm{Ca} }$', fontsize=25, color='b')

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.3])
ax2091.set_xlim(0.0, 0.3)
ax2091.set_ylim(0.0, 2.2)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.5, 1.0, 1.5, 2.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_pmiss_ratio.png')
plt.close()








ifname= 'SRCCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y01 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

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

yBe9 = np.divide(y02, y01)
yB10 = np.divide(y03, y02)
eBe9 = np.sqrt( ( (1/y01) * e2 )**2 + ( (y02/y01**2) * e1)**2)
eB10 = np.sqrt( ( (1/y02) * e2 )**2 + ( (y03/y02**2) * e3)**2)
nBe9 = 1 / (eBe9**2) * yBe9
nB10 = 1 / (eB10**2) * yB10
dBe9 = 1 / (eBe9**2)
dB10 = 1 / (eB10**2)
vBe9 = sum(nBe9) / sum(dBe9)
vB10 = sum(nB10) / sum(dB10)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yBe9, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax2091.plot(x, yB10, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.hlines(y=[vBe9, vB10], xmin=0.3, xmax=0.8, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])

#plt.text(0.7, 1.474, r'$\dfrac{^{48} \mathrm{Ca} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color='b')
#plt.text(0.7, 1.274, r'$\dfrac{^{54} \mathrm{Fe} }{ ^{48} \mathrm{Ca} }$', fontsize=25, color='b')

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax2091.set_xlim(0.35, 0.75)
ax2091.set_ylim(0.95, 1.75)
plt.xticks(fontsize = 22)
plt.yticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_pmiss_ratio.png')
plt.close()









"""
ifname= 'eepD2.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'eepD2.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])
"""
fig, ax105 = plt.subplots(figsize=(10,7))
x=[0.903, 0.909, 0.915, 0.921, 0.927, 0.933, 0.939, 0.945, 0.951, 0.957, 0.963, 0.969, 0.975, 0.981, 0.987, 0.993, 0.999, 1.005, 1.011, 1.017, 1.023, 1.029, 1.035, 1.041, 1.047]
x2=[0.903, 0.909, 0.915, 0.921, 0.927, 0.933, 0.939, 0.945, 0.951, 0.957, 0.963, 0.969, 0.975, 0.981, 0.987, 0.993, 0.999, 1.005, 1.011, 1.017, 1.023, 1.029, 1.035, 1.041, 1.047]

ax105.plot(x, [137.586, 311.86, 858.77, 2546.5, 7462.63, 16292.8, 25956.6, 28975.3, 23386.2, 14215.4, 8123.01, 5089.41, 3748.02, 2969.56, 2536.19, 2219.76, 1876.94, 1662.53, 1464.16, 1043.36, 720.03, 401.29, 198.35, 54.85, 49.30],
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)

ax105.plot(x2, [2.4977, 20.79, 154.27, 1139.87, 4758.67, 13641.3, 25978.4, 32229.7, 27364.5, 17070.8, 8969.5, 5011.25, 2451.34, 2754.06, 2333.49, 2037.13, 1805.85, 1624.42, 1395.75, 1035.13, 778.61, 503.38, 240.69, 56.90, 7.13],
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)



y1=np.array([137.586, 311.86, 858.77, 2546.5, 7462.63, 16292.8, 25956.6, 28975.3, 23386.2, 14215.4, 8123.01, 5089.41, 3748.02, 2969.56, 2536.19, 2219.76, 1876.94, 1662.53, 1464.16, 1043.36, 720.03, 401.29, 198.35, 54.85, 49.30])
error1 = np.array([12.56, 18.91, 31.34, 54.03, 92.50, 136.68, 172.52, 183.27, 163.75, 127.67, 96.51, 76.39, 65.55, 58.35, 53.92, 50.45, 46.39, 43.66, 40.97, 34.59, 28.73, 21.45, 15.08, 9.86, 7.52])
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')



y2=np.array([2.4977, 20.79, 154.27, 1139.87, 4758.67, 13641.3, 25978.4, 32229.7, 27364.5, 17070.8, 8969.5, 5011.25, 2451.34, 2754.06, 2333.49, 2037.13, 1805.85, 1624.42, 1395.75, 1035.13, 778.61, 503.38, 240.69, 56.90, 7.13])
error2 = np.array([0.76, 220, 6.67, 16.54, 33.63, 56.82, 77.93, 87.00, 80.09, 63.08, 45.79, 34.20, 28.44, 25.45, 23.44, 21.93, 20.58, 19.58, 18.15, 15.62, 13.51, 10.73, 7.41, 3.58, 1.21])
plt.plot(x2, y2, linestyle='none')
plt.fill_between(x2, y2-error2, y2+error2, alpha=0.2, edgecolor='b', facecolor='b')




plt.text(0.955, 22000, r'Data', fontsize=25, color='r')
plt.text(0.95, 32000, r'Simulation', fontsize=25, color='b')

ax105.set_title('', fontsize=21)
ax105.set_xlabel(r'Invariant Mass [GeV]', fontsize=25)
ax105.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')



ax105.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.9, 0.95, 1.0, 1.05,])
ax105.set_xlim(0.9, 1.05)
plt.xticks(fontsize = 22)
plt.yticks([0, 10000, 20000, 30000])
ax105.set_ylim(0.0, 35000.0)
plt.yticks(fontsize = 22)

#plt.legend(frameon=False, fontsize=22, loc='upper left')
fig.set_tight_layout(True)

fig.savefig('Woverlay.png')
plt.close()




















fig, ax106 = plt.subplots(figsize=(10,7))
x=[0.903, 0.909, 0.915, 0.921, 0.927, 0.933, 0.939, 0.945, 0.951, 0.957, 0.963, 0.969, 0.975, 0.981, 0.987, 0.993, 0.999, 1.005, 1.011, 1.017, 1.023, 1.029, 1.035, 1.041, 1.047]
x2=[0.903, 0.909, 0.915, 0.921, 0.927, 0.933, 0.939, 0.945, 0.951, 0.957, 0.963, 0.969, 0.975, 0.981, 0.987, 0.993, 0.999, 1.005, 1.011, 1.017, 1.023, 1.029, 1.035, 1.041, 1.047]

ax106.plot(x, [5995.36, 10280.4, 20560.7, 50767.3, 120661, 234053, 348898, 393338, 346354, 258966, 173469, 112394, 78665.9, 59227, 50813.3, 43024.9, 38234.1, 33774.4, 33112.3, 28615.8, 28027.3, 26896.3, 24946.9, 23953.8, 23549.2],
marker='s', markersize=15, alpha=1.0, mfc='r', mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)

ax106.plot(x2, [1221.07, 4278.12, 14568.6, 46338, 128159, 280271, 440234, 482482, 380954, 234302, 131376, 78197.8, 56082.3, 44541.8, 37672, 33039.2, 29525.8, 26941.9, 24518.2, 22975.1, 21660.7, 20497.9, 18978.1, 18842.1, 17319.8],
marker='o', markersize=15, alpha=1.0, mfc='b', mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)



y1=np.array([5995.36, 10280.4, 20560.7, 50767.3, 120661, 234053, 348898, 393338, 346354, 258966, 173469, 112394, 78665.9, 59227, 50813.3, 43024.9, 38234.1, 33774.4, 33112.3, 28615.8, 28027.3, 26896.3, 24946.9, 23953.8, 23549.2])
error1 = np.array([234.80, 307.46, 434.82, 683.25, 1053.3, 1466.77, 1791.07, 1901.82, 1784.52, 1542.83, 1262.85, 1016.62, 850.51, 737.98, 683.56, 628.99, 592.94, 557.29, 551.8, 512.97, 507.67, 497.32, 478.96, 469.33, 465.35])
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-error1, y1+error1, alpha=0.2, edgecolor='r', facecolor='r')



y2=np.array([1221.07, 4278.12, 14568.6, 46338, 128159, 280271, 440234, 482482, 380954, 234302, 131376, 78197.8, 56082.3, 44541.8, 37672, 33039.2, 29525.8, 26941.9, 24518.2, 22975.1, 21660.7, 20497.9, 18978.1, 18842.1, 17319.8])
error2 = np.array([84.38, 151.31, 280.16, 517.19, 880.69, 1326.24, 1682.7, 1761.31, 1564.4, 1211.05, 905.78, 698.86, 597.16, 532.79, 494.52, 463.24, 435.59, 418.95, 398.50, 386.15, 378.22, 367.47, 353.16, 356.55, 340.65])
plt.plot(x2, y2, linestyle='none')
plt.fill_between(x2, y2-error2, y2+error2, alpha=0.2, edgecolor='b', facecolor='b')




plt.text(0.96, 260000, r'Data', fontsize=25, color='r')
plt.text(0.95, 475000, r'Simulation', fontsize=25, color='b')

ax106.set_title('', fontsize=21)
ax106.set_xlabel(r'Invariant Mass [GeV]', fontsize=25)
ax106.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')



ax106.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.9, 0.95, 1.0, 1.05,])
ax106.set_xlim(0.9, 1.05)
plt.xticks(fontsize = 22)
plt.yticks([0, 100000, 200000, 300000, 400000, 500000])
ax106.set_ylim(0.0, 500000.0)
plt.yticks(fontsize = 22)

#plt.legend(frameon=False, fontsize=22, loc='upper left')
fig.set_tight_layout(True)

fig.savefig('Woverlay_singles.png')
plt.close()







ifname= 'SRCC12Ef.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197Ef.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40Ef.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48Ef.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54Ef.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax101.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

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

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Final e$^{-}$ Momentum [GeV/c]', fontsize=25)
ax101.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(9.0, 10.0)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8, 10.0])
ax101.set_ylim(0.0, 600.0)
plt.yticks([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_Ef_Cut.png')
plt.close()



ifname= 'SRCAu197Ef.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40Ef.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48Ef.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc='red', mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.text(9.7, 65.0, r'$^{40}\mathrm{Ca}$', fontsize=25, color='red')
plt.text(9.7, 120.0, r'$^{48}\mathrm{Ca}$', fontsize=25, color='blue')
plt.text(9.82, 130.0, r'$^{54}\mathrm{Fe}$', fontsize=25, color='green')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='blue', facecolor='blue')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='green', facecolor='green')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Final e$^{-}$ Momentum [GeV/c]', fontsize=25)
ax101.set_ylabel('Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(9.0, 10.0)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8, 10.0])
ax101.set_ylim(0.0, 180.0)
plt.yticks([0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_Ef_PN.png')
plt.close()



ifname= 'SRCBe9Ef.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10Ef.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11Ef.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12Ef.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax201 = plt.subplots(figsize=(10,7))
ax201.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax201.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax201.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax201.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax201.set_title('', fontsize=21)
ax201.set_xlabel(r'Final e$^{-}$ Momentum [GeV/c]', fontsize=25)
ax201.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(9.0, 10.0)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8, 10.0])
ax201.set_ylim(0.0, 40.0)
plt.yticks([0.0, 10.0, 20.0, 30.0, 40.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_Ef_Cut.png')
plt.close()










ifname= 'SRCC12Pf.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197Pf.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40Pf.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48Pf.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54Pf.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax101.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

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

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Final Proton Momentum [GeV/c]', fontsize=25)
ax101.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(1.0, 1.5)
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
ax101.set_ylim(0.0, 600.0)
plt.yticks([0.0, 200.0, 400.0, 600.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)
fig.savefig('SRC_Pf_Cut.png')
plt.close()



ifname= 'SRCAu197Pf.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40Pf.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48Pf.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])



fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc='red', mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.text(1.24, 63.0, r'$^{40}\mathrm{Ca}$', fontsize=25, color='red')
plt.text(1.24, 97.0, r'$^{48}\mathrm{Ca}$', fontsize=25, color='blue')
plt.text(1.3, 120.0, r'$^{54}\mathrm{Fe}$', fontsize=25, color='green')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='blue', facecolor='blue')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='green', facecolor='green')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Final Proton Momentum [GeV/c]', fontsize=25)
ax101.set_ylabel('Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(1.0, 1.5)
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
ax101.set_ylim(0.0, 150.0)
plt.yticks([0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_Pf_PN.png')
plt.close()



ifname= 'SRCBe9Pf.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10Pf.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11Pf.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12Pf.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax201 = plt.subplots(figsize=(10,7))
ax201.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax201.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax201.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax201.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax201.set_title('', fontsize=21)
ax201.set_xlabel(r'Final Proton Momentum [GeV/c]', fontsize=25)
ax201.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(1.0, 1.5)
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
ax201.set_ylim(0.0, 30.0)
plt.yticks([0.0, 10.0, 20.0, 30.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_Pf_Cut.png')
plt.close()

















ifname= 'SRCC12the.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197the.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40the.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48the.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54the.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax101.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

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

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Electron Scattering Angle [deg]', fontsize=25)
ax101.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(7.0, 10.0)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
ax101.set_ylim(0.0, 270.0)
plt.yticks([0.0, 50.0, 100.0, 150.0, 200.0, 250])

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_the_Cut.png')
plt.close()



ifname= 'SRCAu197the.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40the.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48the.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc='red', mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.text(8.1, 28.0, r'$^{40}\mathrm{Ca}$', fontsize=25, color='red')
plt.text(8.1, 46.0, r'$^{48}\mathrm{Ca}$', fontsize=25, color='blue')
plt.text(8.5, 60.0, r'$^{54}\mathrm{Fe}$', fontsize=25, color='green')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='blue', facecolor='blue')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='green', facecolor='green')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Electron Scattering Angle [deg]', fontsize=25)
ax101.set_ylabel('Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(7.0, 10.0)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
ax101.set_ylim(0.0, 70.0)
plt.yticks([0.0, 20.0, 40.0, 60.0])

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_the_PN.png')
plt.close()



ifname= 'SRCBe9the.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10the.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11the.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12the.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax201 = plt.subplots(figsize=(10,7))
ax201.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax201.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax201.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax201.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax201.set_title('', fontsize=21)
ax201.set_xlabel(r'Electron Scattering Angle [deg]', fontsize=25)
ax201.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(7.0, 10.0)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
ax201.set_ylim(0.0, 15.0)
plt.yticks([0.0, 5.0, 10.0, 15.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_the_Cut.png')
plt.close()































ifname= 'SRCC12thp.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197thp.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40thp.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48thp.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54thp.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax101.errorbar(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap(rat01), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap(rat1), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap(rat2), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

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

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Proton Scattering Angle [deg]', fontsize=25)
ax101.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([62.0, 64.0, 66.0, 68.0, 70.0])
ax101.set_xlim(62.0, 70.0)
plt.yticks([0.0, 100.0, 200.0, 300.0, 400.0, 500.0])
ax101.set_ylim(0.0, 520.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_thp_Cut.png')
plt.close()



ifname= 'SRCCa40thp.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48thp.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54thp.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

fig, ax101 = plt.subplots(figsize=(10,7))
ax101.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc='red', mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.text(66, 62.0, r'$^{40}\mathrm{Ca}$', fontsize=25, color='red')
plt.text(66, 107.0, r'$^{48}\mathrm{Ca}$', fontsize=25, color='blue')
plt.text(66.5, 144.0, r'$^{54}\mathrm{Fe}$', fontsize=25, color='green')

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='blue', facecolor='blue')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='green', facecolor='green')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Proton Scattering Angle [deg]', fontsize=25)
ax101.set_ylabel('Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([62.0, 64.0, 66.0, 68.0, 70.0])
ax101.set_xlim(62.0, 70.0)
plt.yticks([0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0])
ax101.set_ylim(0.0, 175.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_thp_PN.png')
plt.close()



ifname= 'SRCBe9thp.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10thp.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11thp.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12thp.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

fig, ax201 = plt.subplots(figsize=(10,7))
ax201.errorbar(x, y1, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax201.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax201.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax201.errorbar(x, y4, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax201.set_title('', fontsize=21)
ax201.set_xlabel(r'Proton Scattering Angle [deg]', fontsize=25)
ax201.set_ylabel('Per Nucleus Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([62.0, 64.0, 66.0, 68.0, 70.0])
ax201.set_xlim(62.0, 70.0)
plt.yticks([0.0, 10.0, 20.0, 30.0])
ax201.set_ylim(0.0, 30.0)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_thp_Cut.png')
plt.close()

'''
