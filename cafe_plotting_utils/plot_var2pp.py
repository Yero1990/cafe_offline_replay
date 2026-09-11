
#*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=
# This file when run asks for no inputs from the user but should be in the same folder as the appropriately formated .csv files from the ratiohistos.C code
# this file produces Q2, thrq, xbj, pmiss, emiss, final proton momentum, final electron momentum, theta_e, theta_p, nu, q, 
# pmiss ratio, MF and SRC single ratio without theory, and double ratio without theory plots alongside variations for our papers.
# Do note bools on lines 41 and 42. False removes the 2 sigma dashed lines showing the cut variation. Simple turn it to true and run again to produce plots with the lines
# This file has axis adjusted for per-proton normalization
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

cutvar = False
#cutvar = True




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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

#sm = mpl.cm.ScalarMappable(cmap="YlOrBr", norm=mpl.colors.Normalize(30, 60))
#cbar = plt.colorbar(sm, ax=ax101, ticks=[], label="A")  
#cbar.set_label('A', size=22)

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax101.set_xlim(1.5, 3.0)
ax101.set_ylim(0.0, 11.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax201.set_xlim(1.5, 3.0)
ax201.set_ylim(0.0, 8.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0])
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
ax107.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax107.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax107.set_xlim(1.0, 3.0)
ax107.set_ylim(0.0, 2000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 500.0, 1000.0, 1500.0, 2000.0])
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
ax207.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax207.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.0, 1.5, 2.0, 2.5, 3.0])
ax207.set_xlim(1.0, 3.0)
ax207.set_ylim(0.0, 2000.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 500.0, 1000.0, 1500.0, 2000.0])
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
ax102.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax102.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0])
ax102.set_xlim(0.0, 70.0)
ax102.set_ylim(0.0, 8.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0])
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
ax202.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax202.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0])
ax202.set_xlim(0.0, 70.0)
ax202.set_ylim(0.0, 6.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0])
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
ax103.set_xlabel(r'x$_{B}$', fontsize=25, labelpad=-2)
ax103.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax103.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.5, 1.0, 1.5, 2.0])
ax103.set_xlim(0.5, 2.0)
ax103.set_ylim(0.0, 22.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 5.0, 10.0, 15.0, 20.0])
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
ax203.set_xlabel(r'x$_{B}$', fontsize=25, labelpad=-2)
ax203.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax203.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.5, 1.0, 1.5, 2.0])
ax203.set_xlim(0.5, 2.0)
ax203.set_ylim(0.0, 12.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0])
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
ax104.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax104.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
ax104.set_xlim(0.2, 0.83)
ax104.set_ylim(0.0, 12.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0])
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
ax104.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax104.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
ax104.set_xlim(0.2, 0.83)
ax104.set_ylim(0.0, 12.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0])
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
ax108.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax108.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax108.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax108.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax108.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

if cutvar:
    plt.axvline(x=0.27, color='black', linestyle='-', linewidth=3)
    plt.axvline(x=0.25, color='black', linestyle=':', linewidth=3)
    plt.axvline(x=0.29, color='black', linestyle=':', linewidth=3)
else:
    plt.axvline(x=0.27, color='black', linestyle='--', linewidth=3)

plt.fill_between([0.0,0.27], [0,0], [12000,12000], alpha=0.075, edgecolor='black', facecolor='black')

plt.plot(x, y0, linestyle='none')
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax108.set_title('', fontsize=21)
ax108.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax108.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax108.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.3])
ax108.set_xlim(0.0, 0.32)
ax108.set_ylim(0.0, 600.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0])
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
ax208.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax208.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.3])
ax208.set_xlim(0.0, 0.32)
ax208.set_ylim(0.0, 650.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0])
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
ax109.plot(x, y0, marker='s', markersize=15, alpha=1.0, mfc=cmap(rat3), mec='None', linestyle='None', label=r'$^{40}$Ca', zorder=3)
ax109.plot(x, y1, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax109.plot(x, y2, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax109.errorbar(x, y3, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax109.errorbar(x, y4, marker='P', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

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
plt.fill_between(x, y0-e0, y0+e0, alpha=0.5, edgecolor=cmap(rat3), facecolor=cmap(rat3))
plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, y4, linestyle='none')
plt.fill_between(x, y4-e4, y4+e4, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax109.set_title('', fontsize=21)
ax109.set_xlabel(r'Missing Energy [GeV]', fontsize=25, labelpad=-2)
ax109.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax109.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.1])
ax109.set_xlim(-0.03, 0.13)
ax109.set_ylim(0.0, 1200.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0])
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
ax209.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')


ax209.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([-0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.1])
ax209.set_xlim(-0.03, 0.13)
ax209.set_ylim(0.0, 1300.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LMF_Em_CutVar.png')
else:
    fig.savefig('LMF_Em_Cut.png')
plt.close()





















"""
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

#plt.hlines(y=[vBe9, vB10, vB11], xmin=0.0, xmax=0.275, colors=['black', 'black', 'black'], linestyles=['--', '--', '--'], linewidth=[5,5,5])

plt.text(0.28, 0.73, r'$\dfrac{^{9} \mathrm{Be} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat4))
plt.text(0.28, 1.0, r'$\dfrac{^{10} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat5))
plt.text(0.28, 1.27, r'$\dfrac{^{11} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat6))

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yB11, linestyle='none')
plt.fill_between(x, yB11-eB11, yB11+eB11, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('MF Per Proton Ratio to C', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.3])
ax2091.set_xlim(0.0, 0.31)
ax2091.set_ylim(0.0, 2.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_pmiss_ratio.png')
plt.close()
"""







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
vB10 = sum(nB10) / sum(dB10)
vB11 = sum(nB11) / sum(dB11)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yBe9, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{9}$Be', zorder=3)
ax2091.plot(x, yB10, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{10}$B', zorder=3)
ax2091.plot(x, yB11, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{11}$B', zorder=3)

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yB11, linestyle='none')
plt.fill_between(x, yB11-eB11, yB11+eB11, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))

#plt.hlines(y=[vBe9, vB10, vB11], xmin=0.3, xmax=0.715, colors=['black', 'black', 'black'], linestyles=['--', '--', '--'], linewidth=[5,5,5])

plt.text(0.74, 0.7, r'$\dfrac{^{9} \mathrm{Be} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat4))
plt.text(0.74, 0.88, r'$\dfrac{^{10} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat5))
plt.text(0.74, 1.07, r'$\dfrac{^{11} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat6))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('SRC Per Proton Ratio to C', fontsize=25)
#ax2091.set_ylabel('Per Proton Ratio to C', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7, 0.8])
ax2091.set_xlim(0.35, 0.78)
ax2091.set_ylim(0.0, 1.3)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_pmiss_ratio.png')
plt.close()












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
vB10 = sum(nB10) / sum(dB10)
vB11 = sum(nB11) / sum(dB11)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yBe9, marker='s', markersize=15, alpha=1.0, mfc='red', mec='None', linestyle='None', label=r'$^{9}$Be', zorder=3)
ax2091.plot(x, yB10, marker='o', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{10}$B', zorder=3)
ax2091.plot(x, yB11, marker='^', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{11}$B', zorder=3)

plt.plot(x, yBe9, linestyle='none')
plt.fill_between(x, yBe9-eBe9, yBe9+eBe9, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, yB10, linestyle='none')
plt.fill_between(x, yB10-eB10, yB10+eB10, alpha=0.5, edgecolor='green', facecolor='green')
plt.plot(x, yB11, linestyle='none')
plt.fill_between(x, yB11-eB11, yB11+eB11, alpha=0.5, edgecolor='blue', facecolor='blue')

#plt.hlines(y=[vBe9, vB10, vB11], xmin=0.3, xmax=0.715, colors=['black', 'black', 'black'], linestyles=['--', '--', '--'], linewidth=[5,5,5])

plt.text(0.74, 0.7, r'$\dfrac{^{9} \mathrm{Be} }{ ^{12} \mathrm{C} }$', fontsize=25, color='red')
plt.text(0.74, 0.88, r'$\dfrac{^{10} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='green')
plt.text(0.74, 1.07, r'$\dfrac{^{11} \mathrm{B} }{ ^{12} \mathrm{C} }$', fontsize=25, color='blue')

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('SRC Per Proton Ratio to C', fontsize=25)
#ax2091.set_ylabel('Per Proton Ratio to C', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7, 0.8])
ax2091.set_xlim(0.35, 0.78)
ax2091.set_ylim(0.0, 1.3)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_pmiss_ratio_alt.png')
plt.close()





ifname= 'MFC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y01 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y04 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFAu197Pm.csv'
df = pd.read_csv(ifname, comment='#')
y05 = np.array(df['zcont'])
e5 = np.array(df['zcont_err'])

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
for index, value in enumerate(y05):
    if value == 0:
        y05[index] = 0.1
for index, value in enumerate(e5):
    if value == 0:
        e5[index] = 0.1

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
    y05 = y05[:-1]
    e5 = e5[:-1]

yCa40 = np.divide(y02, y01)
yCa48 = np.divide(y03, y01)
yFe54 = np.divide(y04, y01)
yAu197 = np.divide(y05, y01)
eCa40 = np.sqrt( ( (1/y01) * e2 )**2 + ( (y02/y01**2) * e1)**2)
eCa48 = np.sqrt( ( (1/y01) * e3 )**2 + ( (y03/y01**2) * e1)**2)
eFe54 = np.sqrt( ( (1/y01) * e4 )**2 + ( (y04/y01**2) * e1)**2)
eAu197 = np.sqrt( ( (1/y01) * e5 )**2 + ( (y05/y01**2) * e1)**2)
nCa40 = 1 / (eCa40**2) * yCa40
nCa48 = 1 / (eCa48**2) * yCa48
nFe54 = 1 / (eFe54**2) * yFe54
nAu197 = 1 / (eAu197**2) * yAu197
dCa40 = 1 / (eCa40**2)
dCa48 = 1 / (eCa48**2)
dFe54 = 1 / (eFe54**2)
dAu197 = 1 / (eAu197**2)
vCa40 = sum(nCa40) / sum(dCa40)
vCa48 = sum(nCa48) / sum(dCa48)
vFe54 = sum(nFe54) / sum(dFe54)
vAu197 = sum(nAu197) / sum(dAu197)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yCa40, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax2091.plot(x, yCa48, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax2091.plot(x, yFe54, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax2091.plot(x, yAu197, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

#plt.hlines(y=[vBe9, vB10], xmin=0.0, xmax=0.28, colors=['black', 'black'], linestyles=['--', '--'], linewidth=[5,5])

plt.text(0.282, 0.6, r'$\dfrac{^{40} \mathrm{Ca} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat4))
plt.text(0.282, 0.8, r'$\dfrac{^{48} \mathrm{Ca} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat5))
plt.text(0.282, 1.0, r'$\dfrac{^{54} \mathrm{Fe} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat6))
plt.text(0.282, 1.2, r'$\dfrac{^{197} \mathrm{Au} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat7))

plt.plot(x, yCa40, linestyle='none')
plt.fill_between(x, yCa40-eCa40, yCa40+eCa40, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yCa48, linestyle='none')
plt.fill_between(x, yCa48-eCa48, yCa48+eCa48, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yFe54, linestyle='none')
plt.fill_between(x, yFe54-eFe54, yFe54+eFe54, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))
plt.plot(x, yAu197, linestyle='none')
plt.fill_between(x, yAu197-eAu197, yAu197+eAu197, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat7))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('MF Per Proton Ratio to C', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.30])
ax2091.set_xlim(0.0, 0.32)
ax2091.set_ylim(0.0, 1.5)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_pmiss_ratio.png')
plt.close()








ifname= 'SRCC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y01 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y04 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCAu197Pm.csv'
df = pd.read_csv(ifname, comment='#')
y05 = np.array(df['zcont'])
e5 = np.array(df['zcont_err'])

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
for index, value in enumerate(y05):
    if value == 0:
        y05[index] = 0.1
for index, value in enumerate(e5):
    if value == 0:
        e5[index] = 0.1

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
    y05 = y05[:-1]
    e5 = e5[:-1]

yCa40 = np.divide(y02, y01)
yCa48 = np.divide(y03, y01)
yFe54 = np.divide(y04, y01)
yAu197 = np.divide(y05, y01)
eCa40 = np.sqrt( ( (1/y01) * e2 )**2 + ( (y02/y01**2) * e1)**2)
eCa48 = np.sqrt( ( (1/y01) * e3 )**2 + ( (y03/y01**2) * e1)**2)
eFe54 = np.sqrt( ( (1/y01) * e4 )**2 + ( (y04/y01**2) * e1)**2)
eAu197 = np.sqrt( ( (1/y01) * e5 )**2 + ( (y05/y01**2) * e1)**2)
nCa40 = 1 / (eCa40**2) * yCa40
nCa48 = 1 / (eCa48**2) * yCa48
nFe54 = 1 / (eFe54**2) * yFe54
nAu197 = 1 / (eAu197**2) * yAu197
dCa40 = 1 / (eCa40**2)
dCa48 = 1 / (eCa48**2)
dFe54 = 1 / (eFe54**2)
dAu197 = 1 / (eAu197**2)
vCa40 = sum(nCa40) / sum(dCa40)
vCa48 = sum(nCa48) / sum(dCa48)
vFe54 = sum(nFe54) / sum(dFe54)
vAu197 = sum(nAu197) / sum(dAu197)

fig, ax2091 = plt.subplots(figsize=(10,7))
ax2091.plot(x, yCa40, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{48}$Ca', zorder=3)
ax2091.plot(x, yCa48, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax2091.plot(x, yFe54, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax2091.plot(x, yAu197, marker='v', markersize=15, alpha=1.0, mfc=cmap2(rat7), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.hlines(y=[vCa40, vCa48, vFe54, vAu197], xmin=0.0, xmax=0.28, colors=['black', 'black', 'black', 'black'], linestyles=['--', '--', '--', '--'], linewidth=[5,5,5,5])

plt.text(0.7, 0.8, r'$\dfrac{^{40} \mathrm{Ca} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat4))
plt.text(0.7, 1.09, r'$\dfrac{^{48} \mathrm{Ca} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat5))
plt.text(0.7, 1.38, r'$\dfrac{^{54} \mathrm{Fe} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat6))
plt.text(0.7, 1.67, r'$\dfrac{^{197} \mathrm{Au} }{ ^{12} \mathrm{C} }$', fontsize=25, color=cmap2(rat7))

plt.plot(x, yCa40, linestyle='none')
plt.fill_between(x, yCa40-eCa40, yCa40+eCa40, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yCa48, linestyle='none')
plt.fill_between(x, yCa48-eCa48, yCa48+eCa48, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yFe54, linestyle='none')
plt.fill_between(x, yFe54-eFe54, yFe54+eFe54, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat4))
plt.plot(x, yAu197, linestyle='none')
plt.fill_between(x, yAu197-eAu197, yAu197+eAu197, alpha=0.5, edgecolor=cmap2(rat7), facecolor=cmap2(rat5))

ax2091.set_title('', fontsize=21)
ax2091.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax2091.set_ylabel('SRC Per Proton Ratio to C', fontsize=25)

ax2091.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax2091.set_xlim(0.35, 0.78)
ax2091.set_ylim(0.0, 2.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_pmiss_ratio.png')
plt.close()














ifname= 'MFCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y04 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFAu197Pm.csv'
df = pd.read_csv(ifname, comment='#')
y05 = np.array(df['zcont'])
e5 = np.array(df['zcont_err'])

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
for index, value in enumerate(y05):
    if value == 0:
        y05[index] = 0.1
for index, value in enumerate(e5):
    if value == 0:
        e5[index] = 0.1

for int in range(0,4):
    x = x[:-1]
    y02 = y02[:-1]
    e2 = e2[:-1]
    y03 = y03[:-1]
    e3 = e3[:-1]
    y04 = y04[:-1]
    e4 = e4[:-1]
    y05 = y05[:-1]
    e5 = e5[:-1]

yCa48 = np.divide(y03, y02)
yFe54 = np.divide(y04, y02)
yAu197 = np.divide(y05, y02)
eCa48 = np.sqrt( ( (1/y02) * e3 )**2 + ( (y03/y02**2) * e2)**2)
eFe54 = np.sqrt( ( (1/y02) * e4 )**2 + ( (y04/y02**2) * e2)**2)
eAu197 = np.sqrt( ( (1/y02) * e5 )**2 + ( (y05/y02**2) * e2)**2)
nCa48 = 1 / (eCa48**2) * yCa48
nFe54 = 1 / (eFe54**2) * yFe54
nAu197 = 1 / (eAu197**2) * yAu197
dCa48 = 1 / (eCa48**2)
dFe54 = 1 / (eFe54**2)
dAu197 = 1 / (eAu197**2)
vCa48 = sum(nCa48) / sum(dCa48)
vFe54 = sum(nFe54) / sum(dFe54)
vAu197 = sum(nAu197) / sum(dAu197)

fig, ax291 = plt.subplots(figsize=(10,7))
ax291.plot(x, yCa48, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax291.plot(x, yFe54, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax291.plot(x, yAu197, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

#plt.hlines(y=[vCa48, vFe54, vAu197], xmin=0.0, xmax=0.28, colors=['black', 'black', 'black'], linestyles=['--', '--', '--'], linewidth=[5,5,5])

plt.text(0.282, 0.8, r'$\dfrac{^{48} \mathrm{Ca} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color=cmap2(rat4))
plt.text(0.282, 1.0, r'$\dfrac{^{54} \mathrm{Fe} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color=cmap2(rat5))
plt.text(0.282, 1.2, r'$\dfrac{^{197} \mathrm{Au} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color=cmap2(rat6))

plt.plot(x, yCa48, linestyle='none')
plt.fill_between(x, yCa48-eCa48, yCa48+eCa48, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yFe54, linestyle='none')
plt.fill_between(x, yFe54-eFe54, yFe54+eFe54, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yAu197, linestyle='none')
plt.fill_between(x, yAu197-eAu197, yAu197+eAu197, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))

ax291.set_title('', fontsize=21)
ax291.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax291.set_ylabel(r'MF Per Proton Ratio to $^{40}$Ca', fontsize=25)

ax291.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.0, 0.1, 0.2, 0.30])
ax291.set_xlim(0.0, 0.32)
ax291.set_ylim(0.0, 1.5)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_pmiss_ratio_Ca40.png')
plt.close()








ifname= 'SRCC12Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y01 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa40Pm.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y02 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCCa48Pm.csv'
df = pd.read_csv(ifname, comment='#')
y03 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCFe54Pm.csv'
df = pd.read_csv(ifname, comment='#')
y04 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCAu197Pm.csv'
df = pd.read_csv(ifname, comment='#')
y05 = np.array(df['zcont'])
e5 = np.array(df['zcont_err'])

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
for index, value in enumerate(y05):
    if value == 0:
        y05[index] = 0.1
for index, value in enumerate(e5):
    if value == 0:
        e5[index] = 0.1

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
    y05 = y05[:-1]
    e5 = e5[:-1]



yCa48 = np.divide(y03, y02)
yFe54 = np.divide(y04, y02)
yAu197 = np.divide(y05, y02)
eCa48 = np.sqrt( ( (1/y02) * e3 )**2 + ( (y03/y02**2) * e2)**2)
eFe54 = np.sqrt( ( (1/y02) * e4 )**2 + ( (y04/y02**2) * e2)**2)
eAu197 = np.sqrt( ( (1/y02) * e5 )**2 + ( (y05/y02**2) * e2)**2)
nCa48 = 1 / (eCa48**2) * yCa48
nFe54 = 1 / (eFe54**2) * yFe54
nAu197 = 1 / (eAu197**2) * yAu197
dCa48 = 1 / (eCa48**2)
dFe54 = 1 / (eFe54**2)
dAu197 = 1 / (eAu197**2)
vCa48 = sum(nCa48) / sum(dCa48)
vFe54 = sum(nFe54) / sum(dFe54)
vAu197 = sum(nAu197) / sum(dAu197)

fig, ax292 = plt.subplots(figsize=(10,7))

ax292.plot(x, yCa48, marker='s', markersize=15, alpha=1.0, mfc=cmap2(rat4), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax292.plot(x, yFe54, marker='o', markersize=15, alpha=1.0, mfc=cmap2(rat5), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax292.plot(x, yAu197, marker='^', markersize=15, alpha=1.0, mfc=cmap2(rat6), mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.hlines(y=[vCa48, vFe54, vAu197], xmin=0.0, xmax=0.28, colors=['black', 'black', 'black', 'black'], linestyles=['--', '--', '--', '--'], linewidth=[5,5,5,5])

plt.text(0.7, 1.09, r'$\dfrac{^{48} \mathrm{Ca} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color=cmap2(rat4))
plt.text(0.7, 1.38, r'$\dfrac{^{54} \mathrm{Fe} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color=cmap2(rat5))
plt.text(0.7, 1.67, r'$\dfrac{^{197} \mathrm{Au} }{ ^{40} \mathrm{Ca} }$', fontsize=25, color=cmap2(rat6))

plt.plot(x, yCa48, linestyle='none')
plt.fill_between(x, yCa48-eCa48, yCa48+eCa48, alpha=0.5, edgecolor=cmap2(rat4), facecolor=cmap2(rat4))
plt.plot(x, yFe54, linestyle='none')
plt.fill_between(x, yFe54-eFe54, yFe54+eFe54, alpha=0.5, edgecolor=cmap2(rat5), facecolor=cmap2(rat5))
plt.plot(x, yAu197, linestyle='none')
plt.fill_between(x, yAu197-eAu197, yAu197+eAu197, alpha=0.5, edgecolor=cmap2(rat6), facecolor=cmap2(rat6))

ax292.set_title('', fontsize=21)
ax292.set_xlabel(r'Missing Momentum [GeV/c]', fontsize=25, labelpad=-2)
ax292.set_ylabel(r'SRC Per Proton Ratio to $^{40}$Ca', fontsize=25)

ax292.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.5, 0.6, 0.7])
ax292.set_xlim(0.35, 0.78)
ax292.set_ylim(0.0, 2.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_pmiss_ratio_ca40.png')
plt.close()




















fig, ax292 = plt.subplots(figsize=(10,7))

ax292.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197], 
[0.764, 0.861, 0.996, 1.000, 1.168, 1.330, 1.387, 1.574],
[0.021, 0.023, 0.026, 0.000, 0.037, 0.050, 0.047, 0.094], 
marker='s', markersize=15, alpha=1.0, mfc='black', ecolor='black', 
mec='None', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None',
label='Duer et al', zorder=4
)

plt.xscale('log')

ax292.set_title('', fontsize=25)
ax292.set_xlabel(r'A', fontsize=25)
ax292.set_ylabel('Per Proton Double Ratio to C', fontsize=25)

ax292.tick_params(axis='both', which='major', labelsize=15)
#plt.xticks([0.4, 0.5, 0.6, 0.7])
ax292.set_xlim(6.5, 250)
ax292.set_ylim(0.5, 1.75)
#plt.xticks(fontsize = 22)
plt.yticks([0.5, 0.75, 1.0, 1.25, 1.5, 1.75])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('PP_Double_Ratio_nt.png')
plt.close()











fig, ax292 = plt.subplots(figsize=(10,7))
ax292.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197], 
[1.048, 1.070, 1.005, 1.000, 0.793, 0.764, 0.841, 0.802],
[0.047, 0.031, 0.029, 0.000, 0.036, 0.035, 0.038, 0.046], 
marker='s', markersize=15, alpha=1.0, mfc='black', ecolor='black', 
mec='None', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None',
label='Duer et al', zorder=4
)

plt.xscale('log')

ax292.set_title('', fontsize=25)
ax292.set_xlabel(r'A', fontsize=25)
ax292.set_ylabel('MF Per Proton Ratio to C', fontsize=25)

ax292.tick_params(axis='both', which='major', labelsize=15)
#plt.xticks([0.4, 0.5, 0.6, 0.7])
ax292.set_xlim(6.5, 250)
ax292.set_ylim(0.5, 1.25)
#plt.xticks(fontsize = 22)
plt.yticks([0.5, 0.75, 1.0, 1.25])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('PP_MF_Single_Ratio_nt.png')
plt.close()



















fig, ax292 = plt.subplots(figsize=(10,7))
ax292.errorbar(
[9, 10, 11, 12, 40, 48, 54, 197], 
[0.80, 0.92, 1.00, 1.00, 0.93, 1.02, 1.17, 1.26],
[0.04, 0.04, 0.04, 0.00, 0.05, 0.06, 0.06, 0.10], 
marker='s', markersize=15, alpha=1.0, mfc='black', ecolor='black', 
mec='None', elinewidth=1.5, capsize=4, markeredgewidth=1.2, linestyle='None',
label='Duer et al', zorder=4
)

plt.xscale('log')

ax292.set_title('', fontsize=25)
ax292.set_xlabel(r'A', fontsize=25)
ax292.set_ylabel('SRC Per Proton Ratio to C', fontsize=25)

ax292.tick_params(axis='both', which='major', labelsize=15)
#plt.xticks([0.4, 0.5, 0.6, 0.7])
ax292.set_xlim(6.5, 250)
ax292.set_ylim(0.5, 1.5)
#plt.xticks(fontsize = 22)
plt.yticks([0.5, 0.75, 1.0, 1.25, 1.5])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('PP_SRC_Single_Ratio_nt.png')
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(9.0, 10.0)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8, 10.0])
ax101.set_ylim(0.0, 8.0)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0])
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
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='green', facecolor='green')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='blue', facecolor='blue')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Final e$^{-}$ Momentum [GeV/c]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(9.0, 10.0)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8, 10.0])
ax101.set_ylim(0.0, 7.0)
plt.yticks([0.0, 2.0, 4.0, 6.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_Ef_Cut2.png')
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(9.0, 10.0)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8, 10.0])
ax201.set_ylim(0.0, 6.5)
plt.yticks([0.0, 2.0, 4.0, 6.0])
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(1.0, 1.5)
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
ax101.set_ylim(0.0, 8.0)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0])
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
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='green', facecolor='green')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='blue', facecolor='blue')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Final Proton Momentum [GeV/c]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(1.0, 1.5)
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
ax101.set_ylim(0.0, 7.0)
plt.yticks([0.0, 2.0, 4.0, 6.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_Pf_Cut2.png')
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(1.0, 1.5)
plt.xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
ax201.set_ylim(0.0, 5.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(7.0, 10.0)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
ax101.set_ylim(0.0, 4.0)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0])

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
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='green', facecolor='green')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='blue', facecolor='blue')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Electron Scattering Angle [deg]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(7.0, 10.0)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
ax101.set_ylim(0.0, 3.0)
plt.yticks([0.0, 1.0, 2.0, 3.0])

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_the_Cut2.png')
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(7.0, 10.0)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])
ax201.set_ylim(0.0, 4.0)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_the_Cut.png')
plt.close()




















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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

#sm = mpl.cm.ScalarMappable(cmap="YlOrBr", norm=mpl.colors.Normalize(30, 60))
#cbar = plt.colorbar(sm, ax=ax101, ticks=[], label="A")  
#cbar.set_label('A', size=22)

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax101.set_xlim(1.5, 3.0)
ax101.set_ylim(0.0, 11.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax201.set_xlim(1.5, 3.0)
ax201.set_ylim(0.0, 8.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LSRC_Q2_CutVar.png')
else:
    fig.savefig('LSRC_Q2_Cut.png')
plt.close()



















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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

#sm = mpl.cm.ScalarMappable(cmap="YlOrBr", norm=mpl.colors.Normalize(30, 60))
#cbar = plt.colorbar(sm, ax=ax101, ticks=[], label="A")  
#cbar.set_label('A', size=22)

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax101.set_xlim(1.5, 3.0)
ax101.set_ylim(0.0, 11.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.5, 2.0, 2.5, 3.0])
ax201.set_xlim(1.5, 3.0)
ax201.set_ylim(0.0, 8.0)
plt.xticks(fontsize = 22)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0])
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

if cutvar:
    fig.savefig('LSRC_Q2_CutVar.png')
else:
    fig.savefig('LSRC_Q2_Cut.png')
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([62.0, 64.0, 66.0, 68.0, 70.0])
ax101.set_xlim(62.0, 70.0)
plt.yticks([0.0, 2.0, 4.0, 6.0, 8.0])
ax101.set_ylim(0.0, 8.0)
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
ax101.errorbar(x, y2, marker='o', markersize=15, alpha=1.0, mfc='green', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)
ax101.errorbar(x, y3, marker='^', markersize=15, alpha=1.0, mfc='blue', mec='None', linestyle='None', label=r'$^{54}$Fe', zorder=3)

plt.plot(x, y1, linestyle='none')
plt.fill_between(x, y1-e1, y1+e1, alpha=0.5, edgecolor='red', facecolor='red')
plt.plot(x, y2, linestyle='none')
plt.fill_between(x, y2-e2, y2+e2, alpha=0.5, edgecolor='green', facecolor='green')
plt.plot(x, y3, linestyle='none')
plt.fill_between(x, y3-e3, y3+e3, alpha=0.5, edgecolor='blue', facecolor='blue')

ax101.set_title('', fontsize=21)
ax101.set_xlabel(r'Proton Scattering Angle [deg]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([62.0, 64.0, 66.0, 68.0, 70.0])
ax101.set_xlim(62.0, 70.0)
plt.yticks([0.0, 2.0, 4.0, 6.0])
ax101.set_ylim(0.0, 7.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_thp_Cut2.png')
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([62.0, 64.0, 66.0, 68.0, 70.0])
ax201.set_xlim(62.0, 70.0)
plt.yticks([0.0, 2.0, 4.0, 6.0])
ax201.set_ylim(0.0, 6.0)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_thp_Cut.png')
plt.close()









































ifname= 'SRCC12nu.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197nu.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40nu.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48nu.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54nu.csv'
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
ax101.set_xlabel(r'Energy Transfer ($\omega$) [GeV/c]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.6, 0.8, 1.0, 1.2])
ax101.set_xlim(0.4, 1.2)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
ax101.set_ylim(0.0, 6.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_nu_Cut.png')
plt.close()







ifname= 'SRCBe9nu.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10nu.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11nu.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12nu.csv'
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
ax201.set_xlabel(r'Energy Transfer ($\omega$) [GeV/c]', fontsize=25)
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.4, 0.6, 0.8, 1.0, 1.2])
ax201.set_xlim(0.4, 1.2)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0])
ax201.set_ylim(0.0, 4.5)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_nu_Cut.png')
plt.close()


































ifname= 'SRCC12q.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'SRCAu197q.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'SRCCa40q.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCCa48q.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCFe54q.csv'
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
ax101.set_xlabel(r'3 Momentum Transfer ($\vec{q}$) [GeV/C]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.2, 1.4, 1.6, 1.8, 2.0, 2.2])
ax101.set_xlim(1.2, 2.2)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
ax101.set_ylim(0.0, 5.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('SRC_q_Cut.png')
plt.close()







ifname= 'SRCBe9q.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'SRCB10q.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'SRCB11q.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'SRCC12q.csv'
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
ax201.set_xlabel(r'3 Momentum Transfer ($\vec{q}$) [GeV/C]', fontsize=25)
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.2, 1.4, 1.6, 1.8, 2.0, 2.2])
ax201.set_xlim(1.2, 2.2)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0])
ax201.set_ylim(0.0, 4.5)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LSRC_q_Cut.png')
plt.close()























































ifname= 'MFC12Ef.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197Ef.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40Ef.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48Ef.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54Ef.csv'
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(9.0, 9.8)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8])
ax101.set_ylim(0.0, 1000.0)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800, 1000.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_Ef_Cut.png')
plt.close()







ifname= 'MFBe9Ef.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10Ef.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11Ef.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12Ef.csv'
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(9.0, 9.8)
plt.xticks([9.0, 9.2, 9.4, 9.6, 9.8])
ax201.set_ylim(0.0, 1000.0)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_Ef_Cut.png')
plt.close()





















ifname= 'MFC12Pf.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197Pf.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40Pf.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48Pf.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54Ef.csv'
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
ax201.set_xlabel(r'Final Proton Momentum [GeV/c]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(1.4, 2.2)
plt.xticks([1.4, 1.6, 1.8, 2.0, 2.2])
ax101.set_ylim(0.0, 350.0)
plt.yticks([0.0, 100.0, 200.0, 300.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_Pf_Cut.png')
plt.close()







ifname= 'MFBe9Pf.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10Pf.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11Pf.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12Pf.csv'
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(1.4, 2.2)
plt.xticks([1.4, 1.6, 1.8, 2.0, 2.2])
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0])
ax201.set_ylim(0.0, 900.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_Pf_Cut.png')
plt.close()









"""

ifname= 'MFC12Pf.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197Pf.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40Pf.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48Pf.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54Pf.csv'
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(1.4, 2.2)
plt.xticks([1.4, 1.6, 1.8, 2.0, 2.2])
ax101.set_ylim(0.0, 350.0)
plt.yticks([0.0, 100.0, 200.0, 300.0])
#ax101.set_ylim(0.0, 5.0)
#plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_Pf_Cut.png')
plt.close()








ifname= 'MFBe9Pf.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10Pf.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11Pf.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12Pf.csv'
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(1.4, 2.2)
plt.xticks([1.4, 1.6, 1.8, 2.0, 2.2])
ax201.set_ylim(0.0, 900.0)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0])
ax201.set_ylim(0.0, 900.0)
#plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0])
#ax201.set_ylim(0.0, 5.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_Pf_Cut.png')
plt.close()

"""




















ifname= 'MFC12the.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197the.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40the.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48the.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54the.csv'
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
ax101.set_xlim(7.0, 10.5)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5])
ax101.set_ylim(0.0, 500.0)
plt.yticks([0.0, 100.0, 200.0, 300.0, 400.0, 500.0])

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_the_Cut.png')
plt.close()








ifname= 'MFBe9the.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10the.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11the.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12the.csv'
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
ax201.set_xlim(7.0, 10.5)
plt.xticks([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5])
ax201.set_ylim(0.0, 600.0)
plt.yticks([0.0, 200.0, 400.0, 600.0])
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_the_Cut.png')
plt.close()


































ifname= 'MFC12thp.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197thp.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40thp.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48thp.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54thp.csv'
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
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([44.0, 46.0, 48.0, 50.0, 52.0])
ax101.set_xlim(44.0, 52.0)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0])
ax101.set_ylim(0.0, 900.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_thp_Cut.png')
plt.close()







ifname= 'MFBe9thp.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10thp.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11thp.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12thp.csv'
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
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([44.0, 46.0, 48.0, 50.0, 52.0])
ax201.set_xlim(44.0, 52.0)
plt.yticks([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])
ax201.set_ylim(0.0, 1000.0)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_thp_Cut.png')
plt.close()









































ifname= 'MFC12nu.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197nu.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40nu.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48nu.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54nu.csv'
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
ax101.set_xlabel(r'Energy Transfer ($\omega$) [GeV/c]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.8, 1.0, 1.2, 1.4, 1.6])
ax101.set_xlim(0.8, 1.6)
plt.yticks([0.0, 200.0, 400.0, 600.0])
ax101.set_ylim(0.0, 600.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_nu_Cut.png')
plt.close()







ifname= 'MFBe9nu.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10nu.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11nu.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12nu.csv'
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
ax201.set_xlabel(r'Energy Transfer ($\omega$) [GeV/c]', fontsize=25)
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([0.8, 1.0, 1.2, 1.4, 1.6])
ax201.set_xlim(0.8, 1.6)
plt.yticks([0.0, 200.0, 400.0, 600.0])
ax201.set_ylim(0.0, 700.0)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_nu_Cut.png')
plt.close()


































ifname= 'MFC12q.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y0 = np.array(df['zcont'])
e0 = np.array(df['zcont_err'])

ifname= 'MFAu197q.csv'
df = pd.read_csv(ifname, comment='#')
y4 = np.array(df['zcont'])
e4 = np.array(df['zcont_err'])

ifname= 'MFCa40q.csv'
df = pd.read_csv(ifname, comment='#')
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFCa48q.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFFe54q.csv'
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
ax101.set_xlabel(r'3 Momentum Transfer ($\vec{q}$) [GeV/C]', fontsize=25)
ax101.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax101.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4])
ax101.set_xlim(1.4, 2.4)
plt.yticks([0.0, 200.0, 400.0, 600.0])
ax101.set_ylim(0.0, 700.0)
plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('MF_q_Cut.png')
plt.close()







ifname= 'MFBe9q.csv'
df = pd.read_csv(ifname, comment='#')
x = np.array(df['x0'])
y1 = np.array(df['zcont'])
e1 = np.array(df['zcont_err'])

ifname= 'MFB10q.csv'
df = pd.read_csv(ifname, comment='#')
y2 = np.array(df['zcont'])
e2 = np.array(df['zcont_err'])

ifname= 'MFB11q.csv'
df = pd.read_csv(ifname, comment='#')
y3 = np.array(df['zcont'])
e3 = np.array(df['zcont_err'])

ifname= 'MFC12q.csv'
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
ax201.set_xlabel(r'3 Momentum Transfer ($\vec{q}$) [GeV/C]', fontsize=25)
ax201.set_ylabel('Per Proton Normalized Yield', fontsize=25)#, weight='bold')

ax201.tick_params(axis='both', which='major', labelsize=15)
plt.xticks([1.4, 1.6, 1.8, 2.0, 2.2, 2.4])
ax201.set_xlim(1.4, 2.4)
plt.yticks([0.0, 200.0, 400.0, 600.0])
ax201.set_ylim(0.0, 700.0)

plt.xticks(fontsize = 22)
plt.yticks(fontsize = 22)
fig.set_tight_layout(True)

fig.savefig('LMF_q_Cut.png')
plt.close()