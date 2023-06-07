import numpy as np
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


ifname= 'cafe_ratios_%s.csv' % (npass)

# read input file
df = pd.read_csv(ifname, comment='#')

targ = df['target']

doubleR_per_nucleon  = np.array(df['doubleR_per_nucleon'])
doubleR_per_nucleon_err  = np.array(df['doubleR_per_nucleon_err'])


A = df['A'] 
NoZ = df['NoZ'] 
NmZoA = df['NmZoA'] 

compare_flag = True

if(compare_flag):
    fname= 'cafe_ratios_pass2.csv' 

    # read input file
    df = pd.read_csv(fname, comment='#')
    
    targ_2 = df['target']
    doubleR_per_nucleon_2  = np.array(df['doubleR_per_nucleon'])
    doubleR_per_nucleon_err_2  = np.array(df['doubleR_per_nucleon_err'])
    A_2 = df['A'] 
    NoZ_2 = df['NoZ'] 
    NmZoA_2 = df['NmZoA'] 


#--------------------------
# PLOT Double Ratio vs. A
#--------------------------

fig1= plt.figure()
plt.errorbar(A, doubleR_per_nucleon, doubleR_per_nucleon_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))
if(compare_flag ):
    plt.errorbar(A_2, doubleR_per_nucleon_2, doubleR_per_nucleon_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')
plt.xlabel('A')
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)')
# add target names to plot
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',A[i], doubleR_per_nucleon[i] )
    
    
    if(A[i]>15):
        x = A[i] + 5
        y = doubleR_per_nucleon[i]
    else:
        x = A[i] + 2
        y = doubleR_per_nucleon[i]

    if tgt=="Au197":
        x = A[i] - 70
    if tgt=="B11":
        x = A[i] - 2
        y = doubleR_per_nucleon[i] + 0.05
    plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)


#--------------------------
# PLOT Double Ratio vs. N/Z
#--------------------------
fig2 = plt.figure()
plt.errorbar(NoZ, doubleR_per_nucleon, doubleR_per_nucleon_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NoZ_2, doubleR_per_nucleon_2, doubleR_per_nucleon_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.xlabel('N/Z')
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)')
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',NoZ[i], doubleR_per_nucleon[i] )
    
    x = NoZ[i] + 0.01
    y = doubleR_per_nucleon[i]

    if tgt=="Au197":
        x = NoZ[i] - 0.07
 
    plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)





fig3 = plt.figure()
plt.errorbar(NmZoA, doubleR_per_nucleon, doubleR_per_nucleon_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NmZoA_2, doubleR_per_nucleon_2, doubleR_per_nucleon_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.xlabel('(N-Z)/A')
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)')

for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',NmZoA[i], doubleR_per_nucleon[i] )

    x = NmZoA[i] + 0.005
    y = doubleR_per_nucleon[i]

    if tgt=="Au197":
        x = NmZoA[i] - 0.03
 
    plt.text(x, y, tgt)


plt.legend(frameon=False, fontsize=16)





plt.show()













'''
#print('double_ratio = ', double_ratio)
for i in range(len(src_yield_norm_arr)):

       
        #print('N[i] = ', N[i])
        if (targ[i]=='LD2'):
            continue
        if (targ[i]=='Ca48' or targ[i]=='B10' or targ[i]=='B11'):
            plt.errorbar(N[i]/Z[i], double_ratio_val[i], double_ratio_err[i], marker='o', markersize=10, mfc='gray', ecolor='gray', mec='gray', linestyle='None', label=targ[i])
        else:
            plt.errorbar(N[i]/Z[i], double_ratio_val[i], double_ratio_err[i], marker='o', markersize=10, mfc=tcolor[i], ecolor=tcolor[i], mec='k', linestyle='None', label=targ[i])
            #plt.errorbar((N[i]-Z[i])/Z[i], double_ratio_val[i], double_ratio_err[i], marker='o', markersize=10, mfc=tcolor[i], ecolor=tcolor[i], mec='k', linestyle='None', label=targ[i])
        
        
    plt.ylabel('SRC High Momentum Fraction', fontsize=18)
    plt.xlabel('N/Z', fontsize=18)
    #plt.xlabel('(N-Z)/Z', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.legend()
    plt.show()
'''
