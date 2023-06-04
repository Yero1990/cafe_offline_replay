import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy


ifname= 'double_ratios_pass2.csv'

# read input file
df = pd.read_csv(ifname, comment='#')

targ = df['target']
double_ratio  = np.array(df['doubleR'])
double_ratio_err  = np.array(df['doubleR_err'])
A = df['A'] 
NoZ = df['NoZ'] 
NmZoA = df['NmZoA'] 

fig1 = plt.figure()
plt.errorbar(A, double_ratio, double_ratio_err, marker='o', markersize=10, mfc='k', ecolor='k', mec='k', linestyle='None')
plt.xscale('log')
plt.xlabel('A')
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)')

fig2 = plt.figure()
plt.errorbar(NoZ, double_ratio, double_ratio_err, marker='o', markersize=10, mfc='k', ecolor='k', mec='k', linestyle='None')
plt.xlabel('N/Z')
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)')

fig3 = plt.figure()
plt.errorbar(NmZoA, double_ratio, double_ratio_err, marker='o', markersize=10, mfc='k', ecolor='k', mec='k', linestyle='None')
plt.xlabel('(N-Z)/A')
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)')

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
