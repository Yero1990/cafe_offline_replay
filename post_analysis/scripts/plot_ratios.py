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

# ----------------
# read input file
# ----------------
df = pd.read_csv(ifname, comment='#')

targ = df['target']

singleR_per_nucleon      = np.array(df['singleR_per_nucleon'])
singleR_per_nucleon_err  = np.array(df['singleR_per_nucleon_err'])

singleR_A_c12_src        = np.array(df['singleR_A_c12_src'])
singleR_A_c12_src_err    = np.array(df['singleR_A_c12_src_err'])

doubleR_per_nucleon      = np.array(df['doubleR_per_nucleon'])
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
    singleR_per_nucleon_2      = np.array(df['singleR_per_nucleon'])
    singleR_per_nucleon_err_2  = np.array(df['singleR_per_nucleon_err'])

    singleR_A_c12_src_2        = np.array(df['singleR_A_c12_src'])
    singleR_A_c12_src_err_2    = np.array(df['singleR_A_c12_src_err'])

    doubleR_per_nucleon_2      = np.array(df['doubleR_per_nucleon'])
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

plt.title('CaFe Double Ratio vs. A', fontsize=18)
plt.xlabel('A', fontsize=16)
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)', fontsize=16)
# add target names to plot
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',A[i], doubleR_per_nucleon[i] )

    # standard (x,y) coordinates
    x = A[i] + 2
    y = doubleR_per_nucleon[i]

    if tgt=="Fe54":
        x = A[i] - 15
    elif tgt=="Ca48":
        x = A[i] - 13
    elif tgt=="Ca40":
        x = A[i] - 12
    elif tgt=="Au197":
        x = A[i] - 70
    elif tgt=="B11":
        x = A[i] - 2
        y = doubleR_per_nucleon[i] + 0.05
    elif tgt=="Be9" or tgt=="B10" or tgt=="C12":
        x = A[i] + 1
    else:
        x = A[i] + 2

    plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)




#--------------------------
# PLOT Single Ratio vs. A
#--------------------------

fig1_s= plt.figure()
plt.errorbar(A, singleR_per_nucleon, singleR_per_nucleon_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))
if(compare_flag ):
    plt.errorbar(A_2, singleR_per_nucleon_2, singleR_per_nucleon_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')

plt.title('CaFe Single Ratio vs. A', fontsize=18)
plt.xlabel('A', fontsize=16)
plt.ylabel('A (SRC/MF)', fontsize=16)
# add target names to plot
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',A[i], singleR_per_nucleon[i] )

    # standard (x,y) coordinates
    x = A[i] + 2
    y = singleR_per_nucleon[i]

    if tgt=="Fe54":
        x = A[i] + 15
    elif tgt=="Ca48":
        x = A[i] + 13
    elif tgt=="Ca40":
        x = A[i] + 12
    elif tgt=="Au197":
        x = A[i] - 70
    elif tgt=="B11":
        x = A[i] 
        y = singleR_per_nucleon[i] +0.0005
    elif tgt=="Be9" or tgt=="B10" or tgt=="C12":
        x = A[i] + 1
    else:
        x = A[i] + 2

    plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)

fig1_a2c= plt.figure()
plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))
if(compare_flag ):
    plt.errorbar(A_2, singleR_A_c12_src_2, singleR_A_c12_src_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')
plt.title('CaFe Single Ratio vs. A', fontsize=18)
plt.xlabel('A', fontsize=16)
plt.ylabel('A_SRC / C12_SRC', fontsize=16)

plt.legend(frameon=False, fontsize=16)




#--------------------------
# PLOT Double Ratio vs. N/Z
#--------------------------
fig2 = plt.figure()
plt.errorbar(NoZ, doubleR_per_nucleon, doubleR_per_nucleon_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NoZ_2, doubleR_per_nucleon_2, doubleR_per_nucleon_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.title('CaFe Double Ratio vs. N/Z', fontsize=18)
plt.xlabel('N/Z', fontsize=16)
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)', fontsize=16)
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',NoZ[i], doubleR_per_nucleon[i] )
    
    x = NoZ[i] + 0.01
    y = doubleR_per_nucleon[i]

    if tgt=="Au197":
        x = NoZ[i] - 0.07
 
    plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)




#--------------------------
# PLOT Single Ratio vs. N/Z
#--------------------------
fig2_s = plt.figure()
plt.errorbar(NoZ, singleR_per_nucleon, singleR_per_nucleon_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NoZ_2, singleR_per_nucleon_2, singleR_per_nucleon_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.title('CaFe Single Ratio vs. N/Z', fontsize=18)
plt.xlabel('N/Z', fontsize=16)
plt.ylabel('A (SRC/MF) ', fontsize=16)
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',NoZ[i], doubleR_per_nucleon[i] )
    
    x = NoZ[i] + 0.01
    y = singleR_per_nucleon[i]

    if tgt=="Au197":
        x = NoZ[i] - 0.07
 
    plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)






fig3 = plt.figure()
plt.errorbar(NmZoA, doubleR_per_nucleon, doubleR_per_nucleon_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NmZoA_2, doubleR_per_nucleon_2, doubleR_per_nucleon_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.title('CaFe Double Ratio vs. (N-Z)/A', fontsize=18)
plt.xlabel('(N-Z)/A', fontsize=16)
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)', fontsize=16)

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









