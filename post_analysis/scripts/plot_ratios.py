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

# A_SRC / A_MF
singleR_per_proton      = np.array(df['singleR_per_proton'])
singleR_per_proton_err  = np.array(df['singleR_per_proton_err'])

# A_MF / C12_MF
singleR_A_c12_mf        = np.array(df['singleR_A_c12_mf'])
singleR_A_c12_mf_err    = np.array(df['singleR_A_c12_mf_err'])

# A_SRC / C12_SRC
singleR_A_c12_src        = np.array(df['singleR_A_c12_src'])
singleR_A_c12_src_err    = np.array(df['singleR_A_c12_src_err'])

doubleR     = np.array(df['doubleR'])
doubleR_err  = np.array(df['doubleR_err'])


A = df['A'] 
NoZ = df['NoZ'] 
NmZoA = df['NmZoA'] 


compare_flag = False  # compare to previous pass ?
compare_simc_flag = True  #compare to simc ?

if(compare_simc_flag):

    # read SIMC single_ratio (C12, Fe56, Au197)
    df_simc = pd.read_csv('cafe_simc_summary.csv', comment='#')
    A_simc = df_simc['A']
    singleR_per_proton_simc      = np.array(df_simc['singleR_A_to_C12_MF'])
    singleR_per_proton_err_simc  = np.array(df_simc['singleR_A_to_C12_MF_err'])




if(compare_flag):
    fname= 'cafe_ratios_pass2.csv' 

    # read input file
    df = pd.read_csv(fname, comment='#')
    
    targ_2 = df['target']
    singleR_per_proton_2      = np.array(df['singleR_per_proton'])
    singleR_per_proton_err_2  = np.array(df['singleR_per_proton_err'])

    singleR_A_c12_mf_2        = np.array(df['singleR_A_c12_mf'])
    singleR_A_c12_mf_err_2    = np.array(df['singleR_A_c12_mf_err'])
    
    singleR_A_c12_src_2        = np.array(df['singleR_A_c12_src'])
    singleR_A_c12_src_err_2    = np.array(df['singleR_A_c12_src_err'])

    doubleR_2      = np.array(df['doubleR'])
    doubleR_err_2  = np.array(df['doubleR_err'])
    A_2 = df['A'] 
    NoZ_2 = df['NoZ'] 
    NmZoA_2 = df['NmZoA'] 


#--------------------------
# PLOT Double Ratio vs. A
#--------------------------

fig1= plt.figure()
plt.errorbar(A, doubleR, doubleR_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))
if(compare_flag ):
    plt.errorbar(A_2, doubleR_2, doubleR_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')

plt.title('CaFe Double Ratio vs. A', fontsize=18)
plt.xlabel('A', fontsize=16)
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)', fontsize=16)
# add target names to plot
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',A[i], doubleR[i] )

    # standard (x,y) coordinates
    x = A[i] + 2
    y = doubleR[i]

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
        y = doubleR[i] + 0.05
    elif tgt=="Be9" or tgt=="B10" or tgt=="C12":
        x = A[i] + 1
    else:
        x = A[i] + 2

    plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_doubleR_vs_A.pdf')




#--------------------------
# PLOT Single Ratio vs. A
#--------------------------

#--------------------------
# Single Ratio A_SRC / A_MF 
#--------------------------
fig1_s= plt.figure()
plt.errorbar(A, singleR_per_proton, singleR_per_proton_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))
if(compare_flag ):
    plt.errorbar(A_2, singleR_per_proton_2, singleR_per_proton_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')

plt.title('CaFe Single Ratio (per proton) vs. A', fontsize=18)
plt.xlabel('A', fontsize=16)
plt.ylabel('A (SRC/MF)', fontsize=16)
# add target names to plot
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',A[i], singleR_per_proton[i] )

    # standard (x,y) coordinates
    x = A[i] + 2
    y = singleR_per_proton[i]

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
        y = singleR_per_proton[i] +0.0005
    elif tgt=="Be9" or tgt=="B10" or tgt=="C12":
        x = A[i] + 1
    else:
        x = A[i] + 2

    plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_singleR_vs_A.pdf')

# --------------------------------------------

#-----------------------------
# Single Ratio A_MF / C12_MF 
#-----------------------------
fig1_a2c_mf= plt.figure()
plt.errorbar(A, singleR_A_c12_mf, singleR_A_c12_mf_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))
if(compare_flag ):
    plt.errorbar(A_2, singleR_A_c12_mf_2, singleR_A_c12_mf_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
if(compare_simc_flag ):
    plt.errorbar(A_simc, singleR_per_proton_simc, singleR_per_proton_err_simc, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='SIMC')
        
plt.xscale('log')
plt.title('CaFe Single Ratio (per proton) vs. A', fontsize=18)
plt.xlabel('A', fontsize=16)
plt.ylabel('A_MF / C12_MF', fontsize=16)


# add target names to plot
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',A[i], singleR_A_c12_mf[i] )

    # standard (x,y) coordinates
    x = A[i] + 2
    y = singleR_A_c12_mf[i] 

    
   
    if tgt=="Au197":
        x = A[i] - 60
    elif tgt=="Be9" or tgt=="B10" or tgt=="B11" :
        x = A[i] + 0.8
    elif tgt=="C12":
        x = A[i] + 0.8
        y = y - 0.02
    elif tgt=="Ca40" or tgt=="Ca48" or tgt=="Fe54" :
        x = A[i] + 1.8
        
    plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_MFsingleR_vs_A.pdf')
# --------------------------------------------


#-----------------------------
# Single Ratio A_SRC / C12_SRC 
#-----------------------------
fig1_a2c= plt.figure()
plt.errorbar(A, singleR_A_c12_src, singleR_A_c12_src_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))
if(compare_flag ):
    plt.errorbar(A_2, singleR_A_c12_src_2, singleR_A_c12_src_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')
    
plt.xscale('log')
plt.title('CaFe Single Ratio (per proton) vs. A', fontsize=18)
plt.xlabel('A', fontsize=16)
plt.ylabel('A_SRC / C12_SRC ', fontsize=16)


# add target names to plot
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',A[i], singleR_A_c12_src[i] )

    # standard (x,y) coordinates
    x = A[i] + 2
    y = singleR_A_c12_src[i]

    if tgt=="Fe54":
        x = A[i] + 6
    elif tgt=="Ca48":
        x = A[i] + 5
    elif tgt=="Ca40":
        x = A[i] + 4
    elif tgt=="Au197":
        x = A[i] - 60
    elif tgt=="B11":
        x = A[i] - 3
        y = singleR_A_c12_src[i] +0.01
    elif tgt=="Be9" or tgt=="B10" or tgt=="C12":
        x = A[i] + 2
    else:
        x = A[i] + 2

    plt.text(x, y, tgt)
    
plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_SRCsingleR_vs_A.pdf')
# --------------------------------------------


#--------------------------
# PLOT Single Ratio vs. N/Z
#--------------------------

# A_SRC / A_MF
fig2_s = plt.figure()
plt.errorbar(NoZ, singleR_per_proton, singleR_per_proton_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NoZ_2, singleR_per_proton_2, singleR_per_proton_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.title('CaFe Single Ratio vs. N/Z', fontsize=18)
plt.xlabel('N/Z', fontsize=16)
plt.ylabel('A (SRC/MF) ', fontsize=16)
for i, tgt in enumerate(targ):
    
    x = NoZ[i] + 0.01
    y = singleR_per_proton[i]

    if tgt=="Au197":
        x = NoZ[i] - 0.07
 
    plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_singleR_vs_NoZ.pdf')

# A_SRC / C12_SRC
fig2_s_src = plt.figure()
plt.errorbar(NoZ, singleR_A_c12_src, singleR_A_c12_src_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NoZ_2, singleR_A_c12_src_2, singleR_A_c12_src_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.title('CaFe Single Ratio vs. N/Z', fontsize=18)
plt.xlabel('N/Z', fontsize=16)
plt.ylabel('A_SRC / C12_SRC', fontsize=16)
for i, tgt in enumerate(targ):
    
    x = NoZ[i] + 0.01
    y = singleR_A_c12_src[i]
 
    plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_SRCsingleR_vs_NoZ.pdf')





#--------------------------
# PLOT Double Ratio vs. N/Z
#--------------------------
fig2 = plt.figure()
plt.errorbar(NoZ, doubleR, doubleR_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NoZ_2, doubleR_2, doubleR_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.title('CaFe Double Ratio vs. N/Z', fontsize=18)
plt.xlabel('N/Z', fontsize=16)
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)', fontsize=16)
for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',NoZ[i], doubleR[i] )
    
    x = NoZ[i] + 0.01
    y = doubleR[i]

    if tgt=="Au197":
        x = NoZ[i] - 0.07
 
    plt.text(x, y, tgt)

plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_doubleR_vs_NoZ.pdf')










#-------------------------------
# PLOT Double Ratio vs. (N-Z)/A
#-------------------------------
fig3 = plt.figure()
plt.errorbar(NmZoA, doubleR, doubleR_err, marker='o', markersize=7, mfc='k', ecolor='k', mec='k', linestyle='None', label='%s'%(npass))

if(compare_flag ):
    plt.errorbar(NmZoA_2, doubleR_2, doubleR_err_2, marker='o', markersize=7, mfc='r', ecolor='r', mec='r', linestyle='None', label='pass2')

plt.title('CaFe Double Ratio vs. (N-Z)/A', fontsize=18)
plt.xlabel('(N-Z)/A', fontsize=16)
plt.ylabel('A (SRC/MF) / C12 (SRC/MF)', fontsize=16)

for i, tgt in enumerate(targ):
    print('i, tgt -> ',i, tgt)
    print('(x,y) ->  ',NmZoA[i], doubleR[i] )

    x = NmZoA[i] + 0.005
    y = doubleR[i]

    if tgt=="Au197":
        x = NmZoA[i] - 0.03
 
    plt.text(x, y, tgt)


plt.legend(frameon=False, fontsize=16)
plt.savefig('cafe_doubleR_vs_NmZoA.pdf')





plt.show()









