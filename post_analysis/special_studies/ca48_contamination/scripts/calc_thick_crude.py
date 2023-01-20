'''
C. Yero
October 11, 2022

Brief: this script uses results from a crude calculation using  measured data and 
SIMC quantities from Ca48 MF runs to extract the areal thickness from contamination 
(assuming pure mineral oil) for now, using the alkene components Cn_H_2n+2 (e.g. C15_H32) 
so the areal thickness for each  of these elements is calculated separately by scaling from 
the measured H-contamination. Use a lower, upper bound (C15_H32, C50, H102) 
'''

import sys
import numpy as np

#user input (run number)
run = int(sys.argv[1])

# the N_data counts have been corrected (or scaled up to account for signal lost due to tight cut and
# radiative effects (Y_nora / Y_rad) intergrated over missing energy (these scale factors may have uncertainties
# and will need to be accounted for.

N_data=-1  # data counts (integrated W peak with |Emiss|<20 MeV and Pmiss<30 MeV cuts
N_simc=-1  # charge-normalized SIMC counts (counts / mC)
Q_data=-1  # data charge [mC]

if run==16978:
    N_data = 1156   # ep data elastic counts
    Q_data = 5.485  # data measured charge [mC]
    N_simc = 194.3  # counts/mC
elif run==16979:
    N_data = 12780
    Q_data = 76.223
    N_simc = 154.6
elif run==17093:
    N_data = 1756
    Q_data = 49.725
    N_simc = 32.31
elif run==17094:
    N_data = 377.8
    Q_data = 12.675
    N_simc = 27.18 
elif run==17096:
    N_data = 2914
    Q_data = 98.943
    N_simc = 27.02 
else:
    print('Please enter valid Ca48 MF run')
    sys.exit()

    
e=1.602e-19              # elementary charge C
dsig_dom = 2.409e-32     # ep cross section [cm^2/sr] @ incident e- energy: 10.549 GeV, SHMS: 8.3 deg, 9.438 GeV
omega=0.000378           # SHMS e- solid angle [sr] determined from dOmega_e = dOmega_p * dxptar_e/dxptar_p * dyptar_e/dyptar_p
Hmol = 1.00794           # H molar mass g/mol
Cmol = 12.011            # C molar mass g/mol
Omol = 15.999            # O molar mass g/mol
Na = 6.0221408e23        # Avogadro's number
Ca48_thick = 1.051       # Ca48 target thickness [g/cm2] from D. Meekins

# atoms / cm^2
nx = N_data / ((Q_data*0.001/e)*dsig_dom*omega)   # areal density of target (quantiy we are trying to extract from data)

# H- thick [g/cm^2]
H_thick = nx * Hmol / Na

# H- thick [g/cm^2] corrected (scale previous thickenss by dataYield / simcYield) to match SIMC counts to DATA counts
H_thick_corr =  H_thick *   ( N_data / (Q_data) ) / N_simc

# scale H -thickness to other elements (depends on which chemical composition it is assumed)
# assuming Olive Oil: C88_H164_O10

# C- thick [g/cm^2]
C_thick = H_thick_corr * (Cmol/Hmol) * 88./164.  # do we also scale by number of nucleons ( i.e., 88 atoms * 12 nucleons/atom)

# O- thick [g/cm^2]
O_thick = H_thick_corr * (Omol/Hmol) * 10./164.  # do we also scale by number of nucleons ( i.e., 10 atoms * 16 nucleons/atom)

# relative fraction of H thickness to Ca48 thickness
H_thick_rel = H_thick / Ca48_thick  * 100.

# relative fraction of C-thickness to Ca48 thickness
C_thick_rel = C_thick / Ca48_thick * 100.

# relative fraction of C-thickness to Ca48 thickness
O_thick_rel = O_thick / Ca48_thick * 100.


print('Ca48 Run {}'.format(run))
print('Contamination Chemical Composition:\n(Olive Oil) C88_H164_O10')
print('nx [H-atoms/cm^2]     = {0:.3E}'.format(nx))
print('H_thick [g/cm^2]      = {0:.3E}'.format(H_thick))
print('H_thick_corr [g/cm^2] = {0:.3E}'.format(H_thick_corr))
print('C_thick [g/cm^2]      = {0:.3E}'.format(C_thick))
print('O_thick [g/cm^2]      = {0:.3E}'.format(O_thick))
print('-----------------------------')
print('contamination thickness \n relative to total Ca48 (1.051 g/cm2)')
print('-----------------------------')
print('H_thick_corr/Ca48_thick [%] ={0:.5f}'.format(H_thick_rel))
print('C_thick/Ca48_thick [%]      ={0:.5f}'.format(C_thick_rel))
print('O_thick/Ca48_thick [%]      ={0:.5f}'.format(O_thick_rel))
