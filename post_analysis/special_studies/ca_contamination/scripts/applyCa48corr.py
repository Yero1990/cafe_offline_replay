import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame as pdf

# define exponential function (based on our fit to T2 scalers vs Charge)
def rel_Cntm_src(Q):
    
    # used to parametrize relative drop
    # in overall contamination (units in percent)
    # for now we assume contamination from C_n+H_2n
    A = 3.090*1e-2
    alpha = 6.422*1e-4
    C = 9.728*1e-01

    return (A*np.exp(-alpha*Q) + C)*100 #percent
    

# H,C absolute contamination (used as initial parameter)
# obtained from CaFe Ca48 MF run analysis (in percent)
# assumes pure mineral oil ( mostly consisting of (C_n+H_2n) )
# absolute (H,C)-contam = (H,C)-thickness / Ca48-thickness
def absCntm(run, element=''):
    # in percent
    if(run==16978): 
        if(element=='H'):
            return 0.65 
        elif(element=='C'):
            return 3.8 
    if(run==16979):
        if(element=='H'):
            return 0.51 
        elif(element=='C'):
            return 3.1
    if(run==17093):
        if(element=='H'):
            return 0.11 
        elif(element=='C'):
            return 0.65
    if(run==17094):
        if(element=='H'):
            return 0.091 
        elif(element=='C'):
            return 0.54
    if(run==17096):
        if(element=='H'):
            return 0.09 
        elif(element=='C'):
            return 0.53
    else:
        return np.nan


fname_ca48_MF='../../../summary_files/pass1/cafe_prod_Ca48_MF_report_summary_contaminated.csv'
fname_ca48_SRC='../../../summary_files/pass1/cafe_prod_Ca48_SRC_report_summary.csv'

df_ca48_mf  = pd.read_csv(fname_ca48_MF, comment='#')
df_ca48_src = pd.read_csv(fname_ca48_SRC, comment='#')

run_ca48_mf  = np.array(df_ca48_mf['run'])
run_ca48_src = np.array(df_ca48_src['run'])

Q_ca48_mf  = np.array(df_ca48_mf['charge'])
Q_ca48_src = np.array(df_ca48_src['charge'])


# total ca48 runs in order
run = np.concatenate((run_ca48_mf[:2], run_ca48_src, run_ca48_mf[-3:]))
Qcharge = np.concatenate((Q_ca48_mf[:2], Q_ca48_src, Q_ca48_mf[-3:]))

Qprev = 0  # previous charge (for cumulative sum)

# open file to write factor
myfile = open('ca48_correction.csv', 'w')
myfile.write('# Ca48 Contamination Correction Factors\n'
             '# \n'
             '# Header Definitions: \n'
             '# run      : run number \n'
             '# kin      : kinematics analyzed  \n'
             '# Q        : charge [mC]  \n'
             '# Qsum     : cumulative charge [mC]  \n'
             '# H_absCntm_meas: Hydrogen measured absolute contamination (in percent)\n'
             '# C_absCntm_meas: Carbon measured absolute contamination (in percent)\n'
             '# C_absCntm_calc: Carbon calculated (from fit to scalers)  absolute contamination (in percent)\n'
)

#myfile.write('idx, kin, run, Q, Qsum, H_absCntm_meas, C_absCntm_meas, C_absCntm_calc\n')
myfile.write('run,kin,Q,Qsum,H_absCntm_meas,C_absCntm_meas,C_absCntm_calc\n')

# loop over each Ca48 run and index it
for idx, irun in enumerate(run):

    #get cumulative charge
    Qsum = Qcharge[idx] + Qprev

    # total charge up to run 16979 (basically just add the charge of the first 2 MF Ca48 runs)
    #Qsum_16979 = Qbcm1[0] +  Qbcm1[1]
    Qsum_firstTwo =  Qcharge[0] +  Qcharge[1]
    
    #get absolute contamination from Ca48 MF (for C and H)
    H_absCntm_meas = absCntm(irun, 'H')
    C_absCntm_meas = absCntm(irun, 'C')

    if(irun == 16978 or irun == 16979 or irun == 17093 or irun == 17094 or irun == 17096):
        kin='MF'
    else:
        kin='SRC'
    
    #===== calculating absolute contamination in the SRC kinematics ======
    # ---- original way ----
    # since H-contamination initially is ~0.6 % and last Ca48 MF run ~0.09% (it is negligible in the relative drop)
    #estimate C absolute contamination from Ca48 SRC using parametrized T2scalers / mC expn. fit
    
    #C_absCntm_SRC = absCntm(16979, 'C') - ( rel_Cntm_src(Qsum_firstTwo) - rel_Cntm_src(Qsum) )
    #if(kin[idx].strip()=="MF"):
    #    C_absCntm_SRC = np.nan

    #----- new way  ------
    C_absCntm_calc = (rel_Cntm_src(Qsum) - rel_Cntm_src(1e9))/rel_Cntm_src(1e9) * 100.
    print('rel_Cntm_src(1e9) = ',rel_Cntm_src(1e9))
    
    #print("%i, %s, %i, %.3f, %.3f, %.3f, %.3f, %.3f\n" % (idx, kin[idx], irun, Qbcm1[idx], Qsum, H_absCntm_MF, C_absCntm_MF, C_absCntm_SRC))
    line="%i,%s,%.3f,%.3f,%.3f,%.3f,%.3f\n" % (irun,  kin, Qcharge[idx], Qsum, H_absCntm_meas, C_absCntm_meas, C_absCntm_calc)
    myfile.write(line)

    Qprev = Qsum

