import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy

'''
Author: C. Yero
Date: Jan 02, 2023
Brief: Compilation of scripts for 
handling/combining cafe numerical 
(.csv) summary files
'''

def find_param(param='', fname=''):

    
    with open(fname) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(param) != -1:
                #print('param: ', line)
                param_value = float(line.split(':')[1].strip())
                return param_value


def make_final_summary():

    #brief: function that reads summary files (.csv) to:
    # 1. apply efficiency corrections
    # 2. sum over all efficiency corrected counts
    # 3. normalize  efficiency-corrected counts by charge (mC), transparency and target density (g/cm^2)
    # 4. write combined numerical quantitied to a .csv file
    #  ( for additional subtractions -> boron-carbide subtractions, Ca48 impurity and contamination corrections, and double ratio calculations)
    
    #output file to write summary file
    ofname = 'cafe_final_summary.csv' 
    ofile = open(ofname, 'r+')
    ofile.write('# CaFe Final Summary File (Pass 1) \n')
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# target   : target name analyzed \n'
                '# kin      : kinematics analyzed  \n'
                )
    ofile.write('target,kin,beam_time,avg_current,total_charge,yield_raw,yield_raw_err,yield_corr,yield_corr_err,yield_norm,yield_norm_err,tgt_area_density,T,N,Z,A\n') 

    # target, kin list
    target = ['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54']
    tcolor  = ['c',    'm',   'r',   'g',   'b', 'darkorange', 'violet', 'gold'] #'darkgray']
    kin    = ['MF', 'SRC']
    kmarker =['o',  's']
 
    for idx in np.arange(len(target)):
        
        for jdx in np.arange(len(kin)):
            
            # get generic summary file name
            summary_file_path = 'summary_files/pass1/cafe_prod_%s_%s_report_summary.csv' % (target[idx], kin[jdx])
            
            # read .csv file
            df = pd.read_csv(summary_file_path, comment='#')

            # select specific Ca48 MF runs (ignore contaminated runs)
            #if( target[idx]=='Ca48' and kin[jdx]=='MF' ):
                #df = df[ df['run']==17096 ] # for now, last run is good enough
            #    df = df.tail(3) # select last 3 runs
            #    print(df)
            # read parameters from .csv file
            T                 = find_param('transparency', summary_file_path) # transparency
            tgt_areal_density = find_param('target_areal_density', summary_file_path) # g/cm^2
            N                 = find_param('N:', summary_file_path) # number of neutrons
            Z                 = find_param('Z:', summary_file_path) # number of protons
            A                 = find_param('A:', summary_file_path) # number of nucleons
            
            # read selected data columns (with respective uncertainties)
            run          = df['run']
            charge       = df['charge'] # [mC]
            beam_time    = df['beam_time'] # (beam-on-target time) [sec]
            avg_current  = df['avg_current'] # average beam current [uA]
            real_Yield   = unumpy.uarray(df['real_Yield'],      df['real_Yield_err'])
            hms_trk_eff  = unumpy.uarray(df['hTrkEff'],         df['hTrkEff_err'])
            shms_trk_eff = unumpy.uarray(df['pTrkEff'],         df['pTrkEff_err'])
            total_LT     = unumpy.uarray(df['tLT'],             df['tLT_err_Bi'])
            mult_trk_eff = unumpy.uarray(df['multi_track_eff'], df['multi_track_eff_err'])
            
            T1_scl_rate  = df['T1_scl_rate']  # SHMS  3/4    khZ
            T2_scl_rate  = df['T2_scl_rate']  # SHMS EL-REAL khZ
            T3_scl_rate  = df['T3_scl_rate']  # HMS 3/4      khZ
            T5_scl_rate  = df['T5_scl_rate']  # COIN         khZ

            # calculate variables from data columns
            Qsum = charge.cumsum()  # calculate cumulative charge (may be helpful for plotting)
            T1_scl = T1_scl_rate * 1000 * beam_time   # total scaler counts
            T2_scl = T2_scl_rate * 1000 * beam_time
            T3_scl = T3_scl_rate * 1000 * beam_time
            T5_scl = T5_scl_rate * 1000 * beam_time 
   
            
            # sum over all counts (before applying any corrections, inefficinecy, charge, etc.) 
            real_yield_total = real_Yield.sum()
            
            # apply efficiency corrections to real yield (with uncertasinties included) (run-by-run)
            real_Yield_corr = real_Yield / (hms_trk_eff * shms_trk_eff * total_LT * mult_trk_eff) 

           
             
            # sum over all efficiency-corrected counts 
            real_Yield_corr_total = real_Yield_corr.sum()
            
            # total charge
            total_charge = charge.sum()

            # total beam time
            total_beam_time = beam_time.sum()

            # total average current
            total_avg_current = np.average(avg_current)
            
            # normalize yield by : total charge, transparency and target density
            yield_norm =  real_Yield_corr_total / (total_charge * T * tgt_areal_density)  # counts / (mC * g/cm^2)

                     
            #------------------MAKE RELEVANT PLOTS----------------
            T2_scl_per_Q = T2_scl/ (charge * T * tgt_areal_density)
            real_yield_per_Q  = real_Yield/(charge * tgt_areal_density)
            real_yield_corr_per_Q  = real_Yield_corr/(charge * T * tgt_areal_density)
            
            #print( real_yield_corr_per_Q)
            #print(target[idx],', ', kin[jdx])
            #print(T2_scl_per_Q/T2_scl_per_Q[0])
            
            
            # plot 1: relative scalers vs current (relative to 1st data point) --> woule be best to normalize both SRC, MF to same 1st point, since T2 scalers (shms same location, so should not change)
            #if(kin[jdx]=='MF' and (target[idx]=='Ca48' or target[idx]=='Ca40' or target[idx]=='Fe54')):
            #if(kin[jdx]=='MF' and (target[idx]=='Be9' or target[idx]=='B10' or target[idx]=='B11' or target[idx]=='C12')):
            #if(kin[jdx]=='SRC' and (target[idx]=='Ca48' or target[idx]=='Ca40' or target[idx]=='Fe54')):
            #if(kin[jdx]=='SRC' and (target[idx]=='Be9' or target[idx]=='B10' or target[idx]=='B11' or target[idx]=='C12')):

                # rel yield vs Qsum
                #plt.errorbar(Qsum,  unumpy.nominal_values(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, unumpy.std_devs(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='dashed', label='%s %s'%(target[idx], kin[jdx]))
                #plt.errorbar(Qsum,  avg_current, marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='dashed', label='%s %s'%(target[idx], kin[jdx]))

                # rel yield ( ir T2 scalers) vs. avg current (or rates)
                #plt.errorbar(avg_current,  unumpy.nominal_values(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, unumpy.std_devs(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, marker='o',markersize=8, color=tcolor[idx], mec='k', linestyle='dashed', label='%s %s (yield)'%(target[idx], kin[jdx]))
                #plt.plot(avg_current,  unumpy.nominal_values(T2_scl_per_Q)/T2_scl_per_Q[0], marker='s', color=tcolor[idx], mec='k', linestyle='solid', markersize=8, label='%s %s (T2 e- scalers)'%(target[idx], kin[jdx]))

                #plt.errorbar(T2_scl_rate,  unumpy.nominal_values(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, unumpy.std_devs(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, marker='o',markersize=8, color=tcolor[idx], mec='k', linestyle='dashed', label='%s %s (yield)'%(target[idx], kin[jdx]))
                #plt.plot(T2_scl_rate,  unumpy.nominal_values(T2_scl_per_Q)/T2_scl_per_Q[0], marker='s', color=tcolor[idx], mec='k', linestyle='solid', markersize=8, label='%s %s (T2 e- scalers)'%(target[idx], kin[jdx]))


                # Plot relative T2 scalers vs (current, cumulative charge, or T2 rates)
                #plt.plot(avg_current,  unumpy.nominal_values(T2_scl_per_Q)/T2_scl_per_Q[0], marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
                #plt.plot(Qsum,  unumpy.nominal_values(T2_scl_per_Q)/T2_scl_per_Q[0], marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='dashed', label='%s %s'%(target[idx], kin[jdx]))
            #if(kin[jdx]=='SRC' and (target[idx]=='Ca40')):
            #    plt.plot(Qsum,  unumpy.nominal_values(T2_scl_per_Q)/T2_scl_per_Q[2], marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='dashed', label='%s %s'%(target[idx], kin[jdx]))
                
                #plt.errorbar(avg_current,  unumpy.nominal_values(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, unumpy.std_devs(real_yield_corr_per_Q)/real_yield_corr_per_Q[0].n, marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='', label='%s %s'%(target[idx], kin[jdx]))


                #plt.xlabel('T2 Scaler Rates [kHz]', fontsize=18)
                #plt.xlabel('Cumulative Charge [mC]', fontsize=18)
                #plt.ylabel('Relative Yield (or T2 scalers) / mC ', fontsize=18)
                #plt.ylabel('Average Current [uA]', fontsize=18)
                #plt.xticks(fontsize=14)
                #plt.yticks(fontsize=14)
                #plt.grid(True)
                 
            #if((kin[jdx]=='MF') & (target[idx]=='LD2')):
            #    plt.errorbar(run, unumpy.nominal_values(real_yield_corr_per_Q)/200., unumpy.std_devs(real_yield_corr_per_Q)/200., marker='o', markersize=8, color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
            #if((kin[jdx]=='MF') & (target[idx]!='LD2')):
            #    plt.errorbar(run, unumpy.nominal_values(real_yield_corr_per_Q)/100., unumpy.std_devs(real_yield_corr_per_Q)/100., marker='o', markersize=8, color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
            #if(kin[jdx]=='SRC'):
            #    plt.errorbar(run, unumpy.nominal_values(real_yield_corr_per_Q), unumpy.std_devs(real_yield_corr_per_Q), marker='^', markersize=8, color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
            #if((kin[jdx])=='MF'):
            #   plt.plot(run, avg_current, marker='o', markersize=8, color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
            #if((kin[jdx])=='SRC'):
            #   plt.plot(run, avg_current, marker='^', markersize=8, color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
    
                #plt.errorbar(Qsum, unumpy.nominal_values(real_yield_corr_per_Q), unumpy.std_devs(real_yield_corr_per_Q), marker='^', markersize=8, color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
                #plt.xlabel('Cumulative Charge [mC]', fontsize=18)
            #plt.xlabel('Run Number', fontsize=18)
            #plt.ylabel('Avearge Current', fontsize=18)
            #plt.ylabel('Yield / mC', fontsize=18)
            #plt.xticks(fontsize=14)
            #plt.yticks(fontsize=14)
            #plt.grid(True)
            
            #plt.tick_params(axis='both', which='major', labelsize=18)

            #plot 3: efficiencies
            #plt.errorbar(avg_current, unumpy.nominal_values(hms_trk_eff), unumpy.std_devs(hms_trk_eff), marker=kmarker[jdx], color=tcolor[idx], mec='k',  markersize=8, linestyle='None', label='%s %s '%(target[idx], kin[jdx]))
            #plt.errorbar(T3_scl_rate, unumpy.nominal_values(hms_trk_eff), unumpy.std_devs(hms_trk_eff), marker=kmarker[jdx], color=tcolor[idx], mec='k',  markersize=8, linestyle='None', label='%s %s '%(target[idx], kin[jdx]))

            #plt.errorbar(avg_current, unumpy.nominal_values(shms_trk_eff), unumpy.std_devs(shms_trk_eff), marker=kmarker[jdx], color=tcolor[idx], mec='k',  markersize=8, linestyle='None', label='%s %s '%(target[idx], kin[jdx]))
            #plt.errorbar(T2_scl_rate, unumpy.nominal_values(shms_trk_eff), unumpy.std_devs(shms_trk_eff), marker=kmarker[jdx], color=tcolor[idx], mec='k',  markersize=8, linestyle='None', label='%s %s '%(target[idx], kin[jdx]))

            #plt.errorbar(avg_current, unumpy.nominal_values(total_LT), unumpy.std_devs(total_LT), marker=kmarker[jdx], color=tcolor[idx], mec='k',  markersize=8, linestyle='None', label='%s %s '%(target[idx], kin[jdx]))
            #plt.errorbar(T2_scl_rate, unumpy.nominal_values(total_LT), unumpy.std_devs(total_LT), marker=kmarker[jdx], color=tcolor[idx], mec='k',  markersize=8, linestyle='None', label='%s %s '%(target[idx], kin[jdx]))

            #plt.ylabel('HMS Track Efficiency', fontsize=18)
            #plt.ylabel('Total Live Time Efficiency', fontsize=18)
            #plt.xlabel('Average Current [uA]', fontsize=18)
            #plt.xlabel('T3 Scaler Rate [kHz]', fontsize=18)
            #plt.xticks(fontsize=14)
            #plt.yticks(fontsize=14)
            #plt.grid(True)
                #plt.errorbar(run, unumpy.nominal_values(mult_trk_eff), unumpy.std_devs(mult_trk_eff), marker='D', color=tcolor[idx], mec='k', linestyle='None')

                
                

            #-----------------------------------------------------


            # Write numerical data to final summary file
            ofile.write("%s, %s, %.2f, %.2f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.4f, %.3f, %.1f, %.1f, %.1f\n" % (target[idx].strip(), kin[jdx].strip(), total_beam_time, total_avg_current, total_charge, real_yield_total.n, real_yield_total.s, real_Yield_corr_total.n, real_Yield_corr_total.s, yield_norm.n, yield_norm.s, tgt_areal_density, T, N, Z, A) )


            #----------------------------------------
            # APPLY Boron-Carbide TARGET CORRECTIONS
            #----------------------------------------
            # B4C10, B4C11 targets need to be carbon-subtracted
            # Ca48 needs to be corrected for impurity
            # check if final line in summary file has been reached
            if ((idx == len(target)-1) and (jdx == len(kin)-1)):

                # read final summary file
                df = pd.read_csv(ofname, comment='#')
                
                # select specific data column for carbon-subtraction
                # g/cm2
                b4c10_density = df[(df['target']=='B10')]['tgt_area_density'][:-1]
                b4c11_density = df[(df['target']=='B11')]['tgt_area_density'][:-1]
                c12_density   = df[(df['target']=='C12')]['tgt_area_density'][:-1]

                # Avogadro's number (# atoms / mol)
                Na = 6.0221408e23  

                # isotopic molar mass (g/mol)
                mol_b10 = 10.0129369  
                mol_b11 = 11.0093052
                mol_c12 = 12.0107

                # boron-carbide (4 boron-10 boron-11 atoms + 1 carbon atom) molar mass
                mol_b4c10 = 4*mol_b10 + mol_c12 
                mol_b4c11 = 4*mol_b11 + mol_c12 

                # number of target atoms / cm^2 (# of scatterers) = (atoms/mol) * g/cm2  / (g/mol)  
                targetfac_b4c10 = Na * b4c10_density  /  mol_b10
                targetfac_b4c11 = Na * b4c11_density  /  mol_b11
                targetfac_c12   = Na * c12_density    /  mol_c12

                # get the efficiency corrected yield and total charge for b4c10, b4c11, c12 (mf, src)
                yield_corr_b4c10     = unumpy.uarray(df[(df['target']=='B10')]['yield_corr'], df[(df['target']=='B10')]['yield_corr_err']) #val +/- err
                charge_b4c10         = df[(df['target']=='B10')]['total_charge']

                yield_corr_b4c11     = unumpy.uarray(df[(df['target']=='B11')]['yield_corr'], df[(df['target']=='B11')]['yield_corr_err'])
                charge_b4c11         = df[(df['target']=='B11')]['total_charge']

                yield_corr_c12       = unumpy.uarray(df[(df['target']=='C12')]['yield_corr'], df[(df['target']=='C12')]['yield_corr_err'])
                charge_c12           = df[(df['target']=='C12')]['total_charge']

                # apply carbon subtraction (counts / mC)  
                yield_norm_b410      = (yield_corr_b4c10/charge_b4c10) - (yield_corr_c12/charge_c12) * targetfac_b4c10/targetfac_c12
                yield_norm_b411      = (yield_corr_b4c11/charge_b4c11) - (yield_corr_c12/charge_c12) * targetfac_b4c11/targetfac_c12

                # apply the remaining scale factors (When boron yield is recovered, how to get the areal density (g/cm2) of ONLY boron, for normalizing?)
                yield_norm_b410   = yield_norm_b410 / ( df[(df['target']=='B10')]['T'] )

                

    
    ofile.close()

    

def make_double_ratio():


    filename='cafe_final_summary.csv'
    df = pd.read_csv(filename, comment='#')
    
    src_yield_norm = df[(df['kin']==' SRC')]['yield_norm']
    src_yield_norm_err = df[(df['kin']==' SRC')]['yield_norm_err']

    mf_yield_norm = df[(df['kin']==' MF')]['yield_norm']
    mf_yield_norm_err = df[(df['kin']==' MF')]['yield_norm_err']

    src_yield_norm_C12 = df[(df['kin']==' SRC') & (df['target']=='C12')]['yield_norm']
    src_yield_norm_err_C12 = df[(df['kin']==' SRC') & (df['target']=='C12')]['yield_norm_err']

    mf_yield_norm_C12 = df[(df['kin']==' MF') & (df['target']=='C12')]['yield_norm']
    mf_yield_norm_err_C12 = df[(df['kin']==' MF') & (df['target']=='C12')]['yield_norm_err']

    # put into arrays for error calculation
    src_yield_norm_arr  = unumpy.uarray(src_yield_norm,  src_yield_norm_err)
    mf_yield_norm_arr   = unumpy.uarray(mf_yield_norm,  mf_yield_norm_err)

    src_yield_norm_C12_arr  = unumpy.uarray(src_yield_norm_C12, src_yield_norm_err_C12)
    mf_yield_norm_C12_arr   = unumpy.uarray(mf_yield_norm_C12,  mf_yield_norm_err_C12)

    double_ratio = (src_yield_norm_arr/mf_yield_norm_arr) / (src_yield_norm_C12_arr/mf_yield_norm_C12_arr )
    N = df[(df['kin']==' SRC')]['N']
    Z = df[(df['kin']==' SRC')]['Z']
    A = df[(df['kin']==' SRC')]['A']
    
    print('double_ratio = ', double_ratio)

    plt.errorbar(N/Z, unumpy.nominal_values(double_ratio), unumpy.std_devs(double_ratio), marker='o', markersize=8, color='k', mec='k', linestyle='None', label='double ratio')
    plt.ylabel('SRC High Momentum Fraction', fontsize=18)
    plt.xlabel('N/Z', fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.show()
# Call Functions
    
make_final_summary()
make_double_ratio()
