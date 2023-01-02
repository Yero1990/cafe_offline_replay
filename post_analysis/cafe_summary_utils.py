import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy

'''
Author: C. Yero
Date: Jan 02, 2022
Brief: Compilation of scripts for 
handling/combining cafe numerical 
(.csv) summary files
'''

def find_param(param='', fname=''):

    
    with open(fname) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(param) != -1:
                print('param: ', line)
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
    fname = 'cafe_final_summary.csv' 
    ofile = open(fname, 'w')
    ofile.write('# CaFe Final Summary File (Pass 1) \n')
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# target   : target name analyzed \n'
                '# kin      : kinematics analyzed  \n'
                )
    ofile.write('target,kin,beam_time,avg_current,total_charge,yield_raw,yield_raw_err,yield_corr,yield_corr_err,yield_norm,yield_norm_err,tgt_density,T,N,Z,A\n') 

    # target, kin list
    target = ['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54']
    kin    = ['MF', 'SRC']

    for idx in np.arange(len(target)):
        
        for jdx in np.arange(len(kin)):
            
            # get generic summary file name
            summary_file_path = 'summary_files/pass1/cafe_prod_%s_%s_report_summary.csv' % (target[idx], kin[jdx])
            
            # read .csv file
            df = pd.read_csv(summary_file_path, comment='#')

            # select specific Ca48 MF runs (ignore contaminated runs)
            if( target[idx]=='Ca48' and kin[jdx]=='MF' ):
                df = df[ df['run']==17096 ] # for now, last run is good enough
            
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
            
            real_yield_total = real_Yield.sum()
            
            # apply efficiency corrections to real yield (with uncertasinties included)
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


            # Write numerical data to final summary file
            ofile.write("%s, %s, %.2f, %.2f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.4f, %.3f, %.1f, %.1f, %.1f\n" % (target[idx], kin[jdx], total_beam_time, total_avg_current, total_charge, real_yield_total.n, real_yield_total.s, real_Yield_corr_total.n, real_Yield_corr_total.s, yield_norm.n, yield_norm.s, tgt_areal_density, T, N, Z, A) )

    ofile.close()

    

    
make_final_summary()
