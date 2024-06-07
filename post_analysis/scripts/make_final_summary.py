import numpy as np
import pandas as pd
from itertools import takewhile
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy

'''
Author: C. Yero
Date: April 04, 2023
Brief: Compilation of scripts for 
handling/combining cafe numerical 
(.csv) summary files for plotting
double ratios
'''

# user need to specify which pass to analyze
if (len(sys.argv) != 2):
    print('''\n
    ---------------------------------------------
    Usage: ipython make_final_summary.py <passN>
    <passN>: data pass number to analyze
    (pass1, pass2, ...)
    
    e.g. ipython make_final_summary.py pass3
    ---------------------------------------------
    ''')
    sys.exit()

npass = sys.argv[1] 

#________________________
def get_cntmEff(run=0, tgt=''):
    # method: get Ca48 contamination efficinecy factor on a run-basis
    # tgt = 'ca40' or 'ca48'
    # set filename
    fname='special_studies/ca_contamination/%s_correction_noah.csv'%(tgt)

    # read .csv file
    df = pd.read_csv(fname, comment='#')

    # define condition
    cond = df['run']==run

    # get absolute contamination (in percent) -- to use this, requires that there is not space between comma-sepatated values
    absContamCalc = df[cond].C_absCntm_calc    # calculated from fit
    #absContamMeas = df[cond].C_absCntm_meas    # measured (only for Ca48 MF)

    contam_eff_calc = 1. -  (absContamCalc/100.)
    #contam_eff_meas = 1. -  (absContamMeas/100.)

    # it depends on whether we want to apply calculated or measured contamination eff. for the Ca48 MF
    # for now we use the fit results 
    return contam_eff_calc  # efficiency factor

#________________________
def get_cntmEff(kin='', tgt='', get_err=False):
    # method: returns an array for contamination factor for either Ca48 MF or SRC
    # tgt = 'ca40' or 'ca48'
    # example of appending new row (to append corrected Ca48 yield)
    fname='special_studies/ca_contamination/%s_correction_noah.csv'%(tgt)

    # read .csv file
    df = pd.read_csv(fname, comment='#')

    # define condition (either MF or SRC)
    cond = df['kin']==kin

    # get absolute contamination (in percent) -- to use this, requires that there is not space between comma-sepatated values
    absContamCalc = df[cond].C_absCntm_calc    # calculated from fit (exponential fit func. of relative scalers for Ca48 SRC is used to quantify contamination)
    absContamCalc_err = df[cond].C_absCntm_calc_err # absolute uncertainty of contamination (in percent)
    
    #absContamMeas = df[cond].C_absCntm_meas    # measured (only for Ca48 MF)

    contam_eff_calc = 1. -  (absContamCalc/100.)
    contam_eff_calc_err = absContamCalc_err/100.  # absolute error in correction factor (apply error propagation)

    
    #contam_eff_meas = 1. -  (absContamMeas/100.)

    # it depends on whether we want to apply calculated or measured contamination eff. for the Ca48 MF
    # for now we use the fit results 

    if get_err:
        return np.array(contam_eff_calc_err)  # return abs. error in eff. factor
    else:
        return np.array(contam_eff_calc)  # efficiency factor


#_____________________________________
def find_param(param='', fname=''):
    #brief: function that read parameters from .csv file
    
    with open(fname) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.find(param) != -1:
                #print('param: ', line)
                param_value = float(line.split(':')[1].strip())
                return param_value

#_____________________________________
def make_final_summary():

    #brief: function that reads summary files (.csv) to:
    # 1. apply efficiency corrections on a run-by-run basis
    #    eff. corrections: hms/shms track eff, shms multi-track eff, edtm total live time,  
    # 2. sum over all efficiency corrected counts 
    # 3. normalize  efficiency-corrected counts by charge (mC), transparency and target density (g/cm^2)
    # 4. write combined numerical quantitied to a .csv file
    #  ( for additional subtractions -> boron-carbide subtractions, Ca48 impurity and contamination corrections, and double ratio calculations)

    #C.Y. definition changed:  yield_norm --> sigma_rawg : yield_corr normalized by (charge * g/cm^2 * transparency)
    # sigma_raw_A -->  [area_density(g/cm^2) / molar_mass(g/mol)] * 4 (B10 atoms) * 10 g/mol ??
    
    #output file to write summary file
    ofname = 'cafe_final_summary_%s_testing.csv' %(npass) 
    ofile = open(ofname, 'w+')
    ofile.write('# CaFe Final Summary File (%s) \n'%(npass))
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# target       : target name analyzed \n'
                '# kin          : kinematics analyzed  \n'
                '# beam_time    : beam-on-target time [s] \n'
                '# avg_current  : average beam current [uA] \n'
                '# total_charge : cumulative charge (over all runs) [mC] \n'
                '# yield        : yield (raw counts integrated over Pm) with all data-analysis cuts applied \n'
                '# yield_eff    : yield corrected for inefficiencies (hms/shms tracking + live_time + proton absorption) \n'
                '# yield_corr   : yield_eff corrected for external/internal impurities if any \n'
                '# rad_corr     : radiative correction factor, yield_rad / yield_norad (to be applied when calculating ratios) \n'
                '# sigma_raw_per_nucleon   : corrected yield (yield_corr) normalized by (total_charge|transpacency|area_density(g/cm2) (per nucleon))\n'
                '# sigma_raw_per_proton    : corrected yield (yield_corr) normalized by (total_charge|transpacency|area_density (per proton) e.g. xA/Z )\n'
                '# sigma_raw_per_nucleus   : corrected yield (yield_corr) normalized by (total_charge|transpacency|area_density (per nucleus) e.g xA )\n'
                '# stat_rel_err            : relative statistical error on the total yield calculated as sqrt(yield)/yield \n'
                '# norm_syst_rel_err       : relative systematic error due to normalization correction factors (hms/shms track efficiencies, live time, proton transmission, ca40,48 cntm) \n'
                '# tgt_thick               : target density (g/cm2)\n'
                '# T N Z A                 : transparency (T) # of neutrons (N) protons(Z) and nucleons (A)\n'
                )
    ofile.write('target,kin,beam_time,avg_current,total_charge,yield,yield_eff,yield_corr,rad_corr,sigma_raw_per_nucleon,sigma_raw_per_proton,sigma_raw_per_nucleus,stat_rel_err,norm_syst_rel_err,tgt_thick,tgt_thick_corr,T,N,Z,A\n') 

    # target, kin list
    # target = ['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54']
    target = ['Be9', 'B10', 'B11', 'C12', 'Ca40',       'Ca48',   'Fe54', 'Au197']
    tcolor  = ['m',   'r',   'g',   'b',  'darkorange', 'dodgerblue', 'lime',  'gold'] 
    
    kin    = ['MF', 'SRC']

    # Set pre-defined varaibles to be set and  saved for corrections
    Q_Ca40 = {'MF': 0, 'SRC':0}
    real_Yield_corr_total_Ca40  = {'MF': 0, 'SRC':0}
    ca40_density = 0
    
    # ------ create subplots for quality check plotting -----------

    # efficiencies vesrus T2 scaler rate
    fig1, ax1 = plt.subplots(nrows=2, ncols=3)
    fig1.set_size_inches(14,8, forward=True)

    # T2 scaler / ( Q * tgt_thick )
    fig2, ax2 = plt.subplots(nrows=1, ncols=1)
    fig2.set_size_inches(12,8, forward=True)

    # sigma_raw_per_nucleon vs. run number (on a run-by-run basis) for sanity check
    fig3, ax3 = plt.subplots(nrows=1, ncols=1)
    fig3.set_size_inches(8,8, forward=True)
    
    # -------------------------------------------------------------
    # loop over each target
    for idx in np.arange(len(target)):

        # loop over each kinmeatic (MF, SRC)
        for jdx in np.arange(len(kin)):

            print('target ->', target[idx], ' kin -> ', kin[jdx])

            rad_corr = -1000  # initiate default radiative correction factor  Y_rad / Y_norad (usually < 1)
            
            # set radiative correction factors (based on studies by N. Swan)
            if(kin[jdx]=='MF'):
                if(target[idx]=='Be9' or target[idx]=='B10' or target[idx]=='B11' or target[idx]=='C12'):
                    rad_corr = 0.6
                if(target[idx]=='Ca40' or target[idx]=='Ca48' or target[idx]=='Fe54'):
                    rad_corr = 0.51
                if(target[idx]=='Au197'):
                    rad_corr = 0.38
            
            if(kin[jdx]=='SRC'):
                if(target[idx]=='Be9' or target[idx]=='B10' or target[idx]=='B11' or target[idx]=='C12'):
                    rad_corr = 0.63
                if(target[idx]=='Ca40' or target[idx]=='Ca48' or target[idx]=='Fe54'):
                    rad_corr = 0.52
                if(target[idx]=='Au197'):
                    rad_corr = 0.40

            # set generic summary file name
            summary_file_path = 'summary_files/%s/cafe_prod_%s_%s_report_summary.csv' % (npass,target[idx], kin[jdx])
            
            # read .csv file
            df = pd.read_csv(summary_file_path, comment='#')

  
            # ------- check if target is Ca48 (and read contamination factors) --------

            if(target[idx]=='Ca48'):
                
                # read array of contamination eff. factors
                cntm_eff = get_cntmEff(kin[jdx], 'ca48')
                cntm_eff_err = get_cntmEff(kin[jdx], 'ca48', True)
                
                if (kin[jdx]=='MF'):
                    cntm_eff = cntm_eff[-3:] # read contamination factor corresponding to last 3 Ca48 MF runs
                    cntm_eff_err = cntm_eff_err[-3:]
                    df = df[-3:]  # select last 3 Ca48 MF runs from dataframe (esentially almost no contamination)

            # ----- END: check if target is Ca48 (and read contamination factors) --------

            # ------- check if target is Ca40 (and read contamination factors) --------
            
            if(target[idx]=='Ca40'):
                 
                # read array of contamination eff. factors
                cntm_eff = get_cntmEff(kin[jdx], 'ca40')
                cntm_eff_err = get_cntmEff(kin[jdx], 'ca40', True)
            
                     
            
            # ------ read parameters -------
            #T                 = find_param('transparency', summary_file_path) # transparency
            tgt_thick         = find_param('target_areal_density', summary_file_path) # g/cm^2
            tgt_thick_corr    = tgt_thick  # set corrected thickness to thickness (will be re-defined if impurity is corrected for any target)
            N                 = find_param('N:', summary_file_path) # number of neutrons
            Z                 = find_param('Z:', summary_file_path) # number of protons
            A                 = find_param('A:', summary_file_path) # number of nucleons
            # ---- END: read parameters -----

            # Transparency function: T = c * A ** alpha (Q2), where alpha ~ -0.24 for Q2 >= 2 GeV^2, and c=1, A -> mass number
            # reference: https://arxiv.org/abs/1211.2826  "Color Transparency: past, present and future"
            #alpha=-0.24
            alpha = -1/3.   # May 30, 2024 : LW suggests to use this factor from a previously published paper
            T = A**(alpha)
            
            # read selected data columns (with respective uncertainties)
            avg_current  = df['avg_current'] # average beam current [uA]     
            run          = df['run']
            charge       = df['charge'] # [mC]
            beam_time    = df['beam_time'] # (beam-on-target time) [sec]

            '''
            real_Yield   = unumpy.uarray(df['real_Yield'],      df['real_Yield_err']) 
            hms_trk_eff  = unumpy.uarray(df['hTrkEff'],         df['hTrkEff_err'])
            shms_trk_eff = unumpy.uarray(df['pTrkEff'],         df['pTrkEff_err'])
            total_LT     = unumpy.uarray(df['tLT'],             df['tLT_err_Bi'])
            #mult_trk_eff = unumpy.uarray(df['multi_track_eff'], df['multi_track_eff_err'])
            mult_trk_eff = np.array(df['multi_track_eff'])  # need to figure way to calculate uncertainty on this variable
            '''
            real_Yield   = df['real_Yield']  
            real_Yield_stat_err   = df['real_Yield_err']  # absolute statistical error on raw counts, calculated as sqrt(counts)   

            hms_trk_eff  = df['hTrkEff']           
            hms_trk_eff_err = df['hTrkEff_err']

            shms_trk_eff = df['pTrkEff']
            shms_trk_eff_err = df['pTrkEff_err']

            total_LT     = df['tLT']
            total_LT_err = df['tLT_err_Bi']

            mult_trk_eff = df['multi_track_eff']  # df['multi_track_eff_err']
            mult_trk_eff_err = df['multi_track_eff_err']  # need to figure out how to fix error calculation on multi-track eff. (as it is too large)

            pTrms_eff = 0.952 # Noah Swan's proton transmission factor
            pTrms_eff_err = 0.004
            
            T1_scl_rate  = df['T1_scl_rate']  # SHMS  3/4    kHz
            T2_scl_rate  = df['T2_scl_rate']  # SHMS EL-REAL kHz
            T3_scl_rate  = df['T3_scl_rate']  # HMS 3/4      kHz
            T5_scl_rate  = df['T5_scl_rate']  # COIN         kHz

            
            # calculate total (cumulative) charge (i.e., array with incremental charge)
            Qsum = charge.cumsum()

            # sum over all total charge
            total_charge = charge.sum()

            # total beam time
            total_beam_time = beam_time.sum()

            # total average current
            total_avg_current = np.average(avg_current)
            
            # calculate average hms/shms track efficiency and live time
            # (these may be useful later on, for normalizing yield)
            avg_hms_trk_eff  =  hms_trk_eff.mean()
            avg_shms_trk_eff =  shms_trk_eff.mean()
            avg_total_LT     =  total_LT.mean()
            avg_mult_trk_eff =  mult_trk_eff.mean()

            #------------------------------------------------------------------
            # Define and Calculate the Normalization Systematic Uncertainties
            #------------------------------------------------------------------
            # Y_corr = Y_uncorrr * f1 * f2 * . . . fn,   --> apply correction factors to raw yield (each fn has an uncertainty)
            # where fn = 1/corr_factor, and corr_factor = hms_trk_eff, shms_trk_eff, total_LT, multi_trk_eff, total_charge, proton_transmission
            # (dY_corr / Y_corr)_syst2 = (df1/f1)**2 + (df2/f2)**2 +  . . . (dfn/fn)**2 --> relative norm. syst is sum in quadrature of relative error in corr. factors


            
            # *take the average of all runs rather than run-by-run
            # calculate relative systematic uncertainty on hms_trk_eff
            w_avg = -1
            w_avg_err = -1
            weights = -1
            
            #--------------------------------
            
            # calculate relative systematic uncertainty on hms_trk_eff
            weight =  1. /  hms_trk_eff_err**2
            w_avg = np.average(hms_trk_eff, weights=weight)
            w_avg_err = np.sqrt(1. / np.sum(weight))
            
            f1 = 1./w_avg
            df1 = (1./w_avg**2) * w_avg_err
            df1_f1_sq = np.square(df1 / f1)
            print('df1_f1_sq (shms_trk_eff) = ', df1_f1_sq)
             
            # calculate relative systematic uncertainty on shms_trk_eff
            weight =  1. /  shms_trk_eff_err**2
            w_avg = np.average(shms_trk_eff, weights=weight)
            w_avg_err = np.sqrt(1. / np.sum(weight))

            f2 = 1./w_avg
            df2 = (1./w_avg**2) * w_avg_err
            df2_f2_sq = np.square(df2 / f2)
            print('df2_f2_sq (shms_trk_eff) = ', df2_f2_sq)        
     
            # calculate relative systematic uncertainty on total_LT
            weight =  1. /  total_LT_err**2
            w_avg = np.average(total_LT, weights=weight)
            w_avg_err = np.sqrt(1. / np.sum(weight))
            
            f3 = 1./w_avg
            df3 = (1./w_avg**2) * w_avg_err
            df3_f3_sq = np.square(df3 / f3)
            print('df3_f3_sq (total_LT) = ', df3_f3_sq)

            # calculate relative systematic uncertainty on multi-track eff
            weight =  1. /  mult_trk_eff_err**2
            w_avg = np.average(mult_trk_eff, weights=weight)
            w_avg_err = np.sqrt(1. / np.sum(weight))
            
            f4 = 1./w_avg
            df4 = (1./w_avg**2) * w_avg_err
            df4_f4_sq = np.square(df4 / f4)
            print('df4_f4_sq ( mult_trk_eff) = ', df4_f4_sq)

            # calculate relative systematic uncertainty on proton transmission factor (overall correction, not run-by-run)
            f5 = 1./pTrms_eff  # obtained from Noah Swan's studies of proton absorption
            df5 = pTrms_eff_err 
            df5_f5_sq = 0.0  #(df5 / f5)**2  ----> we should not include this systematic in sing/double ratio as it cancels out
            print('df5_f5_sq (proton_trms_rel_syst) = ', df5_f5_sq)

            # calculate overall relative norm. systematic unc on the relative yield
            norm_syst_rel_err = np.sqrt( df1_f1_sq + df2_f2_sq + df3_f3_sq + df4_f4_sq + df5_f5_sq)

            # --- apply Ca40 contamination systematic correction (the correction is pre-determined for either MF or SRC during get_cntmEff() function call
            if(target[idx]=='Ca40'):
                f6_per_run = cntm_eff
                df6_per_run = cntm_eff_err
                df6_f6_per_run = df6_per_run / f6_per_run
                df6_f6_sq = np.sum(np.square(df6_f6_per_run))  
                print('df6_f6_sq (ca40_cntm_rel_syst) = ', df6_f6_sq)

                norm_syst_rel_err = np.sqrt( df1_f1_sq + df2_f2_sq + df3_f3_sq + df4_f4_sq + df5_f5_sq + df6_f6_sq)

            # --- apply Ca48 contamination systematic correction (the correction is pre-determined for either MF or SRC during get_cntmEff() function call
            if(target[idx]=='Ca48'):
                f7_per_run = cntm_eff
                df7_per_run = cntm_eff_err
                df7_f7_per_run = df7_per_run / f7_per_run
                df7_f7_sq = np.sum(np.square(df7_f7_per_run))  
                print('df7_f7_sq (ca48_cntm_rel_syst) = ', df7_f7_sq)

                norm_syst_rel_err = np.sqrt( df1_f1_sq + df2_f2_sq + df3_f3_sq + df4_f4_sq + df5_f5_sq + df7_f7_sq)

            # calculate total scaler counts
            T1_scl = T1_scl_rate * 1000 * beam_time   
            T2_scl = T2_scl_rate * 1000 * beam_time
            T3_scl = T3_scl_rate * 1000 * beam_time
            T5_scl = T5_scl_rate * 1000 * beam_time 

            # calculate normalized scalers (for sanity check)
            T2_scl_norm = T2_scl / (charge *  tgt_thick_corr)
            T3_scl_norm = T3_scl / (charge *  tgt_thick_corr)
            
            
            # sum over  counts (before applying any corrections, inefficinecy, charge, etc.) for book-keeping 
            real_Yield_total = real_Yield.sum()
            total_stat_err = np.sqrt(real_Yield_total)  # total abs. stats error for a given (target, kin) setting 
            stat_rel_err = total_stat_err / real_Yield_total  # relative statistical error
            
            # apply efficiency corrections to real yield  (run-by-run) and then sum over all counts
            #real_Yield_eff       = real_Yield / (hms_trk_eff * shms_trk_eff * total_LT * mult_trk_eff)     # array of runs         
            real_Yield_eff       = real_Yield * f1 * f2 * f3 * f4
            #print('f1: 1/ hms_trk_eff_avg : hms_htrk_eff_avg = %.3f' % (1./f1) )
            real_Yield_eff_total = real_Yield_eff.sum() * f5                                                # sum of array elements (and apply overall p-transmission factor)

            # re-define efficiency yield as corrected yield  to be used to  correct for any target impurities (either internal contamination, or external)
            # some targets may not need any correction and so real_Yield_eff = real_Yield_corr if no impurities found
            
            real_Yield_corr       = real_Yield_eff
            real_Yield_corr_total = real_Yield_eff_total

            # define raw cross section: corrected yield normalized by  total charge, transparency and target density (g/cm2),
            # for now, do not include systematic error due to (charge, Transparency and tgt_thickness)
            sigma_raw_per_nucleon_per_run = real_Yield_corr / (charge * T * tgt_thick)  # counts / (mC * g/cm^2)
            sigma_raw_per_nucleon =  real_Yield_corr_total / (total_charge * T * tgt_thick)  # counts / (mC * g/cm^2)
            sigma_raw_per_nucleus   =   sigma_raw_per_nucleon * A  # counts / (mC * T * g/cm^2 * nucleus)
            sigma_raw_per_proton   =   sigma_raw_per_nucleon * A / Z  # counts / (mC * T * g/cm^2 * proton)

            
            # --- apply Ca40 contamination
            if(target[idx]=='Ca40'):

                # --- Ca40 oil contamination correction ---
                real_Yield_corr = real_Yield_corr * cntm_eff  # apply oil contamination run-by-run           
                real_Yield_corr_total = real_Yield_corr.sum() # sum the oil-corrected runs 
                
                # define raw cross section: corrected yield normalized by  total charge, transparency 
                sigma_raw_per_nucleon_per_run = real_Yield_corr / (charge * T * tgt_thick)       # counts / (mC * g/cm^2)
                sigma_raw_per_nucleon =  real_Yield_corr_total / (total_charge * T * tgt_thick)  # counts / (mC * g/cm^2)
                sigma_raw_per_nucleus   =  sigma_raw_per_nucleon * A
                sigma_raw_per_proton    =  sigma_raw_per_nucleon * A / Z
            
            # --- apply Ca48 contamination + impurities corrections ----
            if(target[idx]=='Ca48'):

                # --- Ca48 oil contamination correction ---
                print('cntm_eff ---------> ', cntm_eff)
                real_Yield_corr = real_Yield_corr * cntm_eff  # apply oil contamination run-by-run           
                real_Yield_corr_total = real_Yield_corr.sum() # sum the oil-corrected runs 
                

                # --- Ca48 "Ca40 impurity" correction ---
                # brief: Ca48 is only ~90.5 % pure, the remaining is Ca40. The 48Ca target is 1051 mg/cm2.
                # If the purity refers to number of atoms, then the purity by weight is 91.5%.
                # (We need to find out whether purity refers to by mass or weight - ask D. Meekins.)  

                ca48_purity = 0.915   # (91.5% purity, but need to determine if its by mass or weight)

                ca40_cntm         =  (1. - ca48_purity) * tgt_thick # calculate amount of ca40 in ca48 (g/cm2)
                ca48_density_corr =  tgt_thick - ca40_cntm          # corrected ca48 target thickness (g/cm2)
                real_Yield_corr_total = real_Yield_corr_total - real_Yield_corr_total_Ca40[kin[jdx]] * (ca40_cntm / ca40_density) * (total_charge/Q_Ca40[kin[jdx]])

                # define raw cross section: corrected yield normalized by  total charge, transparency and "CORRECTED" target density (g/cm2),
                sigma_raw_per_nucleon_per_run = real_Yield_corr / (charge * T * ca48_density_corr)       # counts / (mC * g/cm^2)
                sigma_raw_per_nucleon =  real_Yield_corr_total / (total_charge * T * ca48_density_corr)  # counts / (mC * g/cm^2)
                sigma_raw_per_nucleus   =  sigma_raw_per_nucleon * A
                sigma_raw_per_proton    =  sigma_raw_per_nucleon * A / Z
                tgt_thick_corr    = ca48_density_corr               # re-define corrected target thickness (to be written to file)
  
            # --- apply B4C10, B4C11 carbon-subtraction ----
            # brief:  B4C10, B4C11 targets need to be carbon-subtracted (4 Boron-10 atoms + 1 C12 atom)
            if(target[idx]=='B10' or target[idx]=='B11'):

                # get necessary C12 information 
                c12_density = find_param('target_areal_density', 'summary_files/%s/cafe_prod_C12_MF_report_summary.csv'%(npass)) #g/cm2

                # define c12 dataframe to get charge and yield from C12 to subtract from B4C-10,11
                df_c12 = pd.read_csv('summary_files/%s/cafe_prod_C12_%s_report_summary.csv'%(npass, kin[jdx]), comment='#') 
                c12_charge       = df_c12['charge'].sum()
                c12_charge_per_run       = df_c12['charge']

                real_Yield_c12       = df_c12['real_Yield']
                real_Yield_c12_err   = df_c12['real_Yield_err']

                hms_trk_eff_c12      = df_c12['hTrkEff']
                hms_trk_eff_c12_err  = df_c12['hTrkEff_err']

                shms_trk_eff_c12     = df_c12['pTrkEff']
                shms_trk_eff_c12_err = df_c12['pTrkEff_err']

                total_LT_c12         = df_c12['tLT']
                total_LT_c12_err     = df_c12['tLT_err_Bi']

                mult_trk_eff_c12     = df_c12['multi_track_eff']
                mult_trk_eff_c12_err = df_c12['multi_track_eff_err']

                real_Yield_corr_c12 = real_Yield_c12 / (hms_trk_eff_c12*shms_trk_eff_c12*total_LT_c12*mult_trk_eff_c12*pTrms_eff) 
                real_Yield_corr_total_c12 = real_Yield_corr_c12.sum()
                
                # calculate the c12 density contribution in b4c (subtracted the boron contribution from the b4c density)
                boron_density_corr = ( 4.*A / (4.*A + 12) ) * tgt_thick   # boron density corrected = b4c density scaled by (nucleons in b4)/ (nucleons in b4c)
                c12_density_b4c    = tgt_thick -  boron_density_corr   # c12 density contribution in b4c 

                # subtract c12 contribution from b4c yield
                real_Yield_corr_total = real_Yield_corr_total - real_Yield_corr_total_c12 * (c12_density_b4c/c12_density) * (total_charge/c12_charge)  # corrected boron yield  (total)

                # define raw cross section: corrected yield normalized by  total charge, transparency and "CORRECTED" target density (g/cm2),
                sigma_raw_per_nucleon      =  real_Yield_corr_total /  (total_charge * T * boron_density_corr)  # counts / (mC * g/cm^2)
                sigma_raw_per_nucleus      =  sigma_raw_per_nucleon * A
                sigma_raw_per_proton       =  sigma_raw_per_nucleon * A / Z
                tgt_thick_corr             = boron_density_corr               # re-define corrected target thickness (to be written to file)            


            show_plots=False
            # ----------- make plots ----------------
            minT2 = T2_scl_rate==min(T2_scl_rate) #condition of minimum scaler rate
            minI  = avg_current==min(avg_current)
            if(show_plots and kin[jdx]=="SRC"):
                fig1.suptitle('%s Kinematics'%(kin[jdx]), fontsize=20)

                ax1[0, 0].errorbar(T2_scl_rate, unumpy.nominal_values(shms_trk_eff),  yerr=unumpy.std_devs(shms_trk_eff), marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1[0, 0].set_title('SHMS Track Eff.',fontsize=14)
                ax1[0, 0].set_xlabel('SHMS (T2) Scaler Rate [kHz]',fontsize=14)
                ax1[0, 0].set_ylabel(r'$\epsilon_{trk,SHMS}$', fontsize=16)
                #y_loc = 0.985 - idx/500
                #plt.text(150, y_loc, '%s'%(target[idx]), color=tcolor[idx], fontsize = 15)

                ax1[0, 1].plot(T2_scl_rate, mult_trk_eff, marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1[0, 1].set_title('SHMS Multi-Track Eff.',fontsize=14)
                ax1[0, 1].set_xlabel('SHMS (T2) Scaler Rate [kHz]', fontsize=14)
                ax1[0, 1].set_ylabel(r'$\epsilon_{multi,SHMS}$',fontsize=16)
                
                ax1[0, 2].errorbar(T3_scl_rate, unumpy.nominal_values(hms_trk_eff),  yerr=unumpy.std_devs(hms_trk_eff), marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1[0, 2].set_title('HMS Track Eff.',fontsize=14)
                ax1[0, 2].set_xlabel('HMS (T3) Scaler Rate [kHz]',fontsize=14)
                ax1[0, 2].set_ylabel(r'$\epsilon_{trk,HMS}$', fontsize=16)
                
                ax1[1, 0].errorbar(T2_scl_rate, unumpy.nominal_values(total_LT),  yerr=unumpy.std_devs(total_LT), marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1[1, 0].set_title('Total EDTM Live Time',fontsize=14)
                ax1[1, 0].set_xlabel('SHMS (T2) Scaler Rate [kHz]',fontsize=14)
                ax1[1, 0].set_ylabel(r'$\epsilon_{tLT}$',fontsize=16)
                ax1[1, 0].legend(target,loc="lower right")
                
                # plot charge-norm yield vs. T2 rates
                ax1[1, 1].errorbar(T2_scl_rate, unumpy.nominal_values( sigma_raw_per_nucleon_per_run)/ unumpy.nominal_values( sigma_raw_per_nucleon_per_run)[minT2],  yerr=unumpy.std_devs( sigma_raw_per_nucleon_per_run)/ unumpy.nominal_values( sigma_raw_per_nucleon_per_run)[minT2], marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                #ax1[1, 1].errorbar(T2_scl_rate, unumpy.nominal_values( sigma_raw_per_nucleon_per_run),  yerr=unumpy.std_devs( sigma_raw_per_nucleon_per_run), marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1[1, 1].set_title('Relative Charge-Normalized Yield', fontsize=14)
                ax1[1, 1].set_xlabel('SHMS (T2) Scaler Rate [kHz]', fontsize=14)
                ax1[1, 1].set_ylabel(r'$Y/Y_{T2,min}$', fontsize=16)
                
                #ax1[1, 2].errorbar(T2_scl_rate, unumpy.nominal_values( sigma_raw_per_nucleon_per_run)/ unumpy.nominal_values( sigma_raw_per_nucleon_per_run)[minT2],  yerr=unumpy.std_devs( sigma_raw_per_nucleon_per_run)/ unumpy.nominal_values( sigma_raw_per_nucleon_per_run)[minT2], marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1[1, 2].errorbar(avg_current, unumpy.nominal_values( sigma_raw_per_nucleon_per_run)/ unumpy.nominal_values( sigma_raw_per_nucleon_per_run)[minI],  yerr=unumpy.std_devs( sigma_raw_per_nucleon_per_run)/ unumpy.nominal_values( sigma_raw_per_nucleon_per_run)[minI], marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1[1, 2].set_title('Relative Charge-Normalized Yield', fontsize=14)
                ax1[1, 2].set_xlabel(r'Average Current [$\mu$A]', fontsize=14)
                ax1[1, 2].set_ylabel(r'$Y/Y_{I,min}$', fontsize=16)
                

                # ---- PLOT T2 (or T3) scsalers / (Q*tgt_thick) -------
                A_lst = np.array([A] * len(run))
                N_lst = np.array([N] * len(run))
                Z_lst = np.array([Z] * len(run))
                NoZ = N_lst /   Z_lst
                
                ax2.errorbar(A_lst, unumpy.nominal_values(T2_scl_norm),  yerr=unumpy.std_devs(T2_scl_norm), marker='o', markersize=10, mec='k', color=tcolor[idx], linestyle='None')
                ax2.legend(target,loc="lower right", fontsize=16)
                ax2.set_title(r'SHMS (%s) T2 Scalers / (Q $\times \sigma_{thick}$) vs A'%(kin[jdx]),fontsize=16)
                ax2.set_xlabel('Mass Number, A',fontsize=16)
                ax2.set_ylabel(r'T2 / (Q $\times$ $\sigma_{thick}$)', fontsize=14)
                ax2.set_xscale('log')

                '''
                ax2[1].errorbar(T3_scl_rate, unumpy.nominal_values(T3_scl_norm),  yerr=unumpy.std_devs(T3_scl_norm), marker='o', markersize=10, mec='k', color=tcolor[idx], linestyle='None')
                ax2[1].legend(target,loc="lower right", fontsize=16)
                ax2[1].set_title(r'SHMS T3 Scalers / (Q $\times \sigma_{thick}$)',fontsize=16)
                ax2[1].set_xlabel('T3 Scaler Rate [kHz]',fontsize=16)
                ax2[1].set_ylabel(r'T3 / (Q $\times$ $\sigma_{thick}$)', fontsize=14)
                '''
                #---- Plot sigma_raw_per_proton_per_run versus A ------
                sigma_raw_per_proton_per_run = sigma_raw_per_nucleon_per_run * A / Z
                ax3.errorbar(NoZ , unumpy.nominal_values(sigma_raw_per_proton_per_run),  yerr=unumpy.std_devs(sigma_raw_per_proton_per_run), marker='o', markersize=8, mec='k', color=tcolor[idx], linestyle='None')
                ax3.legend(target,loc="upper right", fontsize=14)
                ax3.set_title(r'raw %s $\sigma_{xsec}$ (per proton)'%(kin[jdx]),fontsize=16)
                ax3.set_xlabel('N/Z',fontsize=16)
                ax3.set_ylabel(r'Y / $(Q \cdot \epsilon \cdot \sigma_{thick}\cdot T )\times A/Z$', fontsize=14)
                #ax3.set_xscale('log')
                
                
            # SAVE SPECIFIC Ca40 TARGET INFO TO BE USED IN CORRECTIONS IN Ca48
            if(target[idx]=='Ca40'):
                Q_Ca40[kin[jdx]] = total_charge
                real_Yield_corr_total_Ca40[kin[jdx]] = real_Yield_corr_total
                ca40_density = tgt_thick
                
            #-------------------------------------------------
            # WRITE CAFE NUMERICAL DATA (for cross-checking)
            #-------------------------------------------------

            print('----------------')
            print('target: %s %s' % (target[idx].strip(), kin[jdx].strip()) )
            print('----------------')
            print('')
            print('---- raw yield (integrated over Pm) cuts: ALL, corrections: NONE ---')
            print('Y_per_run: ', real_Yield)
            print('Y_sum = SUM_[Y_per_run]: %.1f' % real_Yield_total)
            print('')
            print('---- efficiency corrections ---')
            print('NOTE: the averaged efficiency over all runs is applied on a run-by-run basis')
            print('** The proton transmission factor gets applied as an overall correction factor (after summation of Y_per_run)')
            print('Y_eff_per_run = Y_per_run / (htrk_eff_avg * etrk_eff_avg * tLT_avg * emult_trk_eff_avg) ')
            print('Y_eff = SUM_[ Y_per_run / (htrk_eff_avg * etrk_eff_avg * tLT_avg * emult_trk_eff_avg) ] * (1 / pTransm) = SUM [Y_per_run / %.3f] * 1/%.3f' % ( ((1./f1) * (1./f2) * (1./f3) * (1./f4)), 1./f5 ) )
            print('Y_eff = %.1f' % (real_Yield_eff_total))
            print('')
            print('htrk_eff_avg     : %.3f'%(1/f1))
            print('etrk_eff_avg     : %.3f'%(1/f2))
            print('tLT_avg          : %.3f'%(1/f3))
            print('emult_trk_eff_avg: %.3f'%(1/f4))
            print('pTransm          : %.3f'%(1/f5))
            print('eff_total        : %.3f'%((1./f1) * (1./f2) * (1./f3) * (1./f4) * (1./f5)))
            print('')
            print('--- additional corrections (if needed) for external/internal target impurities ---')
            if(target[idx]=='Ca40'):
                print('Ca-40 Oil Contamination efficiency ')
                print('ca40_cntm_eff_per_run = ', cntm_eff  )
                print('Y_corr_per_run = Y_eff_per_run * ca40_cntm_eff_per_run  ')
                print('Y_corr = SUM_[Y_corr_per_run] * 1/pTransm = %.3f * 1./%.3f' % (real_Yield_corr_total, (1/f5) )  )
                print('Y_corr = %.3f' % (real_Yield_corr_total * f5))
            print('')
            print('----------------')
            # Write numerical data to final summary file
            ofile.write("%s,%s,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.3f,%.1f,%.1f,%.1f\n" % (target[idx].strip(), kin[jdx].strip(), total_beam_time, total_avg_current, total_charge, real_Yield_total, real_Yield_eff_total, real_Yield_corr_total, rad_corr, sigma_raw_per_nucleon, sigma_raw_per_proton, sigma_raw_per_nucleus, stat_rel_err, norm_syst_rel_err, tgt_thick, tgt_thick_corr, T, N, Z, A) )

    if(show_plots):         
        plt.tight_layout()
        plt.show()          
    ofile.close()

    # call function to write singles, double ratios to file for ease of plotting
    write_ratios(ofname, 'cafe_ratios_%s.csv'%(npass))
    

def write_ratios(ifname='', ofname=''):

    # method: reads input file ifname with (cafe final summary file),
    #         write output file ofname with double-ratio numerical value
    #
    #         single_ratios => A1_SRC / A2_SRC,  A1_MF / A2_MF, A_SRC / A_MF
    #         double_ratio => ( A_SRC / A_MF ) / ( C12_SRC / C12_MF )
    
    
    # read input final summary file
    df = pd.read_csv(ifname, comment='#')

    # set output file to write single + double ratio numerical values
    ofile = open(ofname, 'w+')
    ofile.write('# CaFe numerical ratios (%s) \n'%(npass))
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# target       : target A used in single or double ratio \n'
                '# singleR_A_c12_mf      : single ratio of target A_MF / C12_MF (per proton)  \n'
                '# singleR_A_c12_src      : single ratio of target A_SRC / C12_SRC (per proton) \n'
                '# singleR_per_proton      : single ratio of target A(SRC/MF) (per proton) \n'
                '# doubleR      : double ratio of target A(SRC/MF) relative to C12 (SRC/MF) \n'
                '# doubleR_Jmodel : Justin/Andrew model of double ratio of target A(SRC/MF) relative to C12 (SRC/MF) \n'
                '# _av18 and _osu : single and double ratio (AV18 and OSU) models \n'
                '# N: Z: A      : # of neutrons (N): protons(Z): nucleons (A): for target A \n'
                '# NoZ          : N/Z \n'
                '# NmZoA        : (N-Z)/A \n'                
                )
    ofile.write('target,singleR_A_c12_mf,singleR_A_c12_mf_stat_err,singleR_A_c12_mf_norm_syst_err,singleR_A_c12_mf_RC_syst_err,singleR_A_c12_mf_cut_syst_err,singleR_A_c12_mf_syst_err,singleR_A_c12_mf_tot_err,singleR_A_c12_src,singleR_A_c12_src_stat_err,singleR_A_c12_src_norm_syst_err,singleR_A_c12_src_RC_syst_err,singleR_A_c12_src_cut_syst_err,singleR_A_c12_src_syst_err,singleR_A_c12_src_tot_err,singleR_per_proton,singleR_per_proton_stat_err,singleR_per_proton_norm_syst_err,singleR_per_proton_RC_syst_err,singleR_per_proton_syst_err,singleR_per_proton_tot_err,doubleR,doubleR_stat_err,doubleR_norm_syst_err,doubleR_RC_syst_err,doubleR_cut_syst_err,doubleR_syst_err,doubleR_tot_err,doubleR_Jmodel,singleR_A_c12_mf_av18,singleR_A_c12_src_av18,doubleR_av18,singleR_A_c12_mf_osu,singleR_A_c12_src_osu,doubleR_osu,N,Z,A,NoZ,NmZoA\n') 

    # set output file to write CaFe triple [Ca40 Ca48 Fe54] / Ca48 single SRC ratios
    ofile2 = open('cafe_triplet.csv', 'w+')
    ofile2.write('# CaFe Triplet Numerical Ratios (%s) \n'%(npass))
    ofile2.write('# \n'
                '# Header Definitions: \n'
                '# target       : target A used in single ratio \n'
                '# singleR_A_ca48_mf   : single ratio of target A_MF / Ca48_MF (per proton) \n'
                '# singleR_A_ca48_src  : single ratio of target A_SRC / Ca48_SRC (per proton) \n'
                '# doubleR      : double ratio of target A(SRC/MF) relative to Ca48 (SRC/MF) \n'
                '# _stat_err    : absolute statistical uncertainty \n'
                '# _norm_syst_err : absolute systematic error due to normalization correction factors (a) hms/shms track efficiencies (b) live time (c) proton transmission (d) ca 40/48 cntm)\n'
                '# _RC_syst_err   : absolute systematic error due to radiative corrections \n'
                '# _cut_syst_err  : absolute systematic error due to cut sensitivity \n'
                '# _syst_err      : absolute systematic errors added in quadrature\n'
                '# _tot_err       : absolute total errors (statistical + systematics added in quadrature) \n'
                '# N: Z: A      : # of neutrons (N): protons(Z): nucleons (A): for target A \n'
                '# NoZ          : N/Z \n'
                '# NmZoA        : (N-Z)/A \n'                
                )
    ofile2.write('target,singleR_A_ca48_mf,singleR_A_ca48_mf_stat_err,singleR_A_ca48_mf_norm_syst_err,singleR_A_ca48_mf_RC_syst_err,singleR_A_ca48_mf_cut_syst_err,singleR_A_ca48_mf_syst_err,singleR_A_ca48_mf_tot_err,singleR_A_ca48_src,singleR_A_ca48_src_stat_err,singleR_A_ca48_src_norm_syst_err,singleR_A_ca48_src_RC_syst_err,singleR_A_ca48_src_cut_syst_err,singleR_A_ca48_src_syst_err,singleR_A_ca48_src_tot_err,doubleR,doubleR_stat_err,doubleR_norm_syst_err,doubleR_RC_syst_err,doubleR_cut_syst_err,doubleR_syst_err,doubleR_tot_err,N,Z,A,NoZ,NmZoA\n') 

    
    #--------------
    
    # read param
    T = np.array(df[(df['kin']=='SRC')]['T'])
    T_arr = [T] * 8
    N = np.array(df[(df['kin']=='SRC')]['N'])
    Z = np.array(df[(df['kin']=='SRC')]['Z'])
    A = np.array(df[(df['kin']=='SRC')]['A'])
    NoZ = N/Z
    NmZoA = (N-Z) / A

    # define for cafe triplet (Ca40, Ca48, Fe54)
    Nt = np.array([20.0, 28.0, 28.0])
    Zt = np.array([20.0, 20.0, 26.0])
    At = Nt + Zt
    NoZ_t = Nt/Zt
    NmZoA_t = (Nt - Zt ) / At
    
    #--------------------------------------------------------------------------------
    # ----- Read relative stats. and normalization systematic uncertainties
    #--------------------------------------------------------------------------------
    stat_rel_err_src = df[(df['kin']=='SRC')]['stat_rel_err']
    stat_rel_err_mf  = df[(df['kin']=='MF')]['stat_rel_err']

    norm_syst_rel_err_src = df[(df['kin']=='SRC')]['norm_syst_rel_err']
    norm_syst_rel_err_mf  = df[(df['kin']=='MF')]['norm_syst_rel_err']

    stat_rel_err_src_c12 = df[(df['kin']=='SRC') & (df['target']=='C12')]['stat_rel_err']
    stat_rel_err_mf_c12  = df[(df['kin']=='MF') & (df['target']=='C12')]['stat_rel_err']

    norm_syst_rel_err_src_c12 = df[(df['kin']=='SRC') & (df['target']=='C12')]['norm_syst_rel_err']
    norm_syst_rel_err_mf_c12  = df[(df['kin']=='MF') & (df['target']=='C12')]['norm_syst_rel_err']

    # read relative errors for cafe triplet (Ca40,48, Fe54) calculations
    stat_rel_err_mf_ca40 = float(df[(df['kin']=='MF') & (df['target']=='Ca40')]['stat_rel_err'])
    stat_rel_err_mf_ca48 = float(df[(df['kin']=='MF') & (df['target']=='Ca48')]['stat_rel_err'])
    stat_rel_err_mf_fe54 = float(df[(df['kin']=='MF') & (df['target']=='Fe54')]['stat_rel_err'])

    stat_rel_err_src_ca40 = float(df[(df['kin']=='SRC') & (df['target']=='Ca40')]['stat_rel_err'])
    stat_rel_err_src_ca48 = float(df[(df['kin']=='SRC') & (df['target']=='Ca48')]['stat_rel_err'])
    stat_rel_err_src_fe54 = float(df[(df['kin']=='SRC') & (df['target']=='Fe54')]['stat_rel_err'])

    norm_syst_rel_err_mf_ca40 = float(df[(df['kin']=='MF') & (df['target']=='Ca40')]['norm_syst_rel_err'])
    norm_syst_rel_err_mf_ca48 = float(df[(df['kin']=='MF') & (df['target']=='Ca48')]['norm_syst_rel_err'])
    norm_syst_rel_err_mf_fe54 = float(df[(df['kin']=='MF') & (df['target']=='Fe54')]['norm_syst_rel_err'])

    norm_syst_rel_err_src_ca40 = float(df[(df['kin']=='SRC') & (df['target']=='Ca40')]['norm_syst_rel_err'])
    norm_syst_rel_err_src_ca48 = float(df[(df['kin']=='SRC') & (df['target']=='Ca48')]['norm_syst_rel_err'])
    norm_syst_rel_err_src_fe54 = float(df[(df['kin']=='SRC') & (df['target']=='Fe54')]['norm_syst_rel_err'])


    #--------------------------------------------------------------------------------
    # ----- Read charge-normalized corrected yields from the final summary file -----
    #--------------------------------------------------------------------------------

    # NOTE: src_sigma_raw_proton -->  yield / (charge [mC] * transparency * target_thickness [g/cm2] * Z/A)

    # read raw cross sections (raw sigma or yield)
    src_sigma_raw_per_proton          = df[(df['kin']=='SRC')]['sigma_raw_per_proton']
    mf_sigma_raw_per_proton           = df[(df['kin']=='MF')]['sigma_raw_per_proton']
    
    src_sigma_raw_per_proton_C12      = df[(df['kin']=='SRC') & (df['target']=='C12')]['sigma_raw_per_proton']
    mf_sigma_raw_per_proton_C12       = df[(df['kin']=='MF')  & (df['target']=='C12')]['sigma_raw_per_proton']

    # read raw cross sections (raw sigma or yield) for cafe triplet Ca40,48, Fe54
    mf_sigma_raw_per_proton_Ca40      = df[(df['kin']=='MF') & (df['target']=='Ca40')]['sigma_raw_per_proton']
    mf_sigma_raw_per_proton_Ca48      = df[(df['kin']=='MF') & (df['target']=='Ca48')]['sigma_raw_per_proton']
    mf_sigma_raw_per_proton_Fe54      = df[(df['kin']=='MF') & (df['target']=='Fe54')]['sigma_raw_per_proton']

    src_sigma_raw_per_proton_Ca40      = df[(df['kin']=='SRC') & (df['target']=='Ca40')]['sigma_raw_per_proton']
    src_sigma_raw_per_proton_Ca48      = df[(df['kin']=='SRC') & (df['target']=='Ca48')]['sigma_raw_per_proton']
    src_sigma_raw_per_proton_Fe54      = df[(df['kin']=='SRC') & (df['target']=='Fe54')]['sigma_raw_per_proton']


    # read radiative correction factors
    rad_corr_src                      = np.array(df[(df['kin']=='SRC')]['rad_corr'])  # src_rad/src_norad
    rad_corr_mf                       = np.array(df[(df['kin']=='MF')]['rad_corr'])   # mf_rad/mf_norad


    # read radiative correction factors for specified targets to make it easier when taking ratios

    # for c12
    rad_corr_src_c12                  = float(df[(df['kin']=='SRC') & (df['target']=='C12')]['rad_corr'])
    rad_corr_mf_c12                   = float(df[(df['kin']=='MF')  & (df['target']=='C12')]['rad_corr']) 

    # (for cafe triplet)
    rad_corr_mf_ca40                  = float(df[(df['kin']=='MF') & (df['target']=='Ca40')]['rad_corr'])
    rad_corr_mf_ca48                  = float(df[(df['kin']=='MF') & (df['target']=='Ca48')]['rad_corr'])
    rad_corr_mf_fe54                  = float(df[(df['kin']=='MF') & (df['target']=='Fe54')]['rad_corr'])
    
    rad_corr_src_ca40                  = float(df[(df['kin']=='SRC') & (df['target']=='Ca40')]['rad_corr'])
    rad_corr_src_ca48                  = float(df[(df['kin']=='SRC') & (df['target']=='Ca48')]['rad_corr'])
    rad_corr_src_fe54                  = float(df[(df['kin']=='SRC') & (df['target']=='Fe54')]['rad_corr'])
  
    # calculate the absolute stats. and norm. syst error on the raw cross sections
    src_sigma_raw_per_proton_stat_err      = src_sigma_raw_per_proton * stat_rel_err_src
    src_sigma_raw_per_proton_norm_syst_err = src_sigma_raw_per_proton * norm_syst_rel_err_src 
    mf_sigma_raw_per_proton_stat_err       = mf_sigma_raw_per_proton  * stat_rel_err_mf
    mf_sigma_raw_per_proton_norm_syst_err  = mf_sigma_raw_per_proton  * norm_syst_rel_err_mf 

    src_sigma_raw_per_proton_C12_stat_err      = src_sigma_raw_per_proton_C12 * stat_rel_err_src_c12
    src_sigma_raw_per_proton_C12_norm_syst_err = src_sigma_raw_per_proton_C12 * norm_syst_rel_err_src_c12 
    mf_sigma_raw_per_proton_C12_stat_err       = mf_sigma_raw_per_proton_C12  * stat_rel_err_mf_c12
    mf_sigma_raw_per_proton_C12_norm_syst_err  = mf_sigma_raw_per_proton_C12  * norm_syst_rel_err_mf_c12 

    
    
    #-----------------------------------
    # ---- CALCULATE SINGLE RATIOS  ----
    #-----------------------------------

    #--- A_SRC / A_MF (per proton) ---


    # standard error propagation on ratio
    # for R = A/B,  dR = |R| sqrt( (dA/A)**2 + (dB/B)**2 ), assuming A,B are non-correlated

        
    # rad correction factor for SRC/MF single ratio on same nucleus
    # default rad corr relative error (dR/R) on single ratio to 5%,
    # since rad corr between SRC/MF of same nucleus A is very similar
    rad_corr_ratio_rel_err = 0.05

    rad_corr_ratio = (rad_corr_mf/rad_corr_src)

    # re-define rad corr ratios for C12 and Ca48 (to simplify double ratio calculations later on)
    rad_corr_ratio_c12 = (rad_corr_mf_c12/rad_corr_src_c12)
    rad_corr_ratio_ca48 = (rad_corr_mf_ca48/rad_corr_src_ca48)

    # calculate abssolute errors on single SRC / MF ratios
    singleR_per_proton               = (src_sigma_raw_per_proton.values/mf_sigma_raw_per_proton.values) * rad_corr_ratio               # apply radiative correction
    singleR_per_proton_RC_syst_err   = singleR_per_proton * rad_corr_ratio_rel_err                                                     # absolute syst. error due to rad corr ratio
    singleR_per_proton_stat_err      = singleR_per_proton * np.sqrt(stat_rel_err_src.values**2 +  stat_rel_err_mf.values**2)           # absolute stat. error 
    singleR_per_proton_norm_syst_err = singleR_per_proton * np.sqrt(norm_syst_rel_err_src.values**2 +  norm_syst_rel_err_mf.values**2) # absolute norm syst. error 

    # re-define the single ratio for C12 (just to make it easier later on when dividing by C12)
    singleR_per_proton_C12               = (src_sigma_raw_per_proton_C12.values/mf_sigma_raw_per_proton_C12.values) * rad_corr_ratio_c12               # apply radiative correction
    singleR_per_proton_C12_RC_syst_err   = singleR_per_proton_C12 * rad_corr_ratio_rel_err                                                             # absolute syst. error due to rad corr ratio
    singleR_per_proton_C12_stat_err      = singleR_per_proton_C12 * np.sqrt(stat_rel_err_src_c12.values**2 +  stat_rel_err_mf_c12.values**2)           # absolute stat. error 
    singleR_per_proton_C12_norm_syst_err = singleR_per_proton_C12 * np.sqrt(norm_syst_rel_err_src_c12.values**2 +  norm_syst_rel_err_mf_c12.values**2) # absolute norm syst. error 

    # re-define the single ratio for Ca48 (just to make it easier later on when dividing by Ca48 for the triplet)
    singleR_per_proton_Ca48               = (src_sigma_raw_per_proton_Ca48.values/mf_sigma_raw_per_proton_Ca48.values) * rad_corr_ratio_ca48  # apply radiative correction
    singleR_per_proton_Ca48_RC_syst_err   = singleR_per_proton_Ca48 * rad_corr_ratio_rel_err                                                  # absolute syst. error due to rad corr ratio
    singleR_per_proton_Ca48_stat_err      = singleR_per_proton_Ca48 * np.sqrt(stat_rel_err_src_ca48**2 +  stat_rel_err_mf_ca48**2)            # absolute stat. error 
    singleR_per_proton_Ca48_norm_syst_err = singleR_per_proton_Ca48 * np.sqrt(norm_syst_rel_err_src_ca48**2 +  norm_syst_rel_err_mf_ca48**2)  # absolute norm syst. error 

    
    #----- add errors in quadrature -----
    singleR_per_proton_syst_err = np.sqrt(singleR_per_proton_norm_syst_err**2 + singleR_per_proton_RC_syst_err**2)
    singleR_per_proton_tot_err  = np.sqrt(singleR_per_proton_stat_err**2 + singleR_per_proton_syst_err**2)
     
    #--- A_MF / C12_MF (per proton) -should be flat ---

    # cut sensitivity relative errors       be9/c12 b10/c12 b11/c12 c12/c12 ca40/c12 ca48/c12 fe54/c12 au197/c12 
    singleR_mf_cut_syst_rel_err = np.array([0.004,  0.002,  0.001,  0.0,    0.003,   0.005,   0.003,   0.007])  # fractional

    # rad correction factor for MF/C12_MF single ratio
    # has a ratio = 1 for light nuclei --> assign a 5% relative error on the ratio
    # othersiwe, if ratio != 1, such as heavy_A_MF / C12_MF --> assign 20% or correction factor
    # e.g., if ratio = 0.6 / 0.38 = 1.57  --> abs(1-1.57) =  0.57 (57% correction factor)
    # --> assign 0.57 * 0.20 = 0.114 (11%) as relative syst. error
    
    rad_corr_ratio_mf = ( rad_corr_mf_c12 / rad_corr_mf )

    # initialize zero arrays to be filled in depending on relative error value
    singleR_A_c12_mf               = np.zeros([len(rad_corr_ratio_mf)])
    singleR_A_c12_mf_RC_syst_err   = np.zeros([len(rad_corr_ratio_mf)])
    singleR_A_c12_mf_stat_err      = np.zeros([len(rad_corr_ratio_mf)])
    singleR_A_c12_mf_norm_syst_err = np.zeros([len(rad_corr_ratio_mf)])
    singleR_A_c12_mf_cut_syst_err  = np.zeros([len(rad_corr_ratio_mf)])

    
    # loop over rad corr ratio
    for i in range(len(rad_corr_ratio_mf)):
        if rad_corr_ratio_mf[i] == 1:
            rad_corr_ratio_mf_rel_err = 0.05  # default, if corr factor is same (as it is for light nuclei)

        else:
            rad_corr_ratio_mf_rel_err =  np.abs(1.-rad_corr_ratio_mf[i]) * 0.20
            #print('np.abs(1.-rad_corr_ratio_mf[i]) = ',np.abs(1.-rad_corr_ratio_mf[i]))
            #print('rad_corr_ratio_mf_rel_err = ',rad_corr_ratio_mf_rel_err )
        singleR_A_c12_mf[i]               =  ( mf_sigma_raw_per_proton.values[i] / mf_sigma_raw_per_proton_C12 ) * rad_corr_ratio_mf[i]
        singleR_A_c12_mf_RC_syst_err[i]       = singleR_A_c12_mf[i] * rad_corr_ratio_mf_rel_err  # absolute syst. error due to rad corr ratio
        singleR_A_c12_mf_stat_err[i]      = singleR_A_c12_mf[i] * np.sqrt(stat_rel_err_mf.values[i]**2 +  stat_rel_err_mf_c12**2) 
        singleR_A_c12_mf_norm_syst_err[i] = singleR_A_c12_mf[i] * np.sqrt(norm_syst_rel_err_mf.values[i]**2 +  norm_syst_rel_err_mf_c12**2)
        singleR_A_c12_mf_cut_syst_err[i]  = singleR_A_c12_mf[i] * singleR_mf_cut_syst_rel_err[i] 

    #------
    singleR_A_c12_mf_syst_err = np.sqrt(singleR_A_c12_mf_norm_syst_err**2 + singleR_A_c12_mf_RC_syst_err**2 + singleR_A_c12_mf_cut_syst_err**2)
    singleR_A_c12_mf_tot_err  = np.sqrt(singleR_A_c12_mf_stat_err**2      + singleR_A_c12_mf_syst_err**2)
    
    #--- A_SRC / C12_SRC (per proton) ---

    # cut sensitivity relative errors       be9/c12  b10/c12 b11/c12 c12/c12 ca40/c12 ca48/c12 fe54/c12 au197/c12
    singleR_src_cut_syst_rel_err = np.array([0.013,  0.016,  0.009,  0.0,    0.022,   0.031,   0.024,   0.049]) # fractional
    
    # similar treatment of radiative corrections as for the A_MF / C12_MF case above
    rad_corr_ratio_src = ( rad_corr_src_c12 / rad_corr_src )

    # initialize zero arrays to be filled in depending on relative error value
    singleR_A_c12_src               = np.zeros([len(rad_corr_ratio_src)])
    singleR_A_c12_src_RC_syst_err   = np.zeros([len(rad_corr_ratio_src)])
    singleR_A_c12_src_stat_err      = np.zeros([len(rad_corr_ratio_src)])
    singleR_A_c12_src_norm_syst_err = np.zeros([len(rad_corr_ratio_src)])
    singleR_A_c12_src_cut_syst_err  = np.zeros([len(rad_corr_ratio_src)])

    # loop over rad corr ratio
    for i in range(len(rad_corr_ratio_src)):
        if rad_corr_ratio_src[i] == 1:
            rad_corr_ratio_src_rel_err = 0.05  # default, if corr factor is same (as it is for light nuclei)

        else:
            rad_corr_ratio_src_rel_err =  np.abs(1.-rad_corr_ratio_src[i]) * 0.20
            #print(' np.abs(1.-rad_corr_ratio_src[i]) = ', np.abs(1.-rad_corr_ratio_src[i]))
            #print('rad_corr_ratio_src_rel_err = ',rad_corr_ratio_src_rel_err)
        singleR_A_c12_src[i]               =  ( src_sigma_raw_per_proton.values[i] / src_sigma_raw_per_proton_C12 ) * rad_corr_ratio_src[i]
        singleR_A_c12_src_RC_syst_err[i]       = singleR_A_c12_src[i] * rad_corr_ratio_src_rel_err  # absolute syst. error due to rad corr ratio
        singleR_A_c12_src_stat_err[i]      = singleR_A_c12_src[i] * np.sqrt(stat_rel_err_src.values[i]**2 +  stat_rel_err_src_c12**2)
        singleR_A_c12_src_norm_syst_err[i] = singleR_A_c12_src[i] * np.sqrt(norm_syst_rel_err_src.values[i]**2 +  norm_syst_rel_err_src_c12**2)
        singleR_A_c12_src_cut_syst_err[i]  = singleR_A_c12_src[i] * singleR_src_cut_syst_rel_err[i]

    #------
    singleR_A_c12_src_syst_err = np.sqrt(singleR_A_c12_src_norm_syst_err**2 + singleR_A_c12_src_RC_syst_err**2 + singleR_A_c12_src_cut_syst_err**2)
    singleR_A_c12_src_tot_err = np.sqrt(singleR_A_c12_src_stat_err**2 + singleR_A_c12_src_syst_err**2)


    #-------------------------------------------------
    # Calculate CaFe Triplet Ratio relative to Ca48
    #
    #--- A_MF / Ca48_MF (per proton) ---
    #-------------------------------------------------

    # cut sensitivity relative errors               ca40/ca48    ca48/ca48   fe54/ca48        (NOTE: need to be updates, as these are currently dummy placeholders)
    singleR_mf_triplet_cut_syst_rel_err = np.array([0.01,        0.0,        0.01])          # fractional

    # transparency systematic error on triple ratios (As per Larry's suggeation)    
    singleR_mf_triplet_T_syst_rel_err   = np.array([0.01,        0.0,        0.01])  # 1 % relative error on transparency for triplet

    
     # initialize zero arrays to be filled in depending on relative error value
    singleR_A_ca48_mf               = np.zeros([3])
    singleR_A_ca48_mf_RC_syst_err   = np.zeros([3])
    singleR_A_ca48_mf_stat_err      = np.zeros([3])
    singleR_A_ca48_mf_norm_syst_err = np.zeros([3])
    singleR_A_ca48_mf_cut_syst_err  = np.zeros([3])

    singleR_A_ca48_mf_syst_err      = np.zeros([3])
    singleR_A_ca48_mf_tot_err       = np.zeros([3])

    # apply radiative correction factor to single ratio A/Ca48
    '''
    print('mf_sigma_raw_per_proton_Ca40 =', mf_sigma_raw_per_proton_Ca40 )
    print('mf_sigma_raw_per_proton_Ca48 =', mf_sigma_raw_per_proton_Ca48 )
    print('rad_corr_mf_ca48 = ', rad_corr_mf_ca48)
    print('rad_corr_mf_ca40 = ', rad_corr_mf_ca40)
    '''
    singleR_A_ca48_mf[0]               =  ( float(mf_sigma_raw_per_proton_Ca40) / float(mf_sigma_raw_per_proton_Ca48) ) * (rad_corr_mf_ca48 / rad_corr_mf_ca40)
    singleR_A_ca48_mf[1]               =  ( float(mf_sigma_raw_per_proton_Ca48) / float(mf_sigma_raw_per_proton_Ca48) ) * (rad_corr_mf_ca48 / rad_corr_mf_ca48)
    singleR_A_ca48_mf[2]               =  ( float(mf_sigma_raw_per_proton_Fe54) / float(mf_sigma_raw_per_proton_Ca48) ) * (rad_corr_mf_ca48 / rad_corr_mf_fe54)

    # need to figure this out (what is the relative error on the radiative correction factor for (Ca40 Ca48 Fe54) / Ca48 single ratio?
    singleR_A_ca48_mf_RC_syst_err[0]       = 0  # singleR_A_ca48_mf[0] * rad_corr_ratio_mf_rel_err ???
    singleR_A_ca48_mf_RC_syst_err[1]       = 0
    singleR_A_ca48_mf_RC_syst_err[2]       = singleR_A_ca48_mf[2] * 0.025   # (Larry) suggests we add 2.5% systematic uncertainty to Fe for the rad corr. (Is it 2.5 % relative uncertainty on Fe/Ca48) ?


    # calculate absolute error on triplet ratio due to transparency
    singleR_A_ca48_mf_T_syst_err      = singleR_A_ca48_mf  *  singleR_mf_triplet_T_syst_rel_err 

    # calculate
    singleR_A_ca48_mf_stat_err[0]      = singleR_A_ca48_mf[0] * np.sqrt(stat_rel_err_mf_ca40**2 +  stat_rel_err_mf_ca48**2)
    singleR_A_ca48_mf_stat_err[1]      = singleR_A_ca48_mf[1] * np.sqrt(stat_rel_err_mf_ca48**2 +  stat_rel_err_mf_ca48**2)
    singleR_A_ca48_mf_stat_err[2]      = singleR_A_ca48_mf[2] * np.sqrt(stat_rel_err_mf_fe54**2 +  stat_rel_err_mf_ca48**2)

    singleR_A_ca48_mf_norm_syst_err[0]      = singleR_A_ca48_mf[0] * np.sqrt(norm_syst_rel_err_mf_ca40**2 +  norm_syst_rel_err_mf_ca48**2)
    singleR_A_ca48_mf_norm_syst_err[1]      = singleR_A_ca48_mf[1] * np.sqrt(norm_syst_rel_err_mf_ca48**2 +  norm_syst_rel_err_mf_ca48**2)
    singleR_A_ca48_mf_norm_syst_err[2]      = singleR_A_ca48_mf[2] * np.sqrt(norm_syst_rel_err_mf_fe54**2 +  norm_syst_rel_err_mf_ca48**2)

    singleR_A_ca48_mf_cut_syst_err[0]  = singleR_A_ca48_mf[0] * singleR_mf_triplet_cut_syst_rel_err[0]
    singleR_A_ca48_mf_cut_syst_err[1]  = singleR_A_ca48_mf[1] * singleR_mf_triplet_cut_syst_rel_err[1]
    singleR_A_ca48_mf_cut_syst_err[2]  = singleR_A_ca48_mf[2] * singleR_mf_triplet_cut_syst_rel_err[2]

    #------ add errors in quadrature ----
    singleR_A_ca48_mf_syst_err[0] = np.sqrt(singleR_A_ca48_mf_norm_syst_err[0]**2 + singleR_A_ca48_mf_RC_syst_err[0]**2 + singleR_A_ca48_mf_cut_syst_err[0]**2 + singleR_A_ca48_mf_T_syst_err[0]**2)
    singleR_A_ca48_mf_tot_err[0]  = np.sqrt(singleR_A_ca48_mf_stat_err[0]**2 + singleR_A_ca48_mf_syst_err[0]**2)

    singleR_A_ca48_mf_syst_err[1] = np.sqrt(singleR_A_ca48_mf_norm_syst_err[1]**2 + singleR_A_ca48_mf_RC_syst_err[1]**2 + singleR_A_ca48_mf_cut_syst_err[1]**2 + singleR_A_ca48_mf_T_syst_err[1]**2)
    singleR_A_ca48_mf_tot_err[1]  = np.sqrt(singleR_A_ca48_mf_stat_err[1]**2 + singleR_A_ca48_mf_syst_err[1]**2)

    singleR_A_ca48_mf_syst_err[2] = np.sqrt(singleR_A_ca48_mf_norm_syst_err[2]**2 + singleR_A_ca48_mf_RC_syst_err[2]**2 + singleR_A_ca48_mf_cut_syst_err[2]**2 + singleR_A_ca48_mf_T_syst_err[2]**2)
    singleR_A_ca48_mf_tot_err[2]  = np.sqrt(singleR_A_ca48_mf_stat_err[2]**2 + singleR_A_ca48_mf_syst_err[2]**2)

    #----------------------------------
    
    #-------------------------------------------------
    # Calculate CaFe Triplet Ratio relative to Ca48
    #
    #--- A_SRC / Ca48_SRC (per proton) ---
    #-------------------------------------------------

    # cut sensitivity relative errors                ca40/ca48  ca48/ca48 fe54/ca48        (NOTE: need to be updates, as these are currently dummy placeholders)
    singleR_src_triplet_cut_syst_rel_err = np.array([0.01,        0.0,    0.01])          # fractional

    # transparency systematic error on triple ratios (As per Larry's suggeation)    
    singleR_src_triplet_T_syst_rel_err   = np.array([0.01,        0.0,        0.01])  # 1 % relative error on transparency for triplet

     # initialize zero arrays to be filled in depending on relative error value
    singleR_A_ca48_src               = np.zeros([3])
    singleR_A_ca48_src_RC_syst_err   = np.zeros([3])
    singleR_A_ca48_src_stat_err      = np.zeros([3])
    singleR_A_ca48_src_norm_syst_err = np.zeros([3])
    singleR_A_ca48_src_cut_syst_err  = np.zeros([3])

    singleR_A_ca48_src_syst_err      = np.zeros([3])
    singleR_A_ca48_src_tot_err       = np.zeros([3])

    
    # apply radiative correction factor to single ratio A/Ca48
    '''
    print('src_sigma_raw_per_proton_Ca40 =', src_sigma_raw_per_proton_Ca40 )
    print('src_sigma_raw_per_proton_Ca48 =', src_sigma_raw_per_proton_Ca48 )
    print('rad_corr_src_ca48 = ', rad_corr_src_ca48)
    print('rad_corr_src_ca40 = ', rad_corr_src_ca40)
    '''
    
    singleR_A_ca48_src[0]               =  ( float(src_sigma_raw_per_proton_Ca40) / float(src_sigma_raw_per_proton_Ca48) ) * (rad_corr_src_ca48 / rad_corr_src_ca40)
    singleR_A_ca48_src[1]               =  ( float(src_sigma_raw_per_proton_Ca48) / float(src_sigma_raw_per_proton_Ca48) ) * (rad_corr_src_ca48 / rad_corr_src_ca48)
    singleR_A_ca48_src[2]               =  ( float(src_sigma_raw_per_proton_Fe54) / float(src_sigma_raw_per_proton_Ca48) ) * (rad_corr_src_ca48 / rad_corr_src_fe54)

    # need to figure this out (what is the relative error on the radiative correction factor for (Ca40 Ca48 Fe54) / Ca48 single ratio?
    singleR_A_ca48_src_RC_syst_err[0]       = 0  # singleR_A_ca48_src[0] * rad_corr_ratio_src_rel_err ???
    singleR_A_ca48_src_RC_syst_err[1]       = 0
    singleR_A_ca48_src_RC_syst_err[2]       = singleR_A_ca48_src[2] * 0.025   # (Larry) suggests we add 2.5% systematic uncertainty to Fe for the rad corr. (Is it 2.5 % relative uncertainty on Fe/Ca48) ?

    # calculate absolute error on triplet ratio due to transparency
    singleR_A_ca48_src_T_syst_err      = singleR_A_ca48_src  *  singleR_src_triplet_T_syst_rel_err 

    # calculate absolute errors on stats, norm_syst, cut_syst
    singleR_A_ca48_src_stat_err[0]      = singleR_A_ca48_src[0] * np.sqrt(stat_rel_err_src_ca40**2 +  stat_rel_err_src_ca48**2)
    singleR_A_ca48_src_stat_err[1]      = singleR_A_ca48_src[1] * np.sqrt(stat_rel_err_src_ca48**2 +  stat_rel_err_src_ca48**2)
    singleR_A_ca48_src_stat_err[2]      = singleR_A_ca48_src[2] * np.sqrt(stat_rel_err_src_fe54**2 +  stat_rel_err_src_ca48**2)

    singleR_A_ca48_src_norm_syst_err[0]      = singleR_A_ca48_src[0] * np.sqrt(norm_syst_rel_err_src_ca40**2 +  norm_syst_rel_err_src_ca48**2)
    singleR_A_ca48_src_norm_syst_err[1]      = singleR_A_ca48_src[1] * np.sqrt(norm_syst_rel_err_src_ca48**2 +  norm_syst_rel_err_src_ca48**2)
    singleR_A_ca48_src_norm_syst_err[2]      = singleR_A_ca48_src[2] * np.sqrt(norm_syst_rel_err_src_fe54**2 +  norm_syst_rel_err_src_ca48**2)

    singleR_A_ca48_src_cut_syst_err[0]  = singleR_A_ca48_src[0] * singleR_src_triplet_cut_syst_rel_err[0]
    singleR_A_ca48_src_cut_syst_err[1]  = singleR_A_ca48_src[1] * singleR_src_triplet_cut_syst_rel_err[1]
    singleR_A_ca48_src_cut_syst_err[2]  = singleR_A_ca48_src[2] * singleR_src_triplet_cut_syst_rel_err[2]

    #------ add errors in quadrature ----
    singleR_A_ca48_src_syst_err[0] = np.sqrt(singleR_A_ca48_src_norm_syst_err[0]**2 + singleR_A_ca48_src_RC_syst_err[0]**2 + singleR_A_ca48_src_cut_syst_err[0]**2 + singleR_A_ca48_src_T_syst_err[0]**2 )
    singleR_A_ca48_src_tot_err[0]  = np.sqrt(singleR_A_ca48_src_stat_err[0]**2 + singleR_A_ca48_src_syst_err[0]**2)

    singleR_A_ca48_src_syst_err[1] = np.sqrt(singleR_A_ca48_src_norm_syst_err[1]**2 + singleR_A_ca48_src_RC_syst_err[1]**2 + singleR_A_ca48_src_cut_syst_err[1]**2 + singleR_A_ca48_src_T_syst_err[1]**2 )
    singleR_A_ca48_src_tot_err[1]  = np.sqrt(singleR_A_ca48_src_stat_err[1]**2 + singleR_A_ca48_src_syst_err[1]**2)

    singleR_A_ca48_src_syst_err[2] = np.sqrt(singleR_A_ca48_src_norm_syst_err[2]**2 + singleR_A_ca48_src_RC_syst_err[2]**2 + singleR_A_ca48_src_cut_syst_err[2]**2 + singleR_A_ca48_src_T_syst_err[2]**2 )
    singleR_A_ca48_src_tot_err[2]  = np.sqrt(singleR_A_ca48_src_stat_err[2]**2 + singleR_A_ca48_src_syst_err[2]**2)

    #----------------------------------

    #---------------------------------------------------------
    # ---- CALCULATE DOUBLE RATIOS A_SRC/MF / Ca48_SRC/MF ----
    #---------------------------------------------------------

    #                           --------------------
    #        Be9  B10  B11  C12 | Ca40  Ca48  Fe54 | Au197
    # index: 0    1    2    3   | 4     5     6    | 7   
    #                           --------------------
    

    # cut sensitivity relative errors    ca40/ca48   ca48/ca48   fe54/ca48  (Noah needs to give updated numbers on this uncertainty)
    doubleRt_cut_syst_rel_err = np.array([0.01,       0.0,        0.01]) # fractional 

    # radiative correction factor for double ratio
    doubleRt_RC_corr_factor = rad_corr_ratio[4:7] / rad_corr_ratio_ca48

    # initialize zero arrays to be filled in depending on relative error value
    doubleRt               = np.zeros([3])
    doubleRt_stat_err      = np.zeros([3])
    doubleRt_RC_syst_err   = np.zeros([3])
    doubleRt_norm_syst_err = np.zeros([3])
    doubleRt_cut_syst_err  = np.zeros([3])

    doubleRt_RC_corr_factor_rel_err = np.zeros([3])
    
    # loop over rad corr factor and apply relative error on rad corr ( NEED to be UPDATED )
    for i in range(3):
        if doubleRt_RC_corr_factor[i] == 1:
            doubleRt_RC_corr_factor_rel_err[i] = 0.01  # default
            #print('doubleR_RC_corr_factor_rel_err = ',doubleRt_RC_corr_factor_rel_err[i])
        else:
            doubleRt_RC_corr_factor_rel_err[i] =  np.abs(1.-doubleRt_RC_corr_factor[i]) * 0.20
            #print('doubleR_RC_corr_factor_rel_err = ',doubleRt_RC_corr_factor_rel_err[i])
        
    doubleRt                = singleR_per_proton[4:7] / singleR_per_proton_Ca48
    doubleRt_stat_err       = doubleRt * np.sqrt( (singleR_per_proton_stat_err[4:7]/singleR_per_proton[4:7])**2 +  (singleR_per_proton_Ca48_stat_err/singleR_per_proton_Ca48)**2 )
    doubleRt_RC_syst_err    = doubleRt * doubleRt_RC_corr_factor_rel_err
    doubleRt_norm_syst_err  = doubleRt * np.sqrt( (singleR_per_proton_norm_syst_err[4:7]/singleR_per_proton[4:7])**2 + (singleR_per_proton_Ca48_norm_syst_err/singleR_per_proton_Ca48)**2  )
    doubleRt_cut_syst_err   = doubleRt * doubleRt_cut_syst_rel_err

    #print('doubleR_RC_corr_factor = ',doubleRt_RC_corr_factor)
    #print('doubleR_RC_corr_factor_rel_err = ',doubleRt_RC_corr_factor_rel_err)

    # -------- add errors in quadrature
    doubleRt_syst_err = np.sqrt(doubleRt_norm_syst_err**2 + doubleRt_RC_syst_err**2 + doubleRt_cut_syst_err**2)
    doubleRt_tot_err = np.sqrt(doubleRt_stat_err**2 + doubleRt_syst_err**2)
    #-----------------------------------
    
    #-----------------------------------------------
    #
    # Write CaFe Triplet Ratios to Numerical File
    #
    #-----------------------------------------------
    
    # read target varibale (isolate only for a kin setting, since double SRC/MF ratios being taken)
    targ = np.array(['Ca40', 'Ca48', 'Fe54'])
    
    # loop over each target (to write numerical values to file)
    for i in np.arange(len(targ)):

        if(targ[i]=="Ca48"):
            singleR_A_ca48_mf_stat_err[i]=singleR_A_ca48_mf_norm_syst_err[i]=singleR_A_ca48_mf_RC_syst_err[i]=singleR_A_ca48_mf_cut_syst_err[i]=singleR_A_ca48_mf_syst_err[i]=singleR_A_ca48_mf_tot_err[i]=0.0
            singleR_A_ca48_src_stat_err[i]=singleR_A_ca48_src_norm_syst_err[i]=singleR_A_ca48_src_RC_syst_err[i]=singleR_A_ca48_src_cut_syst_err[i]=singleR_A_ca48_src_syst_err[i]=singleR_A_ca48_src_tot_err[i]=0.0
            doubleRt_stat_err[i]=doubleRt_norm_syst_err[i]=doubleRt_RC_syst_err[i]=doubleRt_cut_syst_err[i]=doubleRt_syst_err[i]=doubleRt_tot_err[i]=0.0
        ofile2.write('%s,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.1f,%.1f,%.1f,%.3f,%.3f\n' % (targ[i],
                                                       singleR_A_ca48_mf[i],singleR_A_ca48_mf_stat_err[i],singleR_A_ca48_mf_norm_syst_err[i],singleR_A_ca48_mf_RC_syst_err[i],singleR_A_ca48_mf_cut_syst_err[i],singleR_A_ca48_mf_syst_err[i],singleR_A_ca48_mf_tot_err[i],
                                                       singleR_A_ca48_src[i],singleR_A_ca48_src_stat_err[i],singleR_A_ca48_src_norm_syst_err[i],singleR_A_ca48_src_RC_syst_err[i],singleR_A_ca48_src_cut_syst_err[i],singleR_A_ca48_src_syst_err[i],singleR_A_ca48_src_tot_err[i],
                                                       doubleRt[i],doubleRt_stat_err[i],doubleRt_norm_syst_err[i],doubleRt_RC_syst_err[i],doubleRt_cut_syst_err[i],doubleRt_syst_err[i],doubleRt_tot_err[i],                                                                        
                                                       Nt[i], Zt[i], At[i], NoZ_t[i], NmZoA_t[i]) ) 
        
       
    ofile2.close()

    #-------
    
    #--------------------------------------------------------
    # ---- CALCULATE DOUBLE RATIOS A_SRC/MF / C12_SRC/MF ----
    #--------------------------------------------------------

    #        Be9  B10  B11  C12  Ca40  Ca48  Fe54  Au197
    # index: 0    1    2    3    4     5     6     7   
    
    # cut sensitivity relative errors    be9/c12  b10/c12 b11/c12 c12/c12 ca40/c12 ca48/c12 fe54/c12 au197/c12
    doubleR_cut_syst_rel_err = np.array([0.014,  0.014,  0.009,  0.0,    0.022,   0.031,   0.024,   0.049]) # fractional

    # radiative correction factor for double ratio (this factor was already applied to the src/mf in the single ratios, but here is re-defined
    # for purposes or alculating rad corr factor relative error on double ratio)
    doubleR_RC_corr_factor = rad_corr_ratio / rad_corr_ratio_c12

    # initialize zero arrays to be filled in depending on relative error value
    doubleR               = np.zeros([len(doubleR_RC_corr_factor)])
    doubleR_stat_err      = np.zeros([len(doubleR_RC_corr_factor)])
    doubleR_RC_syst_err   = np.zeros([len(doubleR_RC_corr_factor)])
    doubleR_norm_syst_err = np.zeros([len(doubleR_RC_corr_factor)])
    doubleR_cut_syst_err  = np.zeros([len(doubleR_RC_corr_factor)])

    doubleR_RC_corr_factor_rel_err = np.zeros([len(doubleR_RC_corr_factor)])
    
    # loop over rad corr factor
    for i in range(len(doubleR_RC_corr_factor)):
        if doubleR_RC_corr_factor[i] == 1:
            doubleR_RC_corr_factor_rel_err[i] = 0.01  # default
            #print('doubleR_RC_corr_factor_rel_err = ',doubleR_RC_corr_factor_rel_err[i])
        else:
            doubleR_RC_corr_factor_rel_err[i] =  np.abs(1.-doubleR_RC_corr_factor[i]) * 0.20
            #print('doubleR_RC_corr_factor_rel_err = ',doubleR_RC_corr_factor_rel_err[i])
        
    doubleR                = singleR_per_proton / singleR_per_proton_C12
    doubleR_stat_err       = doubleR * np.sqrt( (singleR_per_proton_stat_err/singleR_per_proton)**2 +  (singleR_per_proton_C12_stat_err/singleR_per_proton_C12)**2 )
    doubleR_RC_syst_err    = doubleR * doubleR_RC_corr_factor_rel_err
    doubleR_norm_syst_err  = doubleR * np.sqrt( (singleR_per_proton_norm_syst_err/singleR_per_proton)**2 + (singleR_per_proton_C12_norm_syst_err/singleR_per_proton_C12)**2  )
    doubleR_cut_syst_err   = doubleR * doubleR_cut_syst_rel_err

    #print('doubleR_RC_corr_factor = ',doubleR_RC_corr_factor)
    #print('doubleR_RC_corr_factor_rel_err = ',doubleR_RC_corr_factor_rel_err)

    # -------- add errors in quadrature
    doubleR_syst_err = np.sqrt(doubleR_norm_syst_err**2 + doubleR_RC_syst_err**2 + doubleR_cut_syst_err**2)
    doubleR_tot_err = np.sqrt(doubleR_stat_err**2 + doubleR_syst_err**2)
    #-----------------------------------
    

    #----------------------------------
    # Justin's model (numerical data)
    #----------------------------------
    #                                   be9/c12  b10/c12  b11/c12  c12/c12     ca40/c12  ca48/c12  fe54/c12  au197/c12
    doubleR_Jmodel           = np.array([0.823,   0.896,   1.056,   np.nan,     1.197,    1.545,    1.374,    1.742])

    # AV18 model
    singleR_A_c12_mf_av18    = np.array([1.029,   1.022,   1.005,   np.nan,     1.021,    np.nan,   np.nan,   np.nan])
    singleR_A_c12_src_av18   = np.array([0.934,   0.896,   1.011,   np.nan,     0.854,    np.nan,   np.nan,   np.nan])
    doubleR_av18             = np.array([0.908,   0.877,   1.006,   np.nan,     0.837,    np.nan,   np.nan,   np.nan])

    # OSU model
    singleR_A_c12_mf_osu     = np.array([np.nan,  np.nan,  np.nan,  np.nan,     0.983,    0.962,    0.974,    np.nan])
    singleR_A_c12_src_osu    = np.array([np.nan,  np.nan,  np.nan,  np.nan,     1.158,    1.347,    1.178,    np.nan])
    doubleR_osu              = np.array([np.nan,  np.nan,  np.nan,  np.nan,     1.179,    1.401,    1.21,     np.nan])


    
    # read target varibale (isolate only for a kin setting, since double SRC/MF ratios being taken)
    targ = np.array(df[(df['kin']=='SRC')]['target'])
    
    # loop over each target (to write numerical values to file)
    for i in np.arange(len(targ)):

        if(targ[i]=="C12"):
            singleR_A_c12_mf_stat_err[i]=singleR_A_c12_mf_norm_syst_err[i]=singleR_A_c12_mf_RC_syst_err[i]=singleR_A_c12_mf_cut_syst_err[i]=singleR_A_c12_mf_syst_err[i]=singleR_A_c12_mf_tot_err[i]=0.0
            singleR_A_c12_src_stat_err[i]=singleR_A_c12_src_norm_syst_err[i]=singleR_A_c12_src_RC_syst_err[i]=singleR_A_c12_src_cut_syst_err[i]=singleR_A_c12_src_syst_err[i]=singleR_A_c12_src_tot_err[i]=0.0
            singleR_per_proton_stat_err[i]=singleR_per_proton_norm_syst_err[i]=singleR_per_proton_RC_syst_err[i]=singleR_per_proton_syst_err[i]=singleR_per_proton_tot_err[i]=0.0
            doubleR_stat_err[i]=doubleR_norm_syst_err[i]=doubleR_RC_syst_err[i]=doubleR_cut_syst_err[i]=doubleR_syst_err[i]=doubleR_tot_err[i]=0.0
            
        ofile.write('%s,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.1f,%.1f,%.1f,%.3f,%.3f\n' % (targ[i],
                                                       singleR_A_c12_mf[i],singleR_A_c12_mf_stat_err[i],singleR_A_c12_mf_norm_syst_err[i],singleR_A_c12_mf_RC_syst_err[i],singleR_A_c12_mf_cut_syst_err[i],singleR_A_c12_mf_syst_err[i],singleR_A_c12_mf_tot_err[i],
                                                       singleR_A_c12_src[i],singleR_A_c12_src_stat_err[i],singleR_A_c12_src_norm_syst_err[i],singleR_A_c12_src_RC_syst_err[i],singleR_A_c12_src_cut_syst_err[i],singleR_A_c12_src_syst_err[i],singleR_A_c12_src_tot_err[i],
                                                       singleR_per_proton[i],singleR_per_proton_stat_err[i],singleR_per_proton_norm_syst_err[i],singleR_per_proton_RC_syst_err[i],singleR_per_proton_syst_err[i],singleR_per_proton_tot_err[i],
                                                       doubleR[i],doubleR_stat_err[i],doubleR_norm_syst_err[i],doubleR_RC_syst_err[i],doubleR_cut_syst_err[i],doubleR_syst_err[i],doubleR_tot_err[i],
                                                       doubleR_Jmodel[i],singleR_A_c12_mf_av18[i],singleR_A_c12_src_av18[i],doubleR_av18[i],singleR_A_c12_mf_osu[i],singleR_A_c12_src_osu[i],doubleR_osu[i],
                                                       N[i], Z[i], A[i], NoZ[i], NmZoA[i]) ) 
        
       
    ofile.close()

    
    

    
make_final_summary()

