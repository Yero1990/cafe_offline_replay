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
    ofname = 'cafe_final_summary_%s.csv' %(npass) 
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
            alpha=-0.24
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
            mult_trk_eff_err = df['multi_track_eff_err'] 

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
            # (dY_corr / Y_corr)_syst = (df1/f1)**2 + (df2/f2)**2 +  . . . (dfn/fn)**2 --> relative norm. syst is sum in quadrature of relative error in corr. factors

            # calculate relative systematic uncertainty on hms_trk_eff
            f1_per_run = 1./hms_trk_eff
            df1_per_run = (1/hms_trk_eff**2) * hms_trk_eff_err
            df1_f1_per_run = df1_per_run / f1_per_run
            df1_f1_sq = np.sum(np.square(df1_f1_per_run))  # sum of the squared of the relative err per run (add in quadrature)

            # calculate relative systematic uncertainty on shms_trk_eff
            f2_per_run = 1./shms_trk_eff
            df2_per_run = (1/shms_trk_eff**2) * shms_trk_eff_err
            df2_f2_per_run = df2_per_run / f2_per_run
            df2_f2_sq = np.sum(np.square(df2_f2_per_run))  

            # calculate relative systematic uncertainty on total_LT
            f3_per_run = 1./total_LT
            df3_per_run = (1/total_LT**2) * total_LT_err
            df3_f3_per_run = df3_per_run / f3_per_run
            df3_f3_sq = np.sum(np.square(df3_f3_per_run))  

            # calculate relative systematic uncertainty on multi-track eff 
            f4_per_run = 1./mult_trk_eff
            df4_per_run = (1/mult_trk_eff**2) * mult_trk_eff_err
            df4_f4_per_run = df4_per_run / f4_per_run
            df4_f4_sq = np.sum(np.square(df4_f4_per_run))  

            # calculate relative systematic uncertainty on proton transmission factor (overall correction, not run-by-run)
            f5 = 1./pTrms_eff  # obtained from Noah Swan's studies of proton absorption
            df5 = pTrms_eff_err
            df5_f5_sq = (df5 / f5)**2
                
            # calculate overall relative norm. systematic unc on the relative yield
            norm_syst_rel_err = np.sqrt( df1_f1_sq + df2_f2_sq + df3_f3_sq + df4_f4_sq + df5_f5_sq)

            # --- apply Ca40 contamination systematic correction (the correction is pre-determined for either MF or SRC during get_cntmEff() function call
            if(target[idx]=='Ca40'):
                f6 = cntm_eff
                df6 = cntm_eff_err
                df6_f6_per_run = df6_per_run / f6_per_run
                df6_f6_sq = np.sum(np.square(df6_f6_per_run))  

                norm_syst_rel_err = np.sqrt( df1_f1_sq + df2_f2_sq + df3_f3_sq + df4_f4_sq + df5_f5_sq + df6_f6_sq)

            # --- apply Ca48 contamination systematic correction (the correction is pre-determined for either MF or SRC during get_cntmEff() function call
            if(target[idx]=='Ca48'):
                f7 = cntm_eff
                df7 = cntm_eff_err
                df7_f7_per_run = df7_per_run / f7_per_run
                df7_f7_sq = np.sum(np.square(df7_f7_per_run))  

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
            real_Yield_eff       = real_Yield * f1_per_run * f2_per_run * f3_per_run * f4_per_run
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


            show_plots=True
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
                
            #----------------------------------------
            # WRITE CAFE NUMERICAL DATA
            #----------------------------------------
     
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
    #         double_ratio = ( A_SRC / A_MF ) / ( C12_SRC / C12_MF )
    
    
    # read input final summary file
    df = pd.read_csv(ifname, comment='#')


    # set output file to write double ratio numerical values
    ofile = open(ofname, 'w+')
    ofile.write('# CaFe numerical ratios (%s) \n'%(npass))
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# target       : target A used in single or double ratio \n'
                '# singleR_A_c12_mf      : single ratio of target A_MF / C12_MF (per proton)  \n'
                '# singleR_A_c12_mf_err  : uncertainty in single ratio \n'
                '# singleR_A_c12_src      : single ratio of target A_SRC / C12_SRC (per proton) \n'
                '# singleR_A_c12_src_err  : uncertainty in single ratio \n'
                '# singleR_per_proton      : single ratio of target A(SRC/MF) (per proton) \n'
                '# singleR_per_proton_err  : uncertainty in single ratio (per proton) \n'
                '# doubleR      : double ratio of target A(SRC/MF) relative to C12 (SRC/MF) \n'
                '# doubleR_err  : uncertainty in double ratio \n'
                '# N: Z: A      : # of neutrons (N): protons(Z): nucleons (A): for target A \n'
                '# NoZ          : N/Z \n'
                '# NmZoA        : (N-Z)/A \n'                
                )
    ofile.write('target,singleR_A_c12_mf,singleR_A_c12_mf_err,singleR_A_c12_src,singleR_A_c12_src_err,singleR_per_proton,singleR_per_proton_err,doubleR,doubleR_err,N,Z,A,NoZ,NmZoA\n') 

    # read param
    T = np.array(df[(df['kin']=='SRC')]['T'])
    T_arr = [T] * 8
    N = np.array(df[(df['kin']=='SRC')]['N'])
    Z = np.array(df[(df['kin']=='SRC')]['Z'])
    A = np.array(df[(df['kin']=='SRC')]['A'])
    NoZ = N/Z
    NmZoA = (N-Z) / A
    
    #--------------------------------------------------------------------------------
    # ----- Read charge-normalized corrected yields from the final summary file -----
    #--------------------------------------------------------------------------------

    # NOTE: src_sigma_raw_proton -->  yield / (charge [mC] * transparency * target_thickness [g/cm2] * Z/A)

    src_sigma_raw_per_proton          = df[(df['kin']=='SRC')]['sigma_raw_per_proton']
    mf_sigma_raw_per_proton           = df[(df['kin']=='MF')]['sigma_raw_per_proton']
    
    src_sigma_raw_per_proton_C12      = df[(df['kin']=='SRC') & (df['target']=='C12')]['sigma_raw_per_proton']
    mf_sigma_raw_per_proton_C12       = df[(df['kin']=='MF')  & (df['target']=='C12')]['sigma_raw_per_proton']

    rad_corr_src                      = df[(df['kin']=='SRC')]['rad_corr']  # src_rad/src_norad
    rad_corr_mf                       = df[(df['kin']=='MF')]['rad_corr']   # mf_rad/mf_norad

    rad_corr_src_c12                  = df[(df['kin']=='SRC') & (df['target']=='C12')]['rad_corr'] 
    rad_corr_mf_c12                   = df[(df['kin']=='MF')  & (df['target']=='C12')]['rad_corr'] 


    #-----------------------------------
    # ---- CALCULATE SINGLE RATIOS  ----
    #-----------------------------------

    #--- A_SRC / A_MF (per proton) ---

    # rad correction factor for SRC/MF single ratio on same nucleus
    # default rad corr relative error (dR/R) on single ratio to 5%,
    # since rad corr between SRC/MF of same nucleus A is very similar

    rad_corr_ratio = (rad_corr_mf/rad_corr_src)
    rad_corr_ratio_rel_err = 0.05  
    singleR_per_proton     = (src_sigma_raw_per_proton/mf_sigma_raw_per_proton) * rad_corr_ratio
    singleR_per_proton_RC_syst = singleR_per_proton * rad_corr_ratio_rel_err  # absolute syst. error due to rad corr ratio

    #--> NEED TO FINISH singleR_per_proton_C12     = (src_sigma_raw_per_proton_C12/mf_sigma_raw_per_proton_C12) * rad_corr_ratio
    # --> NEED TO FINISH singleR_per_proton_RC_syst_C12 = singleR_per_proton * rad_corr_ratio_rel_err  # absolute syst. error due to rad corr ratio
    

    #--- A_MF / C12_MF (per proton) -should be flat ---
    # rad correction factor for MF/C12_MF single ratio
    # has a ratio = 1 for light nuclei --> assign a 5% relative error on the ratio
    # othersiwe, if ratio != 1, such as heavy_A_MF / C12_MF --> assign 20% or correction factor
    # e.g., if ratio = 0.6 / 0.38 = 1.57  --> abs(1-1.57) =  0.57 (57% correction factor)
    # --> assign 0.57 * 0.20 = 0.114 (11%) as relative syst. error
    
    rad_corr_ratio_mf = ( rad_corr_mf_c12 / rad_corr_mf ) 
    rad_corr_ratio_mf_rel_err = 0.05  # default

    if(rad_corr_ratio_mf!=1):
        rad_corr_ratio_mf_rel_err =  np.abs(1.-rad_corr_ratio_mf) * 0.20
    
    singleR_A_c12_mf  =  ( mf_sigma_raw_per_proton / mf_sigma_raw_per_proton_C12 ) * rad_corr_ratio_mf
    singleR_A_c12_mf_RC_syst = singleR_A_c12_mf * rad_corr_ratio_mf_rel_err  # absolute syst. error due to rad corr ratio

    
    #--- A_SRC / C12_SRC (per proton) ---
    # similar treatment of radiative corrections as for the A_MF / C12_MF case above

    rad_corr_ratio_src = ( rad_corr_src_c12 / rad_corr_src ) 
    rad_corr_ratio_src_rel_err = 0.05  # default

    if(rad_corr_ratio_src!=1):
        rad_corr_ratio_src_rel_err =  np.abs(1.-rad_corr_ratio_src) * 0.20
    
    singleR_A_c12_src     =  ( src_sigma_raw_per_proton / src_sigma_raw_per_proton_C12 ) * rad_corr_ratio_src
    singleR_A_c12_src_RC_syst = singleR_A_c12_src * rad_corr_ratio_src_rel_err  # absolute syst. error due to rad corr ratio


    
    #--------------------------------------------------------
    # ---- CALCULATE DOUBLE RATIOS A_SRC/MF / C12_SRC/MF ----
    #--------------------------------------------------------

    doubleR     = singleR_per_proton / singleR_per_proton_C12
    doubleR_val = unumpy.nominal_values(doubleR)
    doubleR_err = unumpy.std_devs(doubleR)
    
    '''
    doubleR_per_nucleon     = singleR_per_nucleon / singleR_per_nucleon_C12
    doubleR_per_nucleon_val = unumpy.nominal_values(doubleR_per_nucleon)
    doubleR_per_nucleon_err = unumpy.std_devs(doubleR_per_nucleon)

    doubleR_per_nucleus     = singleR_per_nucleus / singleR_per_nucleus_C12
    doubleR_per_nucleus_val = unumpy.nominal_values(doubleR_per_nucleus)
    doubleR_per_nucleus_err = unumpy.std_devs(doubleR_per_nucleus)
    '''


    # read target varibale (isolate only for a kin setting, since double SRC/MF ratios being taken)
    targ = np.array(df[(df['kin']=='SRC')]['target'])

    # loop over each target (to write numerical values to file)
    for i in np.arange(len(targ)):

        ofile.write('%s,%.3E,%.3E,%.3E,%.3E,%.3E,%.3E,%.3f,%.3E,%.1f,%.1f,%.1f,%.3f,%.3f\n' % (targ[i], singleR_A_c12_mf_val[i], singleR_A_c12_mf_err[i], singleR_A_c12_src_val[i], singleR_A_c12_src_err[i], singleR_per_proton_val[i], singleR_per_proton_err[i], doubleR_val[i], doubleR_err[i], N[i], Z[i], A[i], NoZ[i], NmZoA[i]) ) 

        
    ofile.close()




    
make_final_summary()

