import numpy as np
import pandas as pd
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

April 5, 2023: Modified
'''

#________________________
def get_cntmEff(run):
    # method: get Ca48 contamination efficinecy factor on a run-basis
    
    # set filename
    fname='special_studies/ca48_contamination/ca48_correction.csv'

    # read .csv file
    df = pd.read_csv(fname, comment='#')

    # define condition
    cond = df['run']==run

    # get absolute contamination (in percent) -- to use this, requires that there is not space between comma-sepatated values
    absContamCalc = df[cond].C_absCntm_calc    # calculated from fit
    absContamMeas = df[cond].C_absCntm_meas    # measured (only for Ca48 MF)

    contam_eff_calc = 1. -  (absContamCalc/100.)
    contam_eff_meas = 1. -  (absContamMeas/100.)

    # it depends on whether we want to apply calculated or measured contamination eff. for the Ca48 MF
    # for now we use the fit results 
    return contam_eff_calc  # efficiency factor

#________________________
def get_cntmEff(kin=''):
    # method: returns an array for contamination factor for either Ca48 MF or SRC
    
    # example of appending new row (to append corrected Ca48 yield)
    fname='special_studies/ca48_contamination/ca48_correction.csv'

    # read .csv file
    df = pd.read_csv(fname, comment='#')

    # define condition
    cond = df['kin']==kin

    # get absolute contamination (in percent) -- to use this, requires that there is not space between comma-sepatated values
    absContamCalc = df[cond].C_absCntm_calc    # calculated from fit
    absContamMeas = df[cond].C_absCntm_meas    # measured (only for Ca48 MF)

    contam_eff_calc = 1. -  (absContamCalc/100.)
    contam_eff_meas = 1. -  (absContamMeas/100.)

    # it depends on whether we want to apply calculated or measured contamination eff. for the Ca48 MF
    # for now we use the fit results 
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
    
    #output file to write summary file
    ofname = 'cafe_final_summary_pass2.csv' 
    ofile = open(ofname, 'w+')
    ofile.write('# CaFe Final Summary File (pass 2) \n')
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# target       : target name analyzed \n'
                '# kin          : kinematics analyzed  \n'
                '# beam_time    : beam-on-target time [s] \n'
                '# avg_current  : average beam current [uA] \n'
                '# total_charge : cumulative charge (over all runs) [mC] \n'
                '# yield        : yield (counts integrated over Pm) with all data-analysis cuts applied \n'
                '# yield_corr   : yield corrected for inefficiencies (hms/shms tracking + live_time) \n'
                '# yield_norm   : corrected yield (yield_corr) normalized by (total_charge transpacency area_density)\n'
                '# tgt_area_density: target density (g/cm2)\n'
                '# T N Z A    : transparency (T) # of neutrons (N) protons(Z) and nucleons (A)\n'
                )
    ofile.write('target,kin,beam_time,avg_current,total_charge,yield,yield_err,yield_corr,yield_corr_err,yield_norm,yield_norm_err,tgt_area_density,T,N,Z,A\n') 

    # target, kin list
    # target = ['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54']
    target = ['Be9', 'B10', 'B11', 'C12', 'Ca40',       'Ca48',   'Fe54', 'Au197']
    tcolor  = ['m',   'r',   'g',   'b',  'darkorange', 'dodgerblue', 'lime',  'gold'] 

    kin    = ['MF', 'SRC']


    # ------ create subplots for quality check plotting -----------

    
    #fig1_mf, ax1_mf  = plt.subplots(2, 3)
    fig1_src,ax1_src = plt.subplots(2, 3)
    
   # fig2, ax2 = plt.subplots()
   # fig3, ax3 = plt.subplots()
   # fig4, ax4 = plt.subplots()

    # -------------------------------------------------------------
    # loop over each target
    for idx in np.arange(len(target)):

        # loop over each kinmeatic (MF, SRC)
        for jdx in np.arange(len(kin)):

            df_empty_flag = 0
            
            # set generic summary file name
            summary_file_path = 'summary_files/pass2/cafe_prod_%s_%s_report_summary.csv' % (target[idx], kin[jdx])
            
            # read .csv file
            df = pd.read_csv(summary_file_path, comment='#')

            
            # ------- check if target is Ca48 (and read contamination factors) --------

            if(target[idx]=='Ca48'):
                
                # read array of contamination eff. factors
                cntm_eff = get_cntmEff(kin[jdx])
                
                if (kin[jdx]=='MF'):
                    cntm_eff = cntm_eff[-3:] # read contamination factor corresponding to last 3 Ca48 MF runs
                    df = df[-3:]  # select last 3 Ca48 MF runs from dataframe (esentially almost no contamination)

            # ----- END: check if target is Ca48 (and read contamination factors) --------

            
            # ------ read parameters -------
            T                 = find_param('transparency', summary_file_path) # transparency
            tgt_areal_density = find_param('target_areal_density', summary_file_path) # g/cm^2
            N                 = find_param('N:', summary_file_path) # number of neutrons
            Z                 = find_param('Z:', summary_file_path) # number of protons
            A                 = find_param('A:', summary_file_path) # number of nucleons
            # ---- END: read parameters -----

            
            # read selected data columns (with respective uncertainties)
            avg_current  = df['avg_current'] # average beam current [uA]     
            run          = df['run']
            charge       = df['charge'] # [mC]
            beam_time    = df['beam_time'] # (beam-on-target time) [sec]

            real_Yield   = unumpy.uarray(df['real_Yield'],      df['real_Yield_err']) 
            hms_trk_eff  = unumpy.uarray(df['hTrkEff'],         df['hTrkEff_err'])
            shms_trk_eff = unumpy.uarray(df['pTrkEff'],         df['pTrkEff_err'])
            total_LT     = unumpy.uarray(df['tLT'],             df['tLT_err_Bi'])
            #mult_trk_eff = unumpy.uarray(df['multi_track_eff'], df['multi_track_eff_err'])
            mult_trk_eff = np.array(df['multi_track_eff'])  # need to figure way to calculate uncertainty on this variable
            
            T1_scl_rate  = df['T1_scl_rate']  # SHMS  3/4    kHz
            T2_scl_rate  = df['T2_scl_rate']  # SHMS EL-REAL kHz
            T3_scl_rate  = df['T3_scl_rate']  # HMS 3/4      kHz
            T5_scl_rate  = df['T5_scl_rate']  # COIN         kHz

            # calculate other variables from data columns
            Qsum = charge.cumsum()  # calculate cumulative charge over all runs of given (target, kin) (may be helpful for plotting)
            T1_scl = T1_scl_rate * 1000 * beam_time   # total scaler counts
            T2_scl = T2_scl_rate * 1000 * beam_time
            T3_scl = T3_scl_rate * 1000 * beam_time
            T5_scl = T5_scl_rate * 1000 * beam_time 
   
            
            # sum over  counts (before applying any corrections, inefficinecy, charge, etc.) for book-keeping 
            real_yield_total = real_Yield.sum()

            
            # apply efficiency corrections to real yield (with uncertainties included) (run-by-run) and then sum over all counts
            real_Yield_corr       = real_Yield / (hms_trk_eff * shms_trk_eff * total_LT * mult_trk_eff)              
            real_Yield_corr_total = real_Yield_corr.sum()

            # --- apply Ca48 contamination correction ----
            if(target[idx]=='Ca48'):
                real_Yield_cntm_corr = real_Yield_corr / cntm_eff              
                real_Yield_cntm_corr_total = real_Yield_cntm_corr.sum()
                
                print('real_Yield_corr_total=',real_Yield_corr_total)
                print('real_Yield_cntm_corr_total=',real_Yield_cntm_corr_total)
                print('cntm_eff = ', cntm_eff)
            
            # sum over all total charge
            total_charge = charge.sum()

            # total beam time
            total_beam_time = beam_time.sum()

            # total average current
            total_avg_current = np.average(avg_current)
            
            # normalize yield by  total charge, transparency and target density
            yield_norm_per_run = real_Yield_corr / (charge * T * tgt_areal_density)  # counts / (mC * g/cm^2)
            yield_norm =  real_Yield_corr_total / (total_charge * T * tgt_areal_density)  # counts / (mC * g/cm^2)

            # normalize  Ca48 contamination-corrected yields 
            if(target[idx]=='Ca48'):
                yield_norm_cntm =  real_Yield_cntm_corr_total / (total_charge * T * tgt_areal_density)  # counts / (mC * g/cm^2)

                
            # ----------- make plots ----------------
            minT2 = T2_scl_rate==min(T2_scl_rate) #condition of minimum scaler rate
            minI  = avg_current==min(avg_current)
            if((kin[jdx]=='MF')):

                ax1_src[0, 0].errorbar(T2_scl_rate, unumpy.nominal_values(shms_trk_eff),  yerr=unumpy.std_devs(shms_trk_eff), marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                #y_loc = 0.985 - idx/500
                #plt.text(150, y_loc, '%s'%(target[idx]), color=tcolor[idx], fontsize = 15)

                ax1_src[0, 1].plot(T2_scl_rate, mult_trk_eff, marker='o', mec='k', color=tcolor[idx], linestyle='None' )

                ax1_src[0, 2].errorbar(T3_scl_rate, unumpy.nominal_values(hms_trk_eff),  yerr=unumpy.std_devs(hms_trk_eff), marker='o', mec='k', color=tcolor[idx], linestyle='None' )

                ax1_src[1, 0].errorbar(T2_scl_rate, unumpy.nominal_values(total_LT),  yerr=unumpy.std_devs(total_LT), marker='o', mec='k', color=tcolor[idx], linestyle='None' )

                # plot charge-norm yield vs. T2 rates
                ax1_src[1, 1].errorbar(T2_scl_rate, unumpy.nominal_values( yield_norm_per_run),  yerr=unumpy.std_devs( yield_norm_per_run), marker='o', mec='k', color=tcolor[idx], linestyle='None' )

                #ax1_src[1, 2].errorbar(T2_scl_rate, unumpy.nominal_values( yield_norm_per_run)/ unumpy.nominal_values( yield_norm_per_run)[minT2],  yerr=unumpy.std_devs( yield_norm_per_run)/ unumpy.nominal_values( yield_norm_per_run)[minT2], marker='o', mec='k', color=tcolor[idx], linestyle='None' )
                ax1_src[1, 2].errorbar(avg_current, unumpy.nominal_values( yield_norm_per_run)/ unumpy.nominal_values( yield_norm_per_run)[minI],  yerr=unumpy.std_devs( yield_norm_per_run)/ unumpy.nominal_values( yield_norm_per_run)[minI], marker='o', mec='k', color=tcolor[idx], linestyle='None' )

                

            #----------------------------------------
     
            # Write numerical data to final summary file
            ofile.write("%s,%s,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.3f,%.1f,%.1f,%.1f\n" % (target[idx].strip(), kin[jdx].strip(), total_beam_time, total_avg_current, total_charge, real_yield_total.n, real_yield_total.s, real_Yield_corr_total.n, real_Yield_corr_total.s, yield_norm.n, yield_norm.s, tgt_areal_density, T, N, Z, A) )

            # write to contamination-corrected values of Ca48 to file
            if(target[idx]=='Ca48'):

                ca48_tgt = 'Ca48_corr'
                # Write numerical data to final summary file
                ofile.write("%s,%s,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.3f,%.1f,%.1f,%.1f\n" % (ca48_tgt.strip(), kin[jdx].strip(), total_beam_time, total_avg_current, total_charge, real_yield_total.n, real_yield_total.s, real_Yield_cntm_corr_total.n, real_Yield_cntm_corr_total.s, yield_norm_cntm.n, yield_norm_cntm.s, tgt_areal_density, T, N, Z, A) )

    plt.show()          
    ofile.close()

    
    # apply b10,11 -carbide corrections + ca48 purity corrections
    applyB4C_correction(ofname)
    
    write_double_ratio(ofname, 'double_ratios_pass2.csv')
    
#_____________________________________
def applyB4C_correction(ofname=''):

    # method: apply boron-10,11 -carbide target corrections + Ca48 purity corrections

    # B4C10, B4C11 targets need to be carbon-subtracted
    # Ca48 is only ~90.5 % pure, the remaining is Ca40 
    # Ca48 is 90.5% pure, assuming that the remainder is Ca40.  The 48Ca target is 1051 mg/cm2.
    # If the purity refers to number of atoms, then the purity by weight is 91.5%.  (We need to find out whether purity refers to mass or to number.)  


    # set filename to be read
    #ofname='cafe_final_summary_pass2.csv'
    
    # read final summary file
    df = pd.read_csv(ofname, comment='#')
    
    # get areal densities for carbon-subtraction
    # g/cm2
    b4c10_density = np.array(df[(df['target']=='B10')]['tgt_area_density'][:-1])
    b4c11_density = np.array(df[(df['target']=='B11')]['tgt_area_density'][:-1])
    c12_density   = np.array(df[(df['target']=='C12')]['tgt_area_density'][:-1])

    # -------------------------------------

    
    # --- set Ca48 purity by weight-----

    ca48_purity = 0.915
    
    # get areal densities for Ca40 subtraction from Ca48 (g/cm2)
    ca48_density = np.array(df[(df['target']=='Ca48')]['tgt_area_density'][:-1])
    ca40_density = np.array(df[(df['target']=='Ca40')]['tgt_area_density'][:-1])

    # get norm. yield for ca40, ca48
    ca48_mf_yield_norm      = np.array(df[(df['target']=='Ca48')  & (df['kin']=='MF') ]['yield_norm'])
    ca48_mf_yield_norm_err  = np.array(df[(df['target']=='Ca48')  & (df['kin']=='MF') ]['yield_norm_err'])

    ca48_src_yield_norm     = np.array(df[(df['target']=='Ca48')  & (df['kin']=='SRC') ]['yield_norm'])
    ca48_src_yield_norm_err = np.array(df[(df['target']=='Ca48')  & (df['kin']=='SRC') ]['yield_norm_err'])

    ca40_mf_yield_norm      = np.array(df[(df['target']=='Ca40')  & (df['kin']=='MF') ]['yield_norm'])
    ca40_mf_yield_norm_err  = np.array(df[(df['target']=='Ca40')  & (df['kin']=='MF') ]['yield_norm_err'])
    
    ca40_src_yield_norm     = np.array(df[(df['target']=='Ca40')  & (df['kin']=='SRC') ]['yield_norm'])
    ca40_src_yield_norm_err = np.array(df[(df['target']=='Ca40')  & (df['kin']=='SRC') ]['yield_norm_err'])

    # add value +/- error for easy uncertainty propagations
    ca48_mf_yield_norm  = ufloat(ca48_mf_yield_norm, ca48_mf_yield_norm_err)
    ca48_src_yield_norm = ufloat(ca48_src_yield_norm, ca48_src_yield_norm_err)

    ca40_mf_yield_norm  = ufloat(ca40_mf_yield_norm, ca40_mf_yield_norm_err)
    ca40_src_yield_norm = ufloat(ca40_src_yield_norm, ca40_src_yield_norm_err)
       
    
    # calculate amount of ca40 in ca48
    ca40_cntm = 1. - (ca48_purity * ca48_density) # g/cm2 (amount of ca40)

    # corrected ca48 target thickness
    ca48_density_corr =  ca48_density - ca40_cntm

    # apply density correction (ca48 corrected density) and  ca40 subtraction from ca48
    ca48_mf_corr  = (ca48_mf_yield_norm  * ca48_density_corr /  ca48_density) - (ca40_cntm /ca40_density) * (ca40_mf_yield_norm  * ca48_density_corr /  ca48_density)
    ca48_src_corr = (ca48_src_yield_norm * ca48_density_corr /  ca48_density) - (ca40_cntm /ca40_density) * (ca40_src_yield_norm * ca48_density_corr /  ca48_density)

    print('ca48 SRC/MF (before Ca40 subtraction): ',  ca48_src_yield_norm/ ca48_mf_yield_norm)
    print('ca48 SRC/MF (after Ca40 subtraction): ',  ca48_src_corr / ca48_mf_corr)
    
    # -------------------------------------

    
    # Avogadro's number (# atoms / mol)
    Na = 6.0221408e23  
    
    # isotopic molar mass (g/mol)
    mol_b10 = 10.0129369  # shouldnt it be number of nucleons A?
    mol_b11 = 11.0093052
    mol_c12 = 12.0107
    
    # boron-carbide (4 boron-10 or boron-11 atoms + 1 carbon atom) molar mass
    mol_b4c10 = 4*mol_b10 + mol_c12 
    mol_b4c11 = 4*mol_b11 + mol_c12 
    
    # number of target atoms / cm^2 (# of scatterers) = (atoms/mol) * g/cm2  / (g/mol)  
    targetfac_b4c10 = Na * b4c10_density  /  mol_b4c10
    targetfac_b4c11 = Na * b4c11_density  /  mol_b4c11
    targetfac_c12   = Na * c12_density    /  mol_c12
    
    # IMPORTANT: Recover the areal density (g/cm2) of B10, B11 (for normalizing after subtraction)
    b10_density = (b4c10_density / mol_b4c10) * 4. * 10.  # [(g/cm^2) / (g/mol)] * 4 (B10 atoms) * 10 g/mol
    b11_density = (b4c11_density / mol_b4c11) * 4. * 11.  # [(g/cm^2) / (g/mol)] * 4 (B11 atoms) * 11 g/mol

    '''
    print('areal densities')
    print('b4c10_density (g/cm2) =', b4c10_density)
    print('b4c11_density (g/cm2) =', b4c11_density)
    print('b10_density (g/cm2) =', b10_density)
    print('b11_density (g/cm2) =', b11_density)
    '''
    
    # set conditions for selecting b10, b11 from dataframe
    cond_b10_mf  =((df['target']=='B10') & (df['kin']=='MF'))
    cond_b10_src =((df['target']=='B10') & (df['kin']=='SRC'))
    cond_b11_mf  =((df['target']=='B11') & (df['kin']=='MF'))
    cond_b11_src =((df['target']=='B11') & (df['kin']=='SRC'))
    cond_c12_mf  =((df['target']=='C12') & (df['kin']=='MF'))
    cond_c12_src =((df['target']=='C12') & (df['kin']=='SRC'))
    
    # get the efficiency corrected yield (not yet normalized by total charge) and total charge for b4c10, b4c11, c12 (mf, src)
    
    # MF
    yield_corr_b4c10_mf     = unumpy.uarray(df[cond_b10_mf]['yield_corr'], df[cond_b10_mf]['yield_corr_err']) #val +/- err
    charge_b4c10_mf         = np.array(df[cond_b10_mf]['total_charge'])
    
    yield_corr_b4c11_mf     = unumpy.uarray(df[cond_b11_mf]['yield_corr'], df[cond_b11_mf]['yield_corr_err'])
    charge_b4c11_mf         = np.array(df[cond_b11_mf]['total_charge'])
    
    yield_corr_c12_mf       = unumpy.uarray(df[cond_c12_mf]['yield_corr'], df[cond_c12_mf]['yield_corr_err'])
    charge_c12_mf           = np.array(df[cond_c12_mf]['total_charge'])
    
    # SRC
    yield_corr_b4c10_src     = unumpy.uarray(df[cond_b10_src]['yield_corr'], df[cond_b10_src]['yield_corr_err']) #val +/- err
    charge_b4c10_src         = np.array(df[cond_b10_src]['total_charge'])
    
    yield_corr_b4c11_src     = unumpy.uarray(df[cond_b11_src]['yield_corr'], df[cond_b11_src]['yield_corr_err'])
    charge_b4c11_src         = np.array(df[cond_b11_src]['total_charge'])
    
    yield_corr_c12_src       = unumpy.uarray(df[cond_c12_src]['yield_corr'], df[cond_c12_src]['yield_corr_err'])
    charge_c12_src           = np.array(df[cond_c12_src]['total_charge'])
    
    
    # --------------------------------------------------------------------------
    # Apply carbon subtraction to charge-normalized b4c10, b4c11 (counts / mC)
    # --------------------------------------------------------------------------
    
    #--- MF ---
    yield_norm_b410_mf      = (yield_corr_b4c10_mf/charge_b4c10_mf) - ((yield_corr_c12_mf/charge_c12_mf) * targetfac_b4c10/targetfac_c12)
    yield_norm_b411_mf      = (yield_corr_b4c11_mf/charge_b4c11_mf) - ((yield_corr_c12_mf/charge_c12_mf) * targetfac_b4c11/targetfac_c12)

    '''
    print('Apply the Carbon subtraction to MF B4C10 and B4C11 Yield/charge')
    print('Boron-10 (MF):')
    print('B4C10_MF Yield/Q = %.3f'%(unumpy.nominal_values(yield_corr_b4c10_mf)/charge_b4c10_mf) )
    print('C12_MF Yield/Q= %.3f' %(unumpy.nominal_values(yield_corr_c12_mf)/charge_c12_mf))
    print('targetfac_b4c10/targetfac_c12 = %.3E' %(targetfac_b4c10/targetfac_c12))
    print('B410_MF Yield/Q = %.3f' % ( unumpy.nominal_values(yield_norm_b410_mf) ))
    print('')
    print('Boron-11 (MF):')
    print('B4C11_MF Yield/Q = %.3f'%(unumpy.nominal_values(yield_corr_b4c11_mf)/charge_b4c11_mf) )
    print('C12_MF Yield/Q= %.3f' %(unumpy.nominal_values(yield_corr_c12_mf)/charge_c12_mf))
    print('targetfac_b4c11/targetfac_c12 = %.3E' %(targetfac_b4c11/targetfac_c12))
    print('B411_MF Yield/Q = %.3f' % ( unumpy.nominal_values(yield_norm_b411_mf) ))
    print('---------------------------------\n')
    '''
    
    # apply the remaining scale factors (transparency and boron 10, 11 densities)
    yield_norm_b10_mf   = (yield_norm_b410_mf / 4.) / ( df[cond_b10_mf]['T'] * b10_density) 
    yield_norm_b11_mf   = (yield_norm_b411_mf / 4.) / ( df[cond_b11_mf]['T'] * b11_density)

    '''
    print('Apply the remaining scaler factors (transparency and areal density) -- not that it matters for double ratio, but just to check')
    print('B10_MF Yield/Q = (B410_MF Yield/Q) / 4 B10 atoms = %.3f' % (unumpy.nominal_values(yield_norm_b410_mf) / 4) )
    print('T, B10_density (g/cm2) = %.3f, %.3f' % (df[cond_b10_mf]['T'],  b10_density) )
    print('B10_MF Yield/(Q * T * B10_density) = %.3f' % ( unumpy.nominal_values(yield_norm_b10_mf) ))
    print('')
    print('B11_MF Yield/Q = (B411_MF Yield/Q) / 4 B11 atoms = %.3f' % (unumpy.nominal_values(yield_norm_b411_mf) / 4) )
    print('T, B11_density (g/cm2) = %.3f, %.3f' % (df[cond_b11_mf]['T'],  b11_density) )
    print('B11_MF Yield/(Q * T * B11_density) = %.3f' % ( unumpy.nominal_values(yield_norm_b11_mf) ))
    print('---------------------------------\n')
    '''
    
    #--- SRC ---
    yield_norm_b410_src      = (yield_corr_b4c10_src/charge_b4c10_src) - (yield_corr_c12_src/charge_c12_src) * targetfac_b4c10/targetfac_c12
    yield_norm_b411_src      = (yield_corr_b4c11_src/charge_b4c11_src) - (yield_corr_c12_src/charge_c12_src) * targetfac_b4c11/targetfac_c12

    '''
    print('Apply the Carbon subtraction to SRC B4C10 and B4C11 Yield/charge')
    print('Boron-10 SRC:')
    print('B4C10_SRC Yield/Q = %.3f'%(unumpy.nominal_values(yield_corr_b4c10_src)/charge_b4c10_src) )
    print('C12_SRC Yield/Q= %.3f' %(unumpy.nominal_values(yield_corr_c12_src)/charge_c12_src))
    print('targetfac_b4c10/targetfac_c12 = %.3E' %(targetfac_b4c10/targetfac_c12))
    print('B410_SRC Yield/Q = %.3f' % ( unumpy.nominal_values(yield_norm_b410_src) ))
    print('')
    print('Boron-11 SRC:')
    print('B4C11_SRC Yield/Q = %.3f'%(unumpy.nominal_values(yield_corr_b4c11_src)/charge_b4c11_src) )
    print('C12_SRC Yield/Q= %.3f' %(unumpy.nominal_values(yield_corr_c12_src)/charge_c12_src))
    print('targetfac_b4c11/targetfac_c12 = %.3E' %(targetfac_b4c11/targetfac_c12))
    print('B411_SRC Yield/Q = %.3f' % ( unumpy.nominal_values(yield_norm_b411_src) ))
    print('---------------------------------\n')
    '''
    
    # apply the remaining scale factors (transparency and boron 10, 11 densities)
    yield_norm_b10_src   = (yield_norm_b410_src / 4. ) / ( df[cond_b10_src]['T'] * b10_density) 
    yield_norm_b11_src   = (yield_norm_b411_src / 4. ) / ( df[cond_b11_src]['T'] * b11_density)

    '''
    print('Apply the remaining scaler factors (transparency and areal density) -- not that it matters for double ratio, but just to check')
    print('B10_SRC Yield/Q = (B410_SRC Yield/Q) / 4 B10 atoms = %.3f' % (unumpy.nominal_values(yield_norm_b410_src) / 4) )
    print('T, B10_density (g/cm2) = %.3f, %.3f' % (df[cond_b10_src]['T'],  b10_density) )
    print('B10_SRC Yield/(Q * T * B10_density) = %.3f' % ( unumpy.nominal_values(yield_norm_b10_src) ))
    print('')
    print('B11_SRC Yield/Q = (B411_SRC Yield/Q) / 4 B11 atoms = %.3f' % (unumpy.nominal_values(yield_norm_b411_src) / 4) )
    print('T, B11_density (g/cm2) = %.3f, %.3f' % (df[cond_b11_src]['T'],  b11_density) )
    print('B11_SRC Yield/(Q * T * B11_density) = %.3f' % ( unumpy.nominal_values(yield_norm_b11_src) ))
    print('---------------------------------\n')
    '''

      
    # recover yield_corr (correctded for inefficiencies, but not yet normalized by Q, T, area_density)
    yield_corr_b10_mf = yield_norm_b10_mf * (b10_density * df[cond_b10_mf]['T'] * df[cond_b10_mf]['total_charge'])
    yield_corr_b11_mf = yield_norm_b11_mf * (b11_density * df[cond_b11_mf]['T'] * df[cond_b11_mf]['total_charge'])

    yield_corr_b10_src = yield_norm_b10_src * (b10_density * df[cond_b10_src]['T'] * df[cond_b10_src]['total_charge'])
    yield_corr_b11_src = yield_norm_b11_src * (b11_density * df[cond_b11_src]['T'] * df[cond_b11_src]['total_charge'])

    # recover back uncorrected yield from inefficiencies
    # Cant do it. The efficiencies were only read in and applied in the previous method. Its ok, This is not as important.

    # this is how to extract the nominal values and error from unumpy (uncertainty package)
    # unumpy.nominal_values(yield_corr_b10_mf),  unumpy.std_devs(yield_corr_b10_mf)
    print('df:', df)
    print('df[2:6]:', df[2:6])

    # add B10, B11 rows to write applied corrections
    df = pd.concat([df,df[2:6]], ignore_index=True)

   
    
    # last index (we need last four indices to update)
    idx_last = df.tail(1).index.item()

    print('unumpy.nominal_values(yield_corr_b10_mf)=',unumpy.nominal_values(yield_corr_b10_mf))
    #Updating B10 MF
    df.loc[idx_last-3, ['target']] = ['B10_corr']
    df.loc[idx_last-3, ['yield_corr']]       = round(unumpy.nominal_values(yield_corr_b10_mf)[0], 3)
    df.loc[idx_last-3, ['yield_corr_err']]   = round(unumpy.std_devs(yield_corr_b10_mf)[0], 3)
    df.loc[idx_last-3, ['yield_norm']]       = round(unumpy.nominal_values(yield_norm_b10_mf)[0], 3)
    df.loc[idx_last-3, ['yield_norm_err']]   = round(unumpy.std_devs(yield_norm_b10_mf)[0], 3)
    df.loc[idx_last-3, ['tgt_area_density']] = round(b10_density[0], 3)
    #Updating B10 SRC
    df.loc[idx_last-2, ['target']] = ['B10_corr']
    df.loc[idx_last-2, ['yield_corr']]       = round(unumpy.nominal_values(yield_corr_b10_src)[0], 3)
    df.loc[idx_last-2, ['yield_corr_err']]   = round(unumpy.std_devs(yield_corr_b10_src)[0], 3)
    df.loc[idx_last-2, ['yield_norm']]       = round(unumpy.nominal_values(yield_norm_b10_src)[0], 3)
    df.loc[idx_last-2, ['yield_norm_err']]   = round(unumpy.std_devs(yield_norm_b10_src)[0], 3)
    df.loc[idx_last-2, ['tgt_area_density']] = round(b10_density[0], 3)
    #Updating B11 MF
    df.loc[idx_last-1, ['target']] = ['B11_corr']
    df.loc[idx_last-1, ['yield_corr']]       = round(unumpy.nominal_values(yield_corr_b11_mf)[0], 3)
    df.loc[idx_last-1, ['yield_corr_err']]   = round(unumpy.std_devs(yield_corr_b11_mf)[0], 3)
    df.loc[idx_last-1, ['yield_norm']]       = round(unumpy.nominal_values(yield_norm_b11_mf)[0], 3)
    df.loc[idx_last-1, ['yield_norm_err']]   = round(unumpy.std_devs(yield_norm_b11_mf)[0], 3)
    df.loc[idx_last-1, ['tgt_area_density']] = round(b11_density[0], 3)
    #Updating B11 SRC
    df.loc[idx_last, ['target']] = ['B11_corr']
    df.loc[idx_last, ['yield_corr']]       = round(unumpy.nominal_values(yield_corr_b11_src)[0], 3)
    df.loc[idx_last, ['yield_corr_err']]   = round(unumpy.std_devs(yield_corr_b11_src)[0], 3)
    df.loc[idx_last, ['yield_norm']]       = round(unumpy.nominal_values(yield_norm_b11_src)[0], 3)
    df.loc[idx_last, ['yield_norm_err']]   = round(unumpy.std_devs(yield_norm_b11_src)[0], 3)
    df.loc[idx_last, ['tgt_area_density']] = round(b11_density[0], 3)

    # Updating Ca48_corr MF 
    # get index value for updating entry
    idx_ca48_corr_mf = df.index[(df['target']=='Ca48_corr') & (df['kin']=='MF')][0]
    df.loc[idx_ca48_corr_mf, ['yield_norm']]       = round(unumpy.nominal_values(ca48_mf_corr)[0], 3)
    df.loc[idx_ca48_corr_mf, ['yield_norm_err']]   = round(unumpy.std_devs(ca48_mf_corr)[0], 3)  
    df.loc[idx_ca48_corr_mf, ['tgt_area_density']] = round(ca48_density_corr[0], 3)
    
    #Updating Ca48_corr SRC 
    # get index value for updating entry
    idx_ca48_corr_src = df.index[(df['target']=='Ca48_corr') & (df['kin']=='SRC')][0]
    df.loc[idx_ca48_corr_src, ['yield_norm']]       = round(unumpy.nominal_values(ca48_src_corr)[0], 3)
    df.loc[idx_ca48_corr_src, ['yield_norm_err']]   = round(unumpy.std_devs(ca48_src_corr)[0], 3)  
    df.loc[idx_ca48_corr_src, ['tgt_area_density']] = round(ca48_density_corr[0], 3)
    
    df.to_csv('cafe_final_summary_pass2.csv', index=False)
    

    

def write_double_ratio(ifname='', ofname=''):

    # method: reads input file ifname with (cafe final summary file),
    #         write output file ofname with double-ratio numerical value
    #
    #         double_ratio = ( A_SRC / A_MF ) / ( C12_SRC / C12_MF )
    
    
    # read input final summary file
    df = pd.read_csv(ifname, comment='#')


    # set output file to write double ratio numerical values
    ofile = open(ofname, 'w+')
    ofile.write('# CaFe numerical double ratios (pass 2) \n')
    ofile.write('# \n'
                '# Header Definitions: \n'
                '# target       : target A used in double ratio \n'
                '# doubleR      : double ratio of target A(SRC/MF) relative to C12 (SRC/MF) \n'
                '# doubleR_err  : uncertainty in double ratio \n'
                '# N: Z: A      : # of neutrons (N): protons(Z): nucleons (A): for target A \n'
                '# NoZ          : N/Z \n'
                '# NmZoA        : (N-Z)/A \n'                
                )
    ofile.write('target,doubleR,doubleR_err,N,Z,A,NoZ,NmZoA\n') 

    
    src_yield_norm = df[(df['kin']=='SRC')]['yield_norm']
    src_yield_norm_err = df[(df['kin']=='SRC')]['yield_norm_err']

    mf_yield_norm = df[(df['kin']=='MF')]['yield_norm']
    mf_yield_norm_err = df[(df['kin']=='MF')]['yield_norm_err']

    src_yield_norm_C12 = df[(df['kin']=='SRC') & (df['target']=='C12')]['yield_norm']
    src_yield_norm_err_C12 = df[(df['kin']=='SRC') & (df['target']=='C12')]['yield_norm_err']

    mf_yield_norm_C12 = df[(df['kin']=='MF') & (df['target']=='C12')]['yield_norm']
    mf_yield_norm_err_C12 = df[(df['kin']=='MF') & (df['target']=='C12')]['yield_norm_err']

    print('src_yield_norm=',src_yield_norm)
    # put into arrays for error calculation
    src_yield_norm_arr  = unumpy.uarray(src_yield_norm,  src_yield_norm_err)
    mf_yield_norm_arr   = unumpy.uarray(mf_yield_norm,  mf_yield_norm_err)

    src_yield_norm_C12_arr  = unumpy.uarray(src_yield_norm_C12, src_yield_norm_err_C12)
    mf_yield_norm_C12_arr   = unumpy.uarray(mf_yield_norm_C12,  mf_yield_norm_err_C12)

    # calculate double-ratio
    double_ratio = (src_yield_norm_arr/mf_yield_norm_arr) / (src_yield_norm_C12_arr/mf_yield_norm_C12_arr )

    double_ratio_val = unumpy.nominal_values(double_ratio)
    double_ratio_err = unumpy.std_devs(double_ratio)
    
    N = np.array(df[(df['kin']=='SRC')]['N'])
    Z = np.array(df[(df['kin']=='SRC')]['Z'])
    A = np.array(df[(df['kin']=='SRC')]['A'])
    NoZ = N/Z
    NmZoA = (N-Z) / A

    # read target varibale (isolate only for a kin setting, since double SRC/MF ratios being taken)
    targ = np.array(df[(df['kin']=='SRC')]['target'])

    # loop over each target (to write numerical values to file)
    for i in np.arange(len(targ)):

        ofile.write('%s,%.3f,%.3E,%.1f,%.1f,%.1f,%.3f,%.3f\n' % (targ[i], double_ratio_val[i], double_ratio_err[i], N[i], Z[i], A[i], NoZ[i], NmZoA[i]) ) 

        
    ofile.close()
    
    '''
    
    #          'LD2', 'Be9', 'B10',   'B11',   'C12', 'Ca40',       'Ca48',  'Ca48_corr' , 'Fe54', 'B10_corr', 'B11_corr'
    # tcolor  = ['c',    'm',   'r',     'g',     'b',  'darkorange', 'violet', 'violet',    'gold',    'r',      'g'     ] 

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

    
make_final_summary()

