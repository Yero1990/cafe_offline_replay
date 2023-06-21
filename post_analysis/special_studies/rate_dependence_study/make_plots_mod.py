import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy


def get_data(header='', fname=''):

    # read any data from the .csv file
    
    # read .csv file
    df = pd.read_csv(fname, comment='#')
    
    # sort from small to larger beam current (for normalizing purposes) 
    df = df.sort_values(by=['avg_current'])

    #print(df)
    data_array = np.array(df[header])

    return data_array

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


def get_relative(name='',  file1=''):

    # brief: function returns the relative charger-normzliced
    #        efficiency-corrected yield and the charge normalized
    #        T2 scalers for a given .csv file input

    # inputs:
    # name: 'rel_yield_nom', 'rel_yield_err', 'rel_T2_nom', 'rel_T2_err' # what the user wants to extract
    # file1:  .csv summary file name to read
    
    # get relevant header info (from file1)
    I                    = get_data('avg_current', file1 )
    Q                    = get_data('charge', file1)
    T2_scl               = get_data('T2_scl_rate', file1) * 1000. * get_data('beam_time', file1) #[kHz] * 1000 Hz/1kHz [sec] = [counts] 

    counts                    = unumpy.uarray(get_data('real_Yield', file1),  get_data('real_Yield_err', file1) )
    hms_trk_eff   = unumpy.uarray(get_data('hTrkEff',    file1),  get_data('hTrkEff_err',    file1) )
    shms_trk_eff  = unumpy.uarray(get_data('pTrkEff',    file1),  get_data('pTrkEff_err',    file1) )
    total_LT      = unumpy.uarray(get_data('tLT',        file1),  get_data('tLT_err_Bi',     file1) )
    mult_trk_eff  = np.array(get_data('multi_track_eff', file1) )

    # ------ read parameters -------

    tgt_thick         = find_param('target_areal_density', file1) # g/cm^2
    tgt_thick_corr    = tgt_thick  # set corrected thickness to thickness (will be re-defined if impurity is corrected for any target)
    N                 = find_param('N:', file1) # number of neutrons
    Z                 = find_param('Z:', file1) # number of protons
    A                 = find_param('A:', file1) # number of nucleons

    # Transparency function: T = c * A ** alpha (Q2), where alpha ~ -0.24 for Q2 >= 2 GeV^2, and c=1, A -> mass number
    # reference: https://arxiv.org/abs/1211.2826  "Color Transparency: past, present and future"
    alpha=-0.24
    T = A**(alpha)
    
    # calcualte charge-normalized yield
    Y     = counts / ( Q * hms_trk_eff * shms_trk_eff * total_LT * mult_trk_eff * T * tgt_thick)              
    relY  = Y / Y[0] 

    Y_nom = unumpy.nominal_values(Y)
    Y_err = unumpy.std_devs(Y)

    relY_nom = unumpy.nominal_values(relY)
    relY_err = unumpy.std_devs(relY)
    
    #calculate T2 scaler counts / charge
    T2    = T2_scl/Q 
    relT2 = T2/T2[0]

    T2_nom = unumpy.nominal_values(T2)
    T2_err = unumpy.std_devs(T2)
    
    relT2_nom = unumpy.nominal_values(relT2)
    relT2_err = unumpy.std_devs(relT2)

    if(name=='rel_yield_nom'):
        return relY_nom
    elif(name=='rel_yield_err'):
        return relY_err
    elif(name=='rel_T2_nom'):
        return relT2_nom
    elif(name=='rel_T2_err'):
        return relT2_err

def compare():
    
    # Brief: this method is used to compare changes in the relative charge-norm yield
    # for a given target, for various studies done. Additional studies can be added, bu
    # then the overlay can get a little bit crowded.  converntion: use file1 for new, and file2 for existing

    # targets: Be9, B10, B11, C12, Ca40, Ca48, Fe54, Au197
    # the phases represent progressions in analysis: phase0 -> baseline, phase1 -> coin. time window cut, phase2-> phase1+tighter ref time cut, phase3 -> phase2+ evt_type>=4 cut


    #--------------------------
    
    # define output file to write correction factors (change in relative yield vs. T2 rates)

    # output file to write summary file
    ofname = 'cafe_relYield_corrections.csv'
    ofile = open(ofname, 'w+')
    ofile.write('# CaFe Relative Yield Correction Factors (using MF kinematics) \n')
    ofile.write('# \n'
                '# Header Definitions: \n'
                #'# phase#       : 0 -> baseline (pruning)   \n'
                #'# phase#       : 1 -> 0 + added the T2,T3 TDC raw time window cuts (tcoin param)   \n'
                #'# phase#       : 2 -> 1 + added tighter reference time cuts on HMS/SHMS to try and recover more lost coincidences    \n'
                #'# phase#       : 3 -> 2 + added an event type cut (g.evtyp>=4) to properly select ALL coincidences   \n'
                '# deltaY       : difference in relative charge-norm, eff.-corrected data yield between lowest and highest T2 rate runs \n'
                '# deltaY_err   : error in deltaY \n'
                '# deltaX       : difference between lowest and highest T2 rates [kHz] \n'
                '# slope        : deltaY/deltaX or drop in relative yield per kHz (to be used as correction factor) \n'
                '# slope_err    : error in slope \n'
                )
    ofile.write('phase,target,deltaY,deltaY_err,deltaX,slope,slope_err\n') 

    target = ['Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54', 'Au197']

    

    # loop over each target
    for t in np.arange(len(target)):

        # read the summary files for each study phase   (only for pass 1)
        # file_arr= ['../../summary_files/rate_dependence_study/phase0/cafe_prod_%s_MF_report_summary.csv'%(target[t]),
        #           '../../summary_files/rate_dependence_study/phase1/cafe_prod_%s_MF_report_summary.csv'%(target[t]),
        #           '../../summary_files/rate_dependence_study/phase2/cafe_prod_%s_MF_report_summary.csv'%(target[t]),
        #           '../../summary_files/rate_dependence_study/phase3/cafe_prod_%s_MF_report_summary.csv'%(target[t]) ]

        # imarker = ['o', 's', '^', 'v']
        # icolor  = ['k', 'gray', 'b', 'g']
        
        # pass 2
        #file_arr = ['../../summary_files/rate_dependence_study/phase0/cafe_prod_%s_MF_report_summary.csv' % (target[t]), '../../summary_files/rate_dependence_study/pass3/cafe_prod_%s_MF_report_summary.csv' % (target[t])]
        file_arr = ['../../summary_files/rate_dependence_study/phase0/cafe_prod_%s_MF_report_summary.csv' % (target[t])]

        imarker = ['s', 'o']
        icolor  = ['gray', 'red']

        # figure to plot relative yields vs. avg current and T2 rates
        fig1, axs1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    
        # loop over each phase study
        for i, ifile in enumerate(file_arr):

            
            print('phase: ', i)
            # get relevant header info (from ifile)
            I                    = get_data('avg_current', ifile )
            Q                    = get_data('charge', ifile)
            T2_scl_rate          = get_data('T2_scl_rate', ifile)
            T2_scl               = get_data('T2_scl_rate', ifile) * 1000. * get_data('beam_time', ifile) #[kHz] * 1000 Hz/1kHz [sec] = [counts] 
            counts                    = unumpy.uarray(get_data('real_Yield', ifile),  get_data('real_Yield_err', ifile) )
            hms_trk_eff   = unumpy.uarray(get_data('hTrkEff',    ifile),  get_data('hTrkEff_err',    ifile) )
            shms_trk_eff  = unumpy.uarray(get_data('pTrkEff',    ifile),  get_data('pTrkEff_err',    ifile) )
            total_LT      = unumpy.uarray(get_data('tLT',        ifile),  get_data('tLT_err_Bi',     ifile) )
            mult_trk_eff  = np.array(get_data('multi_track_eff', ifile) )

            # ------ read parameters -------
            N                 = find_param('N:', ifile) # number of neutrons
            Z                 = find_param('Z:', ifile) # number of protons
            A                 = find_param('A:', ifile) # number of nucleons
            tgt_thick         = find_param('target_areal_density', ifile) # g/cm^2
            
            minT2 = T2_scl_rate==min(T2_scl_rate) #condition of minimum scaler rate
            print(minT2)
            print('file_arr:',  file_arr[i])
            # Transparency function: T = c * A ** alpha (Q2), where alpha ~ -0.24 for Q2 >= 2 GeV^2, and c=1, A -> mass number
            # reference: https://arxiv.org/abs/1211.2826  "Color Transparency: past, present and future"
            alpha=-0.24
            T = A**(alpha)
    
            # calcualte charge-normalized yield
            Y     = counts / ( Q * hms_trk_eff * shms_trk_eff * total_LT * mult_trk_eff * T * tgt_thick)              
            relY  = Y / Y[minT2] 
            
            Y_nom = unumpy.nominal_values(Y)
            Y_err = unumpy.std_devs(Y)
            
            relY_nom = unumpy.nominal_values(relY)
            relY_err = unumpy.std_devs(relY)
            
            #calculate T2 scaler counts / charge
            T2    = T2_scl/Q 
            relT2 = T2/T2[minT2]
            
            T2_nom = unumpy.nominal_values(T2)
            T2_err = unumpy.std_devs(T2)
            
            relT2_nom = unumpy.nominal_values(relT2)
            relT2_err = unumpy.std_devs(relT2)
            
            # calculate the  difference in relative yield (y-axis) and difference in rates (x-axis)
            # this gives the slope:  deltaY / deltaX which gives a correction factor
            deltaY = max(relY) - min(relY)
            deltaY_nom = unumpy.nominal_values(deltaY)  # fractional form (need to multiply by 100 to convert to %)
            deltaY_err = unumpy.std_devs(deltaY)        # fractional form (need to multiply by 100 to convert to %)
            deltaX = max(T2_scl_rate) - min(T2_scl_rate)  # kHz (range of T2 scalers)
            slope= deltaY_nom / deltaX
            slope_err = deltaY_err / deltaX

            # write numerical values to file
            ofile.write("%i, %s, %.3E, %.3E, %.3f, %.3E, %.3E \n" % (i, target[t], deltaY_nom, deltaY_err, deltaX, slope, slope_err))
            
            print('-----------------------------\n')
            print('target=', target)
            print('relY_nom=',relY_nom)
            print('relY_err=',relY_err)
            print('T2_scl_rate=', T2_scl_rate)
            print('deltaY_nom =', deltaY_nom)
            print('deltaY_err =', deltaY_err)
            print('deltaX [kHz]=',deltaX)
            print('slope=', slope)
            print('slope_err=', slope_err)
            print('-----------------------------\n')
            
            
            total_LT_nom = unumpy.nominal_values(total_LT )
            total_LT_err = unumpy.std_devs(total_LT )

            # ---------------------------------------------------------------
            #  PLOT relative yileds vs (avg current or rate) for each phase
            # ---------------------------------------------------------------
            
            #ax1 = axs1[0]
            #ax1.errorbar(I, relY_nom, relY_err,   marker=imarker[i], markersize=8, color=icolor[i], mec='k', linestyle='dashed')
            #ax1.set_ylabel('Relative Yield', fontsize=18)
             #ax1.set_xlabel('Average Current [uA]', fontsize=18)
            #ax1.tick_params(axis='x', labelsize=14)
            #ax1.tick_params(axis='y', labelsize=14)
            #ax1.grid(True)

                       
            axs1.errorbar(  T2_scl_rate , relY_nom, relY_err,   marker=imarker[i], markersize=8, color=icolor[i], mec='k', linestyle='dashed', label='%s MF'%target[t])            
            axs1.set_xlabel('T2 Scaler Rate [kHz]', fontsize=18)
            axs1.set_ylabel('Relative Yield', fontsize=18)
            axs1.tick_params(axis='x', labelsize=14)
            axs1.tick_params(axis='y', labelsize=14)
            axs1.grid(True)
            
            fig1.tight_layout() 

            #----------------------------------------------------------------

            
            
        #plt.legend()
        plt.show()


def plot_correction():

    # input file to read relative yield corr. factors
    ifname = 'cafe_relYield_corrections_pass2.csv'
    
    # read .csv file
    df = pd.read_csv(ifname, comment='#')
    
    # require phase 3 (filter data frame) -- pass1
    #df_phase = df[df['phase']==3]

    # require phase 0 (filter data frame) -- pass2
    df_phase = df[df['phase']==0]
    
    fig, axs = plt.subplots(2)
    fig.suptitle('Rate Dependence Correction')
    axs[0].errorbar(df_phase['deltaX'], df_phase['deltaY'], df_phase['deltaY_err'], marker='o', markersize=8, linestyle='None', color='g', mec='k')
    axs[0].set_ylabel(r'$\Delta$Y', fontsize=18)
    axs[0].set_xlabel(r'$\Delta$X [kHz]', fontsize=18)
    axs[0].tick_params(axis='x', labelsize=14)
    axs[0].tick_params(axis='y', labelsize=14)
    axs[0].grid(True)

    tgt = df_phase['target']
    dX = np.array(df_phase['deltaX'])
    dY = np.array(df_phase['deltaY'])

   
    axs[1].errorbar(tgt, df_phase['slope'], df_phase['slope_err'], marker='o', markersize=8, linestyle='None', color='g', mec='k')
    axs[1].set_ylabel(r'slope', fontsize=18)
    axs[1].set_xlabel(r'target', fontsize=18)
    axs[1].tick_params(axis='x', labelsize=14)
    axs[1].tick_params(axis='y', labelsize=14)
    axs[1].grid(True)

    print(tgt)
    fig.tight_layout() 
    plt.show()


# calling function to compare different phases of the rate-dependece study
compare()

#plot_correction()
