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

    N                    = unumpy.uarray(get_data('real_Yield', file1),  get_data('real_Yield_err', file1) )
    hms_trk_eff   = unumpy.uarray(get_data('hTrkEff',    file1),  get_data('hTrkEff_err',    file1) )
    shms_trk_eff  = unumpy.uarray(get_data('pTrkEff',    file1),  get_data('pTrkEff_err',    file1) )
    total_LT      = unumpy.uarray(get_data('tLT',        file1),  get_data('tLT_err_Bi',     file1) )
    mult_trk_eff  = np.array(get_data('multi_track_eff', file1) )
    
    # calcualte charge-normalized yield
    Y     = N / ( Q * hms_trk_eff * shms_trk_eff * total_LT * mult_trk_eff)              
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
    
    target='Be9'
    file_arr= ['../../summary_files/rate_dependence_study/phase0/cafe_prod_%s_MF_report_summary.csv'%(target),
               '../../summary_files/rate_dependence_study/phase1/cafe_prod_%s_MF_report_summary.csv'%(target),
               '../../summary_files/rate_dependence_study/phase2/cafe_prod_%s_MF_report_summary.csv'%(target),
               '../../summary_files/rate_dependence_study/phase3/cafe_prod_%s_MF_report_summary.csv'%(target) ]

    imarker = ['o', 's', '^', 'v']
    icolor  = ['k', 'gray', 'b', 'g']
    
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))

    for i, ifile in enumerate(file_arr):

        print('phase: ', i)
        # get relevant header info (from ifile)
        I                    = get_data('avg_current', ifile )
        Q                    = get_data('charge', ifile)
        T2_scl_rate          = get_data('T2_scl_rate', ifile)
        T2_scl               = get_data('T2_scl_rate', ifile) * 1000. * get_data('beam_time', ifile) #[kHz] * 1000 Hz/1kHz [sec] = [counts] 
        N                    = unumpy.uarray(get_data('real_Yield', ifile),  get_data('real_Yield_err', ifile) )
        hms_trk_eff   = unumpy.uarray(get_data('hTrkEff',    ifile),  get_data('hTrkEff_err',    ifile) )
        shms_trk_eff  = unumpy.uarray(get_data('pTrkEff',    ifile),  get_data('pTrkEff_err',    ifile) )
        total_LT      = unumpy.uarray(get_data('tLT',        ifile),  get_data('tLT_err_Bi',     ifile) )
        mult_trk_eff  = np.array(get_data('multi_track_eff', ifile) )
        
        # calcualte charge-normalized yield
        Y     = N / ( Q * hms_trk_eff * shms_trk_eff * total_LT * mult_trk_eff)              
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

        # calculate the  difference in relative yield (y-axis) and difference in rates (x-axis)
        # this gives the slope:  deltaY / deltaX which gives a correction factor
        deltaY = max(relY) - min(relY)
        deltaY_nom = unumpy.nominal_values(deltaY)  # fractional form (need to multiply by 100 to convert to %)
        deltaY_err = unumpy.std_devs(deltaY)        # fractional form (need to multiply by 100 to convert to %)
        deltaX = max(T2_scl_rate) - min(T2_scl_rate)  # kHz
        slope= deltaY_nom / deltaX
        slope_err = deltaY_err / deltaX
        print('-----------------------------\n')
        print('relY_nom=',relY_nom)
        print('relY_err=',relY_err)
        print('T2_scl_rate=', T2_scl_rate)
        print('deltaY_nom=', deltaY_nom)
        print('deltaY_err=', deltaY_err)
        print('deltaX=',deltaX)
        print('slope=', slope)
        print('slope_err=', slope_err)
        print('-----------------------------\n')

        
        total_LT_nom = unumpy.nominal_values(total_LT )
        total_LT_err = unumpy.std_devs(total_LT )

        ax = axs[0]
        #plt.subplot(1, 2, 1)
        #plot relative yield for files 1, 2
        ax.errorbar(I, relY_nom, relY_err,   marker=imarker[i], markersize=8, color=icolor[i], mec='k', linestyle='dashed')
        #plt.errorbar(I, total_LT_nom, total_LT_err,   marker=imarker[i], markersize=8, color=icolor[i], mec='k', linestyle='dashed')
 
        ax.set_ylabel('Relative Yield', fontsize=18)
        #plt.ylabel('EDTM Live Time', fontsize=18)

        ax.set_xlabel('Average Current [uA]', fontsize=18)
        #ax.set_xticks(fontsize=14)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        #plt.yticks(fontsize=14)
        #plt.grid(True)
        ax.grid(True)
        ax = axs[1]
        #plt.subplot(1, 2, 2)
        #plot relative T2 scalers for files 1, 2

        ax.errorbar(  T2_scl_rate , relY_nom, relY_err,   marker=imarker[i], markersize=8, color=icolor[i], mec='k', linestyle='dashed')
        #plt.errorbar(  T2_scl_rate , total_LT_nom, total_LT_err,   marker=imarker[i], markersize=8, color=icolor[i], mec='k', linestyle='dashed')
     
    
        ax.set_xlabel('T2 Scaler Rate [kHz]', fontsize=18)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.grid(True)
                 
    fig.tight_layout() 
    #plt.legend()
    plt.show()




# calling the relevant functions
compare()
