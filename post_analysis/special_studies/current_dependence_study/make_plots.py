import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy


def get_data(target='', kin='', header=''):

    # read any data from the .csv file
    
    # define summary file
    summary_file_path = '../../summary_files/pass1_with_tcoin/tight_current_cut/cafe_prod_%s_%s_report_summary.csv' % (target, kin)

    # read .csv file
    df = pd.read_csv(summary_file_path, comment='#')

    #print(df)
    data_array = np.array(df[header])

    return data_array
    

def make_plots():

    print('Calling make_plots()')

    I_be9_mf             = get_data('Be9', 'MF', 'avg_current')
    Q_be9_mf             = get_data('Be9', 'MF', 'charge')
    N_be9_mf             = unumpy.uarray(get_data('Be9', 'MF', 'real_Yield'),  get_data('Be9', 'MF', 'real_Yield_err') )
    hms_trk_eff_be9_mf   = unumpy.uarray(get_data('Be9', 'MF', 'hTrkEff'),     get_data('Be9', 'MF', 'hTrkEff_err') )
    shms_trk_eff_be9_mf  = unumpy.uarray(get_data('Be9', 'MF', 'pTrkEff'),     get_data('Be9', 'MF', 'pTrkEff_err') )
    total_LT_be9_mf      = unumpy.uarray(get_data('Be9', 'MF', 'tLT'),         get_data('Be9', 'MF', 'tLT_err_Bi') )
    mult_trk_eff_be9_mf  = np.array(get_data('Be9', 'MF', 'multi_track_eff') )
    Y_be9_mf = N_be9_mf / ( Q_be9_mf * hms_trk_eff_be9_mf * shms_trk_eff_be9_mf * total_LT_be9_mf * mult_trk_eff_be9_mf)              
    relY_be9_mf  = Y_be9_mf / Y_be9_mf[0] 

    I_b10_mf             = get_data('B10', 'MF', 'avg_current')
    Q_b10_mf             = get_data('B10', 'MF', 'charge')
    N_b10_mf             = unumpy.uarray(get_data('B10', 'MF', 'real_Yield'),  get_data('B10', 'MF', 'real_Yield_err') )
    hms_trk_eff_b10_mf   = unumpy.uarray(get_data('B10', 'MF', 'hTrkEff'),     get_data('B10', 'MF', 'hTrkEff_err') )
    shms_trk_eff_b10_mf  = unumpy.uarray(get_data('B10', 'MF', 'pTrkEff'),     get_data('B10', 'MF', 'pTrkEff_err') )
    total_LT_b10_mf      = unumpy.uarray(get_data('B10', 'MF', 'tLT'),         get_data('B10', 'MF', 'tLT_err_Bi') )
    mult_trk_eff_b10_mf  = np.array(get_data('B10', 'MF', 'multi_track_eff') )
    Y_b10_mf = N_b10_mf / ( Q_b10_mf * hms_trk_eff_b10_mf * shms_trk_eff_b10_mf * total_LT_b10_mf * mult_trk_eff_b10_mf)              
    relY_b10_mf  = Y_b10_mf / Y_b10_mf[0]

    I_b11_mf             = get_data('B11', 'MF', 'avg_current')
    Q_b11_mf             = get_data('B11', 'MF', 'charge')
    N_b11_mf             = unumpy.uarray(get_data('B11', 'MF', 'real_Yield'),  get_data('B11', 'MF', 'real_Yield_err') )
    hms_trk_eff_b11_mf   = unumpy.uarray(get_data('B11', 'MF', 'hTrkEff'),     get_data('B11', 'MF', 'hTrkEff_err') )
    shms_trk_eff_b11_mf  = unumpy.uarray(get_data('B11', 'MF', 'pTrkEff'),     get_data('B11', 'MF', 'pTrkEff_err') )
    total_LT_b11_mf      = unumpy.uarray(get_data('B11', 'MF', 'tLT'),         get_data('B11', 'MF', 'tLT_err_Bi') )
    mult_trk_eff_b11_mf  = np.array(get_data('B11', 'MF', 'multi_track_eff') )
    Y_b11_mf = N_b11_mf / ( Q_b11_mf * hms_trk_eff_b11_mf * shms_trk_eff_b11_mf * total_LT_b11_mf * mult_trk_eff_b11_mf)              
    relY_b11_mf  = Y_b11_mf / Y_b11_mf[1]

    I_c12_mf             = get_data('C12', 'MF', 'avg_current')
    Q_c12_mf             = get_data('C12', 'MF', 'charge')
    N_c12_mf             = unumpy.uarray(get_data('C12', 'MF', 'real_Yield'),  get_data('C12', 'MF', 'real_Yield_err') )
    hms_trk_eff_c12_mf   = unumpy.uarray(get_data('C12', 'MF', 'hTrkEff'),     get_data('C12', 'MF', 'hTrkEff_err') )
    shms_trk_eff_c12_mf  = unumpy.uarray(get_data('C12', 'MF', 'pTrkEff'),     get_data('C12', 'MF', 'pTrkEff_err') )
    total_LT_c12_mf      = unumpy.uarray(get_data('C12', 'MF', 'tLT'),         get_data('C12', 'MF', 'tLT_err_Bi') )
    mult_trk_eff_c12_mf  = np.array(get_data('C12', 'MF', 'multi_track_eff') )
    Y_c12_mf = N_c12_mf / ( Q_c12_mf * hms_trk_eff_c12_mf * shms_trk_eff_c12_mf * total_LT_c12_mf * mult_trk_eff_c12_mf)              
    relY_c12_mf  = Y_c12_mf[-4:] / (Y_c12_mf[-4:])[1]
    
    # -------------
    # heavy nuclei
    #--------------
    I_ca40_mf             = get_data('Ca40', 'MF', 'avg_current')
    Q_ca40_mf             = get_data('Ca40', 'MF', 'charge')
    N_ca40_mf             = unumpy.uarray(get_data('Ca40', 'MF', 'real_Yield'),  get_data('Ca40', 'MF', 'real_Yield_err') )
    hms_trk_eff_ca40_mf   = unumpy.uarray(get_data('Ca40', 'MF', 'hTrkEff'),     get_data('Ca40', 'MF', 'hTrkEff_err') )
    shms_trk_eff_ca40_mf  = unumpy.uarray(get_data('Ca40', 'MF', 'pTrkEff'),     get_data('Ca40', 'MF', 'pTrkEff_err') )
    total_LT_ca40_mf      = unumpy.uarray(get_data('Ca40', 'MF', 'tLT'),         get_data('Ca40', 'MF', 'tLT_err_Bi') )
    mult_trk_eff_ca40_mf  = np.array(get_data('Ca40', 'MF', 'multi_track_eff') )
    Y_ca40_mf = N_ca40_mf / ( Q_ca40_mf * hms_trk_eff_ca40_mf * shms_trk_eff_ca40_mf * total_LT_ca40_mf * mult_trk_eff_ca40_mf)              
    relY_ca40_mf  = Y_ca40_mf / Y_ca40_mf[0]  

    I_ca48_mf             = get_data('Ca48', 'MF', 'avg_current')
    Q_ca48_mf             = get_data('Ca48', 'MF', 'charge')
    N_ca48_mf             = unumpy.uarray(get_data('Ca48', 'MF', 'real_Yield'),  get_data('Ca48', 'MF', 'real_Yield_err') )
    hms_trk_eff_ca48_mf   = unumpy.uarray(get_data('Ca48', 'MF', 'hTrkEff'),     get_data('Ca48', 'MF', 'hTrkEff_err') )
    shms_trk_eff_ca48_mf  = unumpy.uarray(get_data('Ca48', 'MF', 'pTrkEff'),     get_data('Ca48', 'MF', 'pTrkEff_err') )
    total_LT_ca48_mf      = unumpy.uarray(get_data('Ca48', 'MF', 'tLT'),         get_data('Ca48', 'MF', 'tLT_err_Bi') )
    mult_trk_eff_ca48_mf  = np.array(get_data('Ca48', 'MF', 'multi_track_eff') )
    Y_ca48_mf = N_ca48_mf / ( Q_ca48_mf * hms_trk_eff_ca48_mf * shms_trk_eff_ca48_mf * total_LT_ca48_mf * mult_trk_eff_ca48_mf)              
    relY_ca48_mf  = Y_ca48_mf[-3:] / (Y_ca48_mf[-3:])[0]  

    I_fe54_mf             = get_data('Fe54', 'MF', 'avg_current')
    Q_fe54_mf             = get_data('Fe54', 'MF', 'charge')
    N_fe54_mf             = unumpy.uarray(get_data('Fe54', 'MF', 'real_Yield'),  get_data('Fe54', 'MF', 'real_Yield_err') )
    hms_trk_eff_fe54_mf   = unumpy.uarray(get_data('Fe54', 'MF', 'hTrkEff'),     get_data('Fe54', 'MF', 'hTrkEff_err') )
    shms_trk_eff_fe54_mf  = unumpy.uarray(get_data('Fe54', 'MF', 'pTrkEff'),     get_data('Fe54', 'MF', 'pTrkEff_err') )
    total_LT_fe54_mf      = unumpy.uarray(get_data('Fe54', 'MF', 'tLT'),         get_data('Fe54', 'MF', 'tLT_err_Bi') )
    mult_trk_eff_fe54_mf  = np.array(get_data('Fe54', 'MF', 'multi_track_eff') )
    Y_fe54_mf = N_fe54_mf / ( Q_fe54_mf * hms_trk_eff_fe54_mf * shms_trk_eff_fe54_mf * total_LT_fe54_mf * mult_trk_eff_fe54_mf)              
    relY_fe54_mf  = Y_fe54_mf / Y_fe54_mf[0]  

    plt.subplot(1, 2, 1)
    plt.errorbar(I_be9_mf,       unumpy.nominal_values(Y_be9_mf)/Y_be9_mf[0].n,             unumpy.std_devs(Y_be9_mf)/Y_be9_mf[0].n,             marker='o',markersize=8, color='r', mec='k', linestyle='dashed', label='Be9 MF (yield)')
    plt.errorbar(I_b10_mf,       unumpy.nominal_values(Y_b10_mf)/Y_b10_mf[0].n,             unumpy.std_devs(Y_b10_mf)/Y_b10_mf[0].n,             marker='o',markersize=8, color='b', mec='k', linestyle='dashed', label='B10 MF (yield)')
    plt.errorbar(I_b11_mf,       unumpy.nominal_values(Y_b11_mf)/Y_b11_mf[1].n,             unumpy.std_devs(Y_b11_mf)/Y_b11_mf[1].n,             marker='o',markersize=8, color='g', mec='k', linestyle='dashed', label='B11 MF (yield)')
    plt.errorbar(I_c12_mf[-4:],       unumpy.nominal_values(Y_c12_mf[-4:])/(Y_c12_mf[-4:])[1].n,             unumpy.std_devs(Y_c12_mf[-4:])/(Y_c12_mf[-4:])[1].n,             marker='o',markersize=8, color='darkorange', mec='k', linestyle='dashed', label='C12 MF (yield)')

    print('I_c12_mf -->',I_c12_mf)
    print('Y_c12_mf --->', unumpy.nominal_values(Y_c12_mf))
    
    plt.legend()

    
    plt.subplot(1, 2, 2)
    plt.errorbar(I_ca48_mf[-3:],  unumpy.nominal_values(Y_ca48_mf[-3:])/(Y_ca48_mf[-3:])[0].n, unumpy.std_devs(Y_ca48_mf[-3:])/(Y_ca48_mf[-3:])[0].n, marker='o',markersize=8, color='b', mec='k', linestyle='dashed', label='Ca48 MF (yield)')
    plt.errorbar(I_ca40_mf,       unumpy.nominal_values(Y_ca40_mf)/Y_ca40_mf[0].n,             unumpy.std_devs(Y_ca40_mf)/Y_ca40_mf[0].n,             marker='o',markersize=8, color='g', mec='k', linestyle='dashed', label='Ca40 MF (yield)')
    plt.errorbar(I_fe54_mf,       unumpy.nominal_values(Y_fe54_mf)/Y_fe54_mf[1].n,             unumpy.std_devs(Y_fe54_mf)/Y_fe54_mf[1].n,             marker='o',markersize=8, color='m', mec='k', linestyle='dashed', label='Fe54 MF (yield)')

    plt.legend()
    plt.show()
    #print(relY_ca48_mf)
    #print(Y[-3:]/(Y[-3:])[0])
    #print(unumpy.nominal_values(Y)[-3:])
    #print(unumpy.std_devs(Y)[-3:])

    
make_plots()
