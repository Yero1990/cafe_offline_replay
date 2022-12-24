import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame as pdf
import sys

#user arguments
#target=sys.argv[1]  # LH2, LD2,Be9, B10, B11, C12, Ca40, Ca48, Fe54
#kin=sys.argv[2]     # MF, SRC heep_coin, heep_singles, optics

def plot_data(xkey='', ykey='', ykey_err='', ylo=0, yhi=0, target=[], tcolor=[], kin=[], kmarker=[], calc_ratio=False):

    plt_err_flg = False
    if len(ykey_err)>0:
        plt_err_flg = True
        
    for idx in np.arange(len(target)):
    
        for jdx in np.arange(len(kin)):
            
            filename='summary_files/pass1/cafe_prod_%s_%s_report_summary.csv' % (target[idx], kin[jdx])
            print('Reading filename: %s' % filename)

            df = pd.read_csv(filename, comment='#')

            
            # get any to column
            x_data = df[xkey]
            y_data = df[ykey]

            run = df['run']
            
            # assumes that ykey is the numerator and xkey is denominator (from the user input)
            if(calc_ratio == True):                
                ratio = y_data / x_data

                # will need to put a requirement here (only applies to yield/mC)
                if( (kin[jdx] == "MF") and (target[idx] == "LD2") ):
                    ratio = ratio / 1000.
                if( (kin[jdx] == "MF") and (target[idx] != "LD2") ):
                    ratio = ratio / 100.
                
            if(plt_err_flg == True):
                y_data_err = df[ykey_err]
                
                if(calc_ratio == True):                    
                    ratio_err =  y_data_err / x_data


                    if( (kin[jdx] == "MF") and (target[idx] == "LD2") ):
                        ratio_err = ratio_err / 1000.
                    if( (kin[jdx] == "MF") and (target[idx] != "LD2") ):
                        ratio_err = ratio_err / 100.
                    
            if(ylo==0 and yhi==0):
                print('No y-range provided, will set to default')
                # set y-axis range
                y_lo = min(y_data) - 0.25*min(y_data)
                y_hi = max(y_data) + 0.25*max(y_data)
            else:
                y_lo = ylo
                y_hi = yhi
            plt.ylim(y_lo, y_hi)



            # plot the data
            if(calc_ratio == False):
                
                if(plt_err_flg == False):
                    plt.plot(x_data, y_data, marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
                if(plt_err_flg == True):
                    plt.errorbar(x_data, y_data, y_data_err, marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))

                plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')

            elif(calc_ratio == True):

                if(plt_err_flg == False):
                    plt.plot(run, ratio, marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
                if(plt_err_flg == True):
                    plt.grid(True)
                    plt.errorbar(run, ratio, ratio_err, marker=kmarker[jdx], color=tcolor[idx], mec='k', linestyle='None', label='%s %s'%(target[idx], kin[jdx]))
                    

            plt.legend()

    plt.figure(1)
    plt.show()


#--------------------------
# plot quality-check data
#--------------------------

'''
# total live time
plot_data(xkey='run', ykey='tLT', ykey_err='tLT_err_Bi', ylo=0.96, yhi=1.05, target=['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54'], tcolor=['c', 'm', 'r', 'g', 'b', 'darkorange', 'violet', 'gold'], kin=['MF', 'SRC'], kmarker=['o','v'])

# hms tracking efficiency
plot_data(xkey='run', ykey='hTrkEff', ykey_err='hTrkEff_err', ylo=0.96, yhi=1.05, target=['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54'], tcolor=['c', 'm', 'r', 'g', 'b', 'darkorange', 'violet', 'gold'], kin=['MF', 'SRC'], kmarker=['o','v'])

# shms tracking efficiency
plot_data(xkey='run', ykey='pTrkEff', ykey_err='pTrkEff_err', ylo=0.96, yhi=1.05, target=['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54'], tcolor=['c', 'm', 'r', 'g', 'b', 'darkorange', 'violet', 'gold'], kin=['MF', 'SRC'], kmarker=['o','v'])
'''

#--------------------------
# plot quality-check data
#--------------------------
#plot_data(xkey='charge', ykey='real_Yield', ykey_err='real_Yield_err', ylo=0., yhi=10, target=['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54'], tcolor=['c', 'chocolate', 'r', 'g', 'b', 'darkorange', 'violet', 'gold'], kin=['MF', 'SRC'], kmarker=['o','v'], calc_ratio=True)

plot_data(xkey='charge', ykey='T2_scl_rate', ykey_err='', ylo=0., yhi=10, target=['LD2', 'Be9', 'B10', 'B11', 'C12', 'Ca40', 'Ca48', 'Fe54'], tcolor=['c', 'chocolate', 'r', 'g', 'b', 'darkorange', 'violet', 'gold'], kin=['MF', 'SRC'], kmarker=['o','v'], calc_ratio=True)















#get_data('run', 'Be9', 'SRC')

#plt.plot( get_data('run', 'Be9', 'SRC'),
#              get_data('ctime_sigma', 'Be9', 'SRC'),
#              marker='o', color='r'
#)
#plt.show()

'''
#generate filename from user argument
filename='cafe_prod_%s_%s_report_summary.csv' % (target, kin)

# read .csv file into datraframe
df = pd.read_csv(filename, comment='#')

run = df['run']
charge = df['charge']  # mC

Yield = df['real_Yield']
Yield_err = df['real_Yield_err']

htrk_eff = df['hTrkEff']
htrk_eff_err = df['hTrkEff_err']

ptrk_eff = df['pTrkEff']
ptrk_eff_err = df['pTrkEff_err']

tLT = df['tLT']
tLT_err = df['tLT_err_Bi']

PS1 = df['PS1']
PS2 = df['PS2']
PS3 = df['PS3']
PS5 = df['PS5']

ctime_offset = df['ctime_offset']
ctime_sigma = df['ctime_sigma']

hbeta_mean = df['hbeta_mean']
hbeta_sigma = df['hbeta_sigma']

pbeta_mean = df['pbeta_mean']
pbeta_sigma = df['pbeta_sigma']

pcal_mean = df['pcal_mean']
pcal_sigma = df['pcal_sigma']

#                     rows, col
# hdc_res[:, col] --> array of all rows for each column
hdc_res = np.zeros((len(run),12))
hdc_res_err = np.zeros(12)

pdc_res = np.zeros(12)
pdc_res_err = np.zeros(12)

#define HMS/SHMS plane names in order (beam ----> | | | .  . . )
hdc_pl_name = ["1u1", "1u2", "1x1", "1x2", "1v2", "1v1", "2v1", "2v2", "2x2", "2x1", "2u2", "2u1"]
pdc_pl_name = ["1u1", "1u2", "1x1", "1x2", "1v1", "1v2", "2v2", "2v1", "2x2", "2x1", "2u2", "2u1"]

for idx, i in enumerate(hdc_res):

    print(idx, i)
    # define key names for dirft chamber residuals/errors
    hdc_key_name_1 = 'hdc'+hdc_pl_name[idx]+'_res'
    hdc_key_name_2 = 'hdc'+hdc_pl_name[idx]+'_res_err'

    pdc_key_name_1 = 'pdc'+pdc_pl_name[idx]+'_res'
    pdc_key_name_2 = 'pdc'+pdc_pl_name[idx]+'_res_err'
    print(hdc_key_name_1 )
    # read dc residuals/erros from dafaframe
    hdc_res[idx]     = df[hdc_key_name_1][idx]
    hdc_res_err[idx] = df[hdc_key_name_2][idx]

    print('hdc_res = ', hdc_res[idx])
    
    pdc_res[idx]     = df[pdc_key_name_1][idx]
    pdc_res_err[idx] = df[pdc_key_name_2][idx]

    

    

    
# Make Quality Check Plots
#plt.errorbar(run, htrk_eff, yerr = htrk_eff_err, fmt ='o', color='blue', ecolor='blue')
#plt.show()

# create a new file and wirte header and lines (example)
#myfile = open('ca48_correction.csv', 'w')
#myfile.write('idx, kin, run, Q, Qsum, H_absCntm_MF, C_absCntm_MF, C_absCntm_SRC\n')
# line="%i, %s, %i, %.3f, %.3f, %.3f, %.3f, %.3f\n" % (idx, kin[idx], irun, Qbcm1[idx], Qsum, H_absCntm_MF, C_absCntm_MF, C_absCntm_SRC)
#    myfile.write(line)
'''
