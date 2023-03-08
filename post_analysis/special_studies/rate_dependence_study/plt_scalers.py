import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
df = pd.read_csv('scl_rates.csv', comment='#')

df_h = df[df['tgt']=='h']
df_d = df[df['tgt']=='d']
df_c = df[df['tgt']=='c']

#scaler per charge
T2scl_perQ_h = np.array(df_h['T2'] /  df_h['charge'])
T2scl_perQ_d = np.array(df_d['T2'] /  df_d['charge'])
T2scl_perQ_c = np.array(df_c['T2'] /  df_c['charge'])


T2scl_perQ_h_norm = T2scl_perQ_h / T2scl_perQ_h[0]
T2scl_perQ_d_norm = T2scl_perQ_d / T2scl_perQ_d[0] 
T2scl_perQ_c_norm = T2scl_perQ_c / T2scl_perQ_c[0]


#cumulative charge
Qsum_h = np.array(df_h['charge'].cumsum())
Qsum_d = np.array(df_d['charge'].cumsum())
Qsum_c = np.array(df_c['charge'].cumsum())

Ib_h =  np.array(df_h['current'])
Ib_d =  np.array(df_d['current'])
Ib_c =  np.array(df_c['current'])



fig1, axs1 = plt.subplots(2)
axs1[0].title.set_text('Relative T2 scalers per Charge')
axs1[0].plot(Qsum_h, T2scl_perQ_h_norm, marker='o', color='r', mec='k', label='LH2')
axs1[0].plot(Qsum_d, T2scl_perQ_d_norm, marker='^', color='b', mec='k', label='LD2')
axs1[0].plot(Qsum_c, T2scl_perQ_c_norm, marker='s', color='g', mec='k', label='C12')

axs1[0].set_ylabel(r'(T2_scalers/Q) / (T2_scalers/Q)$_{0}$ ')
axs1[0].set_xlabel('Cumulative Charge [mC]')


axs1[1].title.set_text('Relative T2 scalers per Charge')
axs1[1].plot(Qsum_h, Ib_h, marker='o', color='r', mec='k', label='LH2')
axs1[1].plot(Qsum_d, Ib_d, marker='^', color='b', mec='k', label='LD2')
axs1[1].plot(Qsum_c, Ib_c, marker='s', color='g', mec='k', label='C12')

axs1[1].set_ylabel(r'Avg Current ')
axs1[1].set_xlabel('Cumulative Charge [mC]')

plt.tight_layout()

plt.show()
