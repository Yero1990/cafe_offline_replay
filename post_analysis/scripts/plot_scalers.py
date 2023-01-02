import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame as pdf
import sys

# heavy
df_ca40_src = pd.read_csv('cafe_prod_Ca40_SRC_report_summary.csv', comment='#')
df_ca48_src = pd.read_csv('cafe_prod_Ca48_SRC_report_summary.csv', comment='#')
df_fe54_src = pd.read_csv('cafe_prod_Fe54_SRC_report_summary.csv', comment='#')

# light
df_be9_src = pd.read_csv('cafe_prod_Be9_SRC_report_summary.csv', comment='#')
df_b10_src = pd.read_csv('cafe_prod_B10_SRC_report_summary.csv', comment='#')
df_b11_src = pd.read_csv('cafe_prod_B11_SRC_report_summary.csv', comment='#')
df_c12_src = pd.read_csv('cafe_prod_C12_SRC_report_summary.csv', comment='#')



I40 = df_ca40_src['avg_current']
I48 = df_ca48_src['avg_current']
I54 = df_fe54_src['avg_current']

I9  = df_be9_src['avg_current']
I10 = df_b10_src['avg_current']
I11 = df_b11_src['avg_current']
I12 = df_c12_src['avg_current']


Q40 = df_ca40_src['charge']
Q48 = df_ca48_src['charge']
Q54 = df_fe54_src['charge']

Q9  = df_be9_src['charge']
Q10 = df_b10_src['charge']
Q11 = df_b11_src['charge']
Q12 = df_c12_src['charge']

Q9csum  = Q9.cumsum()
Q10csum = Q10.cumsum() 
Q11csum = Q11.cumsum() 
Q12csum = Q12.cumsum() 

Q40csum =  Q40.cumsum()
Q48csum =  Q48.cumsum()
Q54csum =  Q54.cumsum()



T2_scl_ca40 = df_ca40_src['T2_scl_rate']*1000.*df_ca40_src['beam_time']
T2_scl_ca48 = df_ca48_src['T2_scl_rate']*1000.*df_ca48_src['beam_time']
T2_scl_fe54 = df_fe54_src['T2_scl_rate']*1000.*df_fe54_src['beam_time']

T2_scl_be9 = df_be9_src['T2_scl_rate']*1000.*df_be9_src['beam_time']
T2_scl_b10 = df_b10_src['T2_scl_rate']*1000.*df_b10_src['beam_time']
T2_scl_b11 = df_b11_src['T2_scl_rate']*1000.*df_b11_src['beam_time']
T2_scl_c12 = df_c12_src['T2_scl_rate']*1000.*df_c12_src['beam_time']


T2_scl_per_Q_ca48 = T2_scl_ca48 / Q48
T2_scl_per_Q_ca40 = T2_scl_ca40 / Q40
T2_scl_per_Q_fe54 = T2_scl_fe54 / Q54

T2_scl_per_Q_be9 = T2_scl_be9 / Q9
T2_scl_per_Q_b10 = T2_scl_b10 / Q10
T2_scl_per_Q_b11 = T2_scl_b11 / Q11
T2_scl_per_Q_c12 = T2_scl_c12 / Q12


#heavy
fig1, axs1 = plt.subplots(2)
axs1[0].title.set_text('Relative T2 scalers per Charge')
axs1[0].plot(Q48csum, T2_scl_per_Q_ca48/T2_scl_per_Q_ca48[0], marker='^', color='r', mec='k', label='Ca48 SRC')
axs1[0].plot(Q40csum, T2_scl_per_Q_ca40/T2_scl_per_Q_ca40[0], marker='o', color='b', mec='k', label='Ca40 SRC')
axs1[0].plot(Q54csum, T2_scl_per_Q_fe54/T2_scl_per_Q_fe54[0], marker='s', color='g', mec='k', label='Fe54 SRC')

axs1[0].set_ylabel(r'(T2_scalers/Q) / (T2_scalers/Q)$_{0}$ ')
axs1[0].set_xlabel('Cumulative Charge [mC]')
plt.legend()

axs1[1].title.set_text('Avg. Beam Current vs. Charge')
axs1[1].plot(Q48csum,  I48, marker='^', color='r', mec='k', label=r'Ca48 SRC')
axs1[1].plot(Q40csum,  I40, marker='o', color='b', mec='k', label=r'Ca40 SRC')
axs1[1].plot(Q54csum,  I54, marker='s', color='g', mec='k', label=r'Fe54 SRC')

axs1[1].set_ylabel('Avg Beam Current')
axs1[1].set_xlabel('Cumulative Charge [mC]')
plt.legend()

fig1.tight_layout()
plt.show()


#light

fig2, axs2 = plt.subplots(2)
axs2[0].title.set_text('Relative T2 scalers per Charge')
axs2[0].plot(Q9csum , T2_scl_per_Q_be9/T2_scl_per_Q_be9[0], marker='^', color='r', mec='k', label='Be9 SRC')
axs2[0].plot(Q10csum, T2_scl_per_Q_b10/T2_scl_per_Q_b10[0], marker='o', color='b', mec='k', label='B10 SRC')
axs2[0].plot(Q11csum, T2_scl_per_Q_b11/T2_scl_per_Q_b11[0], marker='s', color='g', mec='k', label='B11 SRC')
axs2[0].plot(Q12csum, T2_scl_per_Q_c12/T2_scl_per_Q_c12[0], marker='D', color='orange', mec='k', label='C12 SRC')

axs2[0].set_ylabel(r'(T2_scalers/Q) / (T2_scalers/Q)$_{0}$ ')
axs2[0].set_xlabel('Cumulative Charge [mC]')
plt.legend()

axs2[1].title.set_text('Avg. Beam Current vs. Charge')
axs2[1].plot(Q9csum , I9 , marker='^', color='r', mec='k', label=r'Be9 SRC')
axs2[1].plot(Q10csum, I10, marker='o', color='b', mec='k', label=r'B10 SRC')
axs2[1].plot(Q11csum, I11, marker='s', color='g', mec='k', label=r'B11 SRC')
axs2[1].plot(Q12csum, I12, marker='D', color='orange', mec='k', label=r'C12 SRC')

axs2[1].set_ylabel('Avg Beam Current')
axs2[1].set_xlabel('Cumulative Charge [mC]')

plt.legend()

fig2.tight_layout()

plt.show()
