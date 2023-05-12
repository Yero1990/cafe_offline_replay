import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


# description: this script reads the data yield
# over the systematically-varied cuts for each 
# (target, kinematics) as well as the ratios of
# SRC/MF for target A to a reference nucleus, etc.
# to study the cut sensitivity on ratios as well
# NOTE: When studying cut sensitivity on ratios,
# the ratio of yields should be taken on an entry-by-entry basis
# i.e., ratio of numerical arrays

fname='cafe_systematics_Be9_SRC.csv'
df_data = pd.read_csv(fname, comment='#')


# find index of minimum/maximum yield
idx_min = df_data[['real_yield']].idxmin()
idx_max = df_data[['real_yield']].idxmax()

real_yield = df_data['real_yield']  #integrated missing momentum counts


plt.hist(real_yield, bins=80, density=True, alpha=0.6, color='b')
mu, std = norm.fit(real_yield) 
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
y = norm.pdf(x, mu, std)
plt.plot(x, y, 'r', linewidth=2)
title = r'Be9 SRC: $\mu$={:.2f}, $\sigma$= {:.2f}'.format(mu, std)

plt.ylabel('Frequency')
plt.xlabel('Integrated Missing Momentum')
plt.title(title)

plt.show()

