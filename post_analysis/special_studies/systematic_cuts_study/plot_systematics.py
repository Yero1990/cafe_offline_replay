import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

fname='cafe_systematics_Be9_MF.csv'
df_data = pd.read_csv(fname, comment='#')


# find index of minimum/maximum yield
idx_min = df_data[['real_yield']].idxmin()
idx_max = df_data[['real_yield']].idxmax()

print('(idx, real_yield_min): (%i, %.1f)' % (idx_min, df_data['real_yield'].min() ) )
print('(idx, real_yield_max): (%i, %.1f)' % (idx_max, df_data['real_yield'].max() ) )

'''
#Covert Panda dataframe to numpy
df_numpy=df_data.to_numpy(dtype ='float32')


real_yield = df_numpy[:,1]  #integrated missing momentum counts


plt.hist(real_yield, bins=80, density=True, alpha=0.6, color='b')
mu, std = norm.fit(real_yield) 
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
y = norm.pdf(x, mu, std)
plt.plot(x, y, 'r', linewidth=2)
title = r'Be9 MF: $\mu$={:.2f}, $\sigma$= {:.2f}'.format(mu, std)

plt.ylabel('Frequency')
plt.xlabel('Integrated Missing Momentum')
plt.title(title)

plt.show()
'''
