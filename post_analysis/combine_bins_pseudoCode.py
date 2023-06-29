import numpy as np
import pandas as pd
from itertools import takewhile
import matplotlib.pyplot as plt
import sys
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy

df = pd.read_csv('realYield_Pm_bins_Ca40_SRC_17008.csv', comment='#')

run = np.array([3288, 3290, 3290])
real_Yield_per_bin_total = 0

for i in np.arange(len(run)):
    print(i, run[i])

    # read the real yield per bin per run
    real_Yield_per_bin =  unumpy.uarray(df['ycont'], df['ycont_err']) 

    # apply efficiency correction (scale factor) per bin per run
    real_Yield_per_bin_eff = real_Yield_per_bin / (0.99 * 0.98 *0.94 )
    
    # add the eff. corr.  yield bin by bin per run to get total yield per bin over all runs
    real_Yield_per_bin_total = real_Yield_per_bin_total + real_Yield_per_bin
    
    print(real_Yield_per_bin_total)
    print(real_Yield_per_bin_total.sum())

# this array represents the bin-by-bin eff. corrected counts
real_Yield_per_bin_total 
