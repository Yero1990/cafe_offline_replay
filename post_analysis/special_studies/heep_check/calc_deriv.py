import numpy as np


dtr = np.pi/180.

Mp=0.938272

# central kinematic settings
the = 8.3
Ef =  8.52264
Ei = 10.549

# W derivatives 
dW_dEb = Ef / Ei
dW_dEf = -1.* Ei / Ef
dW_dthe = -2.*Ei*Ef*np.sin(the*dtr/2.)*np.cos(the*dtr/2.) / Mp

# measured W difference between data and simc
dW = -0.009   # W_data - W_simc [GeV]
#dEf = 0.001
#dthe = 0.000356  # offset in radians (e- angle) from Pmx fit

# by how much would Eb, Ef or the need to be changed to fully account for dW
dEb = dW / dW_dEb
dEf = dW / dW_dEf
dthe = dW / dW_dthe

#dW = dEf * dW_dEf
#dW = dthe * dW_dthe

#print('dW = %.3f'% dW)

print('dW_dEb = %.3f \n' % (dW_dEb),
      'dW_dEf = %.3f \n' % (dW_dEf),
      'dW_dthe = %.3f \n' % (dW_dthe),
      'dEb = %.4f \n' % (dEb),
      'dEf = %.4f \n' % (dEf),
      'dthe = %.4f \n' % (dthe),
      )
