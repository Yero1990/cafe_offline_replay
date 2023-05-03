import numpy as np


dtr = np.pi/180.

Mp=0.938272

# central kinematic settings
the = 8.3
Ef =  8.52264
Ei = 10.553

# W derivatives 
dW_dEb = Ef / Ei
dW_dEf = -1.* Ei / Ef
dW_dthe = -2.*Ei*Ef*np.sin(the*dtr/2.)*np.cos(the*dtr/2.) / Mp

# measured W difference between data and simc
# W_data - W_simc [GeV]
#dEf = 0.001
#dthe = 0.000356  # offset in radians (e- angle) from Pmx fit

# by how much would Eb, Ef or the need to be changed to fully account for dW
#dEb = dW / dW_dEb
#dEf = dW / dW_dEf
#dthe = dW / dW_dthe


dW_1 = -0.001
dW_2 = 0.001

dEb = dW_1 / dW_dEb
dEf = dW_2 / dW_dEf

delta = dEb - dEf

dthe = 0.00052982939
dW = dW_dthe * dthe
print('dW_dEb = %.3f \n' % (dW_dEb),
      'dW_dEf = %.3f \n' % (dW_dEf),
      'dW_dthe = %.3f \n' % (dW_dthe),
      'dth_e = %.5f [rad]-> dW = %.5f [GeV]\n' % (dthe, dW),
      'dEb = %.4f-> dW = %.4f [GeV] \n' % (dEb, dW_1),
      'dEf = %.4f-> dW = %.4f [GeV] \n' % (dEf, dW_2),
      'dEb - dEf = %.4f \n' % (delta)
      )
      ##      'dEb = %.4f \n' % (dEb),
#      'dEf = %.4f \n' % (dEf),
#      'dthe = %.4f \n' % (dthe),
#      )
