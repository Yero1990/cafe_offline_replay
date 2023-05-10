import numpy as np


dtr = np.pi/180.

Mp=0.938272

# central kinematic settings
the = 8.295
Ef =  8.55
Ei = 10.549

# W derivatives 
dW_dEb = Ef / Ei
dW_dEf = -1.* Ei / Ef
dW_dthe = -2.*Ei*Ef*np.sin(the*dtr/2.)*np.cos(the*dtr/2.) / Mp

# variations in kinmeatics dE, dEf dthe (to find by how much would W shift)
dthe = 0.001  # rad
dW = dW_dthe * dthe

# by how much would Eb, Ef or the need to be changed to fully account for dW
#dEb = dW / dW_dEb
#dEf = dW / dW_dEf
#dthe = dW / dW_dthe



print('--- W derivarives---\n'
    'dW_dEb = %.4f \n' % (dW_dEb),
    'dW_dEf = %.4f \n' % (dW_dEf),
    'dW_dthe = %.4f \n' % (dW_dthe),
)

print('--- dW variations ---\n'
    'dthe = %.4f [rad] -> dW = %.4f ' % (dthe, dW)
)
      ##      'dEb = %.4f \n' % (dEb),
#      'dEf = %.4f \n' % (dEf),
#      'dthe = %.4f \n' % (dthe),
#      )