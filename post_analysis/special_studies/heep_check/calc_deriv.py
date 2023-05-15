import numpy as np
import matplotlib.pyplot as plt

dtr = np.pi/180.

Mp=0.938272

# central kinematic settings
the = 8.295
Ef =  8.55 * 0.9967
Pf = 1.820 * 0.999
Ei = 10.549
Ep = np.sqrt(Mp**2 + Pf**2)

# W derivatives 
dW_dEb = Ef / Ei
dW_dEf = -1.* Ei / Ef
dW_dthe = -2.*Ei*Ef*np.sin(the*dtr/2.)*np.cos(the*dtr/2.) / Mp


# vary Eb, Ef by same amount to keep Em fixed, but to re-align W
dEb = 0.004
dEf = -0.004

dW = dW_dthe * 0.000627
dW_fromEb = dW_dEb * dEb
dW_fromEf = dW_dEf * dEf

print('dW_fromthe = %.4f'%dW)
print('dW_fromEb = %.4f'%dW_fromEb)
print('dW_fromEf = %.4f'%dW_fromEf)


# Emiss derivatives
dEm_dEb = 1
dEm_dEf = -1
dEm_dPf = -Pf/Ep

#plt.plot(Ei, dW_dEb, marker='None', linestyle='-', color='b', label='dW_dEb')
#plt.plot(Ei, dEm_dEb, marker='None', linestyle='-', color='b', label='dEm_dEb')

#plt.plot(Ei/1000., dW_dEf, marker='None', linestyle='-', color='b', label='dW_dEf')

#plt.plot(Ei/1000., dW_dthe, marker='None', linestyle='-', color='r', label='dW_dthe')


#plt.legend()
#plt.show()

#print('--- W derivarives---\n'
#    'dW_dEb = %.4f \n' % (dW_dEb),
#    'dW_dEf = %.4f \n' % (dW_dEf),
#    'dW_dthe = %.4f \n' % (dW_dthe),
#)

#print('--- dW variations ---\n'
#    'dthe = %.4f [rad] -> dW = %.4f ' % (dthe, dW)
#)
      ##      'dEb = %.4f \n' % (dEb),
#      'dEf = %.4f \n' % (dEf),
#      'dthe = %.4f \n' % (dthe),
#      )
