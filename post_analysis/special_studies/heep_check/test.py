import numpy as np
import matplotlib.pyplot as plt
th_e = np.linspace(0,0.05, 1000) # radians
Eb = 10.549
kfc = 8.55
kf = 9.45

dthp_dthe_c = kfc*(kfc-Eb*np.cos(th_e)) / ( kfc**2 - 2.*Eb*np.cos(th_e) + Eb**2 )
dthp_dthe = kf*(kf-Eb*np.cos(th_e)) / ( kf**2 - 2.*Eb*np.cos(th_e) + Eb**2 )

#plt.plot(th_e, dthp_dthe, linestyle='-', color='k')
#plt.plot(th_e, dthp_dthe_c, linestyle='-', color='b')

# dthp / dthe = f(Eb, kf, the)

# dthe = dthp / f(Eb,kf,the)

th_e = 8.3 * np.pi/180.
dthe = -0.000699855 / ( kf*(kf-Eb*np.cos(th_e)) / ( kf**2 - 2.*Eb*np.cos(th_e) + Eb**2 )  )
print('dthe=',dthe)

plt.show()

