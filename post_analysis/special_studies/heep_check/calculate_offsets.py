import numpy as np

# this script will solve a system of linear equations
# to determine the optimal angular offsets to align Pmx

#dPmx = E'*cos(theta_e)*dtheta_e - pf*cos(theta_p)dtheta_p
#dPmz = E'*sin(theta_e)*dtheta_e + pf*sin(theta_p)dtheta_p


a_00 = 3. # Ef * np.cos(theta_e)
a_01 = 2 #-pf * np.cos(theta_p)
a_10 = 0 # Ef * np.sin(theta_e)
a_11 = 6 # pf * np.sin(theta_p)


A = np.array([[a_00, a_01], 
              [a_10, a_11]])

A_inv = np.linalg.inv(A)

#dth = [ dtheta_e, dtheta_p ]
#dPm = [ dPmx, dPmz ]


result = A_inv.dot(x)
print('result=',result)

