############################################################
# Figura 3 del paper de correia 2016

import numpy as np
import matplotlib.pyplot as plt
from librerias import Units
from numpy.linalg import norm
from scipy.integrate import odeint

unidades_canonicas = Units.units(uM=1.9885e30,uL=0.1*149.6e9)
uM = unidades_canonicas[0]
uL = unidades_canonicas[1]
uT = unidades_canonicas[2]

G = 1.0

m_0 = 10.0*1.98885e30 / uM
m_1 = 0.2*1.98885e30 / uM
m_2 = 0.001*1.98885e30 / uM

a_1 = 0.1*149.6e9 / uL 
a_2 = 1.5*149.6e9 / uL

e_1 = 0.5
e_2 = 0.2

beta_1 = m_0*m_1 / (m_0 + m_1)
beta_2 = m_2*(m_0 + m_1) / (m_0 + m_1 + m_2)
mu_1 = G*(m_0 + m_1)
mu_2 = G*(m_0 + m_1 + m_2)

gamma_2 = 3.0*G*m_2*beta_1*a_1**2 / ( 4.0*a_2**3*(1.0-e_2**2)**1.5)
gamma_3 = (15.0/64.0)*( G*m_2*beta_1*a_1**3 / (a_2**4*(1.0-e_2**2))**2.5 )*( (m_0 - m_1)/(m_0 + m_1) )


def eta(y,t):
    rho_1 = y[0:3]
    rho_2 = y[3:6]
    v_rho_1 = y[6:9]
    v_rho_2 = y[9:12]

    return [ v_rho_1[0] , v_rho_1[1] , v_rho_1[2] ,\
             v_rho_2[0] , v_rho_2[1] , v_rho_2[2] ,\
             -G*(m_0+m_1)*(rho_1[0]/norm(rho_1)**3 + G*m_2*( rho_2[0] - m_0/(m_0+m_1)*rho_1[0] )/norm (rho_2 - m_0/(m_0+m_1)*rho_1)**3 \
                           - (rho_2[0] + m_2/(m_0+m_1)*rho_1[0]) / norm(rho_2 + m_2/(m_0+m_1)*rho_1)**3 ) , \
             -G*(m_0+m_1)*(rho_1[1]/norm(rho_1)**3 + G*m_2*( rho_2[1] - m_0/(m_0+m_1)*rho_1[1] )/norm (rho_2 - m_0/(m_0+m_1)*rho_1)**3 \
                           - (rho_2[1] + m_2/(m_0+m_1)*rho_1[1]) / norm(rho_2 + m_2/(m_0+m_1)*rho_1)**3 ) , \
             -G*(m_0+m_1)*(rho_1[2]/norm(rho_1)**3 + G*m_2*( rho_2[2] - m_0/(m_0+m_1)*rho_1[2] )/norm (rho_2 - m_0/(m_0+m_1)*rho_1)**3 \
                           - (rho_2[2] + m_2/(m_0+m_1)*rho_1[2]) / norm(rho_2 + m_2/(m_0+m_1)*rho_1)**3 ) , \
             -G*(m_0+m_1+m_2)/(m_0+m_1)*( m_1*(rho_2[0] - m_0/(m_0+m_1)*rho_1[0])/norm(rho_2 - m_0/(m_0+m_1)*rho_1)**3 \
                                          +m_0*(rho_2[0] + m_1/(m_0+m_1)*rho_1[0])/norm(rho_2 + m_1/(m_0+m_1)*rho_1)**3 ) , \
             -G*(m_0+m_1+m_2)/(m_0+m_1)*( m_1*(rho_2[1] - m_0/(m_0+m_1)*rho_1[1])/norm(rho_2 - m_0/(m_0+m_1)*rho_1)**3 \
                                          +m_0*(rho_2[1] + m_1/(m_0+m_1)*rho_1[1])/norm(rho_2 + m_1/(m_0+m_1)*rho_1)**3 ) , \
             -G*(m_0+m_1+m_2)/(m_0+m_1)*( m_1*(rho_2[2] - m_0/(m_0+m_1)*rho_1[2])/norm(rho_2 - m_0/(m_0+m_1)*rho_1)**3 \
                                          +m_0*(rho_2[2] + m_1/(m_0+m_1)*rho_1[2])/norm(rho_2 + m_1/(m_0+m_1)*rho_1)**3 ) ]


eta_0 = [ 1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.001, 0.0 ]
time = np.linspace(0,100,1000)
sol = odeint(eta,eta_0,time,mxstep=10000)

file=open("jacobi.dat","w")
for i in range (0,sol.shape[0]):
    file.write("%e   %e  %e   %e  %e   %e  \n" % (sol[i][0],sol[i][1],sol[i][2],sol[i][3],sol[i][4],sol[i][5] ))
print sol

             
             
