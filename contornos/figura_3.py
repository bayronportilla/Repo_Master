############################################################
# Figura 3 del paper de correia 2016

import numpy as np
import matplotlib.pyplot as plt
from librerias import Units

unidades_canonicas = Units.units(uM=1.9885e30,uL=0.1*149.6e9)
uM = unidades_canonicas[0]
uL = unidades_canonicas[1]
uT = unidades_canonicas[2]

G = 1.0

m_0 = 1.0*1.98885e30 / uM
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


x = np.arange(-np.pi, np.pi, 0.01)
y = np.arange(0.0, 0.3, 0.01)
X,Y = np.meshgrid(x,y)

U_av = -gamma_2/3.0 * (1.0 + 1.5*e_1**2) + 4.0*gamma_3*(1.0+3.0/4.0*e_1**2)*e_1*Y*np.cos(X) + 0.5*(e_1**2 + X**2)

cs = plt.contour(X,Y,U_av)

"""
plt.figure()
plt.contour(X,Y,U_av)
plt.clabel(cs, inline=1, fontsize=10)
plt.savefig("contorno_1.png")
"""

plt.figure()
CS=plt.contour(X,Y,U_av)
plt.clabel(CS, inline=1, fontsize=10)
plt.savefig("contorno_1.png")