#!/usr/bin/env python

############################################################
# Integrator of the equations for the orbital evolution of
# 'e' and 'I' in the quadrupolar approximation.
# @bayronportilla - 2016


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

InDeg = 180/np.pi
InRad = np.pi/180
G = 1.0


############################################################
# Initial values
m_0 = 1.0
m_1 = 0.6
m_2 = 0.001

a_1     = 0.1
a_2     = 1.5
e_1_ini = 0.5
e_2     = 0.2
I_ini   = 10*InRad
z_ini   = np.cos(I_ini) 
varpi   = 45*InRad


############################################################
# Derived properties
beta_1 = m_0*m_1 / (m_0 + m_1)
beta_2 = m_2*(m_0 + m_1) / (m_0 + m_1 + m_2)

mu_1 = G*(m_0 + m_1)
mu_2 = G*(m_0 + m_1 + m_2)

gamma_2 = 3.0*G*m_2*beta_1*a_1**2 / ( 4.0*a_2**3*(1.0-e_2**2)**1.5)
gamma_3 = (15.0/64.0)*( G*m_2*beta_1*a_1**3 / (a_2**4*(1.0-e_2**2))**2.5 )*( (m_0 - m_1)/(m_0 + m_1) )



############################################################
# System of first order differential equations. The solution
# will be y = [e(t),I(t)]

def diff(y,t):
    e = y[0]
    z = y[1]
    dedt = 2.5*gamma_2*(1.0-e*e)*e*(1.0-z**2)*np.sin(2*varpi) / np.sqrt( a_1*beta_1**2*mu_1*(1-e**2) )
    dzdt = e*dedt*z / (1-e**2)

    return [dedt, dzdt]



############################################################
# Integration
time = np.linspace(0.0,1000000,1000)
eta_ini = [e_1_ini, z_ini]
sol = odeint(diff,eta_ini,time)



############################################################
# Plotting

fig1 = plt.figure()
ax1 = plt.subplot()
ax1.plot(time , np.arccos(sol[:,1]))
fig1.savefig("i_vs_t.png")

fig2 = plt.figure()
ax2 = plt.subplot()
ax2.plot(time , sol[:,0])
fig2.savefig("e_vs_t.png")










