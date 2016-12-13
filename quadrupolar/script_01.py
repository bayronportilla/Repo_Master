#!/usr/bin/env python

############################################################
# Integrator of the equations for the orbital evolution of
# 'e' and 'I' in the quadrupolar approximation.
# @bayronportilla - 2016


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from Units import units

InDeg = 180/np.pi
InRad = np.pi/180
G = 1.0


############################################################
# Define canonical units


Canonical_Units = units(uM=1.989e30,uL=149.6e9)

uM = Canonical_Units[0]
uL = Canonical_Units[1]
uT = Canonical_Units[2]


############################################################
# Constant values
m_0 = 1.0
m_1 = 0.2
m_2 = 0.001


############################################################
# Initial values
a_1_ini   = 0.1
a_2_ini   = 1.5
e_1_ini   = 0.5
e_2_ini   = 0.2
I_ini     = 10*InRad
z_ini     = np.cos(I_ini) 
varpi_ini = 45*InRad



############################################################
# Derived properties
beta_1 = m_0*m_1 / (m_0 + m_1)
beta_2 = m_2*(m_0 + m_1) / (m_0 + m_1 + m_2)

mu_1 = G*(m_0 + m_1)
mu_2 = G*(m_0 + m_1 + m_2)

gamma_2 = 3*G*m_2*beta_1*a_1_ini**2 / ( 4.0*a_2_ini**3*(1.0-e_2_ini**2)**1.5)
gamma_3 = (15.0/64.0)*( G*m_2*beta_1*a_1_ini**3 / (a_2_ini**4*(1.0-e_2_ini**2))**2.5 )*( (m_0 - m_1)/(m_0 + m_1) )



############################################################
# System of first order differential equations. The solution
# will be y = [e(t),I(t)]

def diff(y,t):
    e = y[0]
    z = y[1]
    dedt = 2.5*gamma_2*(1.0-e**2)*e*(1.0-z**2)*np.sin(2*varpi_ini) / np.sqrt( a_1_ini*beta_1**2*mu_1*(1-e**2) )
    dzdt = e * dedt * z / (1-e**2)

    return [dedt, dzdt]



############################################################
# Integration
t_ini = ( 0 )      *86400*365.25/uT   
t_end = ( 100000 ) *86400*365.25/uT
time = np.linspace(t_ini, t_end, 1000)

eta_ini = [e_1_ini, z_ini]

sol = odeint(diff,eta_ini,time)



############################################################
# Lidov-Kozai cycles

h_value = np.zeros(len(time))



for i in range(0,len(h_value)):
    h_value[i] = np.sqrt(1-sol[i][0]**2)*sol[i][1]




############################################################
# Plotting

fig1 = plt.figure()
ax1 = plt.subplot()
ax1.grid()
ax1.plot(time*uT/(86400*365.25)/1000 , np.arccos(sol[:,1])*InDeg)
ax1.set_xlabel("time [kyr]")
ax1.set_ylabel("$I$ [degrees]")
fig1.savefig("i_vs_t.png")

fig2 = plt.figure()
ax2 = plt.subplot()
ax2.grid()
ax2.plot(time*uT/(86400*365.25)/1000 , sol[:,0])
ax2.set_xlabel("time [kyr]")
ax2.set_ylabel("$e_1$")
fig2.savefig("e_vs_t.png")

fig3 = plt.figure()
ax3 = plt.subplot()
ax3.plot(time*uT/(86400*365.25)/1000 , h_value)
ax3.set_xlabel("time [kyr]")
ax3.set_ylabel("$\sqrt{1-e_1^2}\cos I$")
ax3.set_ylim(0.85286820 , 0.85286880)
ax3.grid()
fig3.savefig("h_value.png")










