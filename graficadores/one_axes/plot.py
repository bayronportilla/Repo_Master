import numpy as np
import matplotlib.pyplot as plt


RS        = 6.957000e8
Year      = 365.25*86400
Myr       = Year*1.0e6 
Day       = 86400.0
InDeg     = 180.0/np.pi
AU = 149.6e9

fig_file  = "test_01.png"



############################################################
# Fancy style
plt.rc('font', **{'family': 'family', 'family': ['serif']})
plt.rc('text', usetex=True) #uses Latex instead of Tex to compile axes labels


fig=plt.figure()
ax=plt.axes()


############################################################
# Plot on each axes

ax.plot(time_2*uT/Myr,R_A_2*uL/RS)

############################################################
# ticks params


############################################################
# labels

ax.set_xlabel("t (Myr)")
ax.set_ylabel("$R \, (R_{\odot})$")


fig.savefig("radio.png")
