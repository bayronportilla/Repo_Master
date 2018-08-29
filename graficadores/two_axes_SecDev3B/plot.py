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
# Data from simulations

sim_name_1  = "test"
sim_name_2  = "test_SE"

data_file_1 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b.dat"
data_file_2 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b_SE.dat"

info_file_1 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b.log"
info_file_2 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b_SE.log"


data_1      = np.loadtxt(data_file_1)
data_2      = np.loadtxt(data_file_2)

info_1      = np.genfromtxt(info_file_1,dtype=None)
#info_2      = np.genfromtxt(info_file_2,dtype=None)

uT        = float(info_1[34][2])
uL        = float(info_1[33][2])

time_1     = data_1[:,0]
a_in_1     = data_1[:,1:2]
a_out_1    = data_1[:,2:3]
e_in_1     = data_1[:,3:4]
e_out_1    = data_1[:,4:5]
I_in_1     = data_1[:,5:6]
I_out_1    = data_1[:,6:7]
W_in_1     = data_1[:,7:8]
W_out_1    = data_1[:,8:9]
w_in_1     = data_1[:,9:10]
w_out_1    = data_1[:,10:11]
Om_Ax_in_1 = data_1[:,11:12]
Om_Ay_in_1 = data_1[:,12:13]
Om_Az_in_1 = data_1[:,13:14]
Om_Bx_in_1 = data_1[:,14:15]
Om_By_in_1 = data_1[:,15:16]
Om_Bz_in_1 = data_1[:,16:17]
R_A_1      = data_1[:,19:20]
R_B_1      = data_1[:,20:21]

time_2     = data_2[:,0]
a_in_2     = data_2[:,1:2]
a_out_2    = data_2[:,2:3]
e_in_2     = data_2[:,3:4]
e_out_2    = data_2[:,4:5]
I_in_2     = data_2[:,5:6]
I_out_2    = data_2[:,6:7]
W_in_2     = data_2[:,7:8]
W_out_2    = data_2[:,8:9]
w_in_2     = data_2[:,9:10]
w_out_2    = data_2[:,10:11]
Om_Ax_in_2 = data_2[:,11:12]
Om_Ay_in_2 = data_2[:,12:13]
Om_Az_in_2 = data_2[:,13:14]
Om_Bx_in_2 = data_2[:,14:15]
Om_By_in_2 = data_2[:,15:16]
Om_Bz_in_2 = data_2[:,16:17]
R_A_2      = data_2[:,19:20]
R_B_2      = data_2[:,20:21]


# axes rect in relative 0,1 coords left, bottom, width, height.  Turn
# off xtick labels on all but the lower plot

############################################################
# Fancy style
plt.rc('font', **{'family': 'family', 'family': ['serif']})
plt.rc('text', usetex=True) #uses Latex instead of Tex to compile axes labels

fig,((ax_1,ax_2),(ax_3,ax_4)) = plt.subplots(2, 2,figsize=(10,8),sharex=True)
plt.subplots_adjust(hspace=0.0)

############################################################
# Plot on each axes

ax_1.plot(time_1*uT/Myr, a_in_1*uL/AU)
ax_2.plot(time_2*uT/Myr, a_in_2*uL/AU)
ax_3.plot(time_1*uT/Myr, e_in_1)
ax_4.plot(time_2*uT/Myr, e_in_2)

############################################################
# ticks params


############################################################
# labels

ax_3.set_xlabel("t [Myr]")
ax_3.set_ylabel("$e$")
ax_4.set_xlabel("t [Myr]")
ax_4.set_ylabel("$e_1$")
ax_1.set_ylabel("$a_1$ [AU]")
ax_2.set_ylabel("$a_1$ [AU]")


#xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
#plt.setp(xticklabels, visible=False)

fig.savefig("prueba.png")
