import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


RS        = 6.957000e8
Year      = 365.25*86400
Myr       = Year*1.0e6
Gyr       = Year*1.0e9
Day       = 86400.0
InDeg     = 180.0/np.pi
AU = 149.6e9

radius      = "/Users/bportilla/SecDev3B_V2.0/radius.dat"

data_file_1 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b.dat"
data_file_2 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b_SE.dat"

info_file_1 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b.log"
info_file_2 = "/Users/bportilla/SecDev3B_V2.0/HD_41004_A_b_SE.log"


data_radius = np.loadtxt(radius)
time_radius = data_radius[:,0:1]
radius_valu = data_radius[:,1:2]

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



############################################################
# Fancy style
plt.rc('font', **{'family': 'family', 'family': ['serif']})
plt.rc('text', usetex=True) #uses Latex instead of Tex to compile axes labels


############################################################
# Creating axes

fig=plt.figure(figsize=(8,10))

gs_1 = gridspec.GridSpec(3,2,hspace=0.4)
ax_1=plt.subplot(gs_1[0, :])

gs_2 = gridspec.GridSpec(3,2,hspace=0.0,wspace=0.0)
ax_2=plt.subplot(gs_2[1, 0:1])
ax_3=plt.subplot(gs_2[1, 1:2])
ax_4=plt.subplot(gs_2[2, 0:1])
ax_5=plt.subplot(gs_2[2, 1:2])


############################################################
# Labels

labels=15
ax_1.set_xlabel("$t$ [Myr]",fontsize=labels)
ax_1.set_ylabel("$R \, [R_{\odot}]$",fontsize=labels)
ax_2.set_ylabel("$a$ [AU]",fontsize=labels)
ax_4.set_ylabel("$e_1$",fontsize=labels)
ax_4.set_xlabel("$t$ [Myr]",fontsize=labels)
ax_5.set_xlabel("$t$ [Myr]",fontsize=labels)


############################################################
# Create line of constant radius

radius_cte = []
for i in range(0,len(time_radius)):
    radius_cte.append(0.696)
    

############################################################
# Plots

ax_1.plot(time_radius*Gyr/Myr,radius_valu,linewidth=0.5,color='red')
ax_1.plot(time_radius*Gyr/Myr,radius_cte,linewidth=0.5,color='blue')
ax_2.plot(time_1*uT/Myr,a_in_1,linewidth=0.5,color='blue')
ax_3.plot(time_2*uT/Myr,a_in_2,linewidth=0.5,color='red')
ax_4.plot(time_1*uT/Myr,e_in_1,linewidth=0.5,color='blue')
ax_5.plot(time_2*uT/Myr,e_in_2,linewidth=0.5,color='red')



############################################################
# Create inset plot

x_inset_min = min(time_2)*uT/Myr
x_inset_max = max(time_2)*uT/Myr
y_inset_min = 1
y_inset_max = 1.8

ax_1_inset = zoomed_inset_axes(ax_1, 10.5, loc=1) # zoom-factor: 2.5, location: upper-left
ax_1_inset.plot(time_2*uT/Myr,R_A_2*uL/RS,color='red')
ax_1_inset.set_xlim(x_inset_min,x_inset_max)
ax_1_inset.tick_params(direction='in')
ax_1_inset.set_xlabel("$t$ [Myr]")
ax_1_inset.set_ylabel("$R \, [R_{\odot}]$")

mark_inset(ax_1, ax_1_inset, loc1=2, loc2=3, fc="none", ec="0.5")

############################################################
# Set scales

ax_1.set_xscale('log')
ax_1.set_yscale('log')


############################################################
# Ticks

ax_1.tick_params(direction='in')
ax_2.tick_params(direction='in')
ax_4.tick_params(direction='in')
ax_3.tick_params(direction='in')
ax_5.tick_params(direction='in')
ax_2.tick_params(right=False)
ax_4.tick_params(right=False)
ax_3.tick_params(labelleft=False)
ax_5.tick_params(labelleft=False)


############################################################
# Format of tick labels


y_tick_labels = ['0.6','1.0','2.0']
y_ticks       = [0.6,1.0,2.0]

ax_1.set_yticklabels(y_tick_labels)
ax_1.set_yticks(y_ticks)

ax_1.yaxis.set_minor_formatter(ticker.FormatStrFormatter('%0.1f'))
ax_1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax_1.minorticks_off()




############################################################
# Limits

ax_3.set_xlim(ax_2.get_xlim())
ax_5.set_xlim(ax_4.get_xlim())
ax_1.set_ylim(0.5,2.5)


############################################################
# axvspan

ax_3.axvspan(max(time_2*uT/Myr),max(time_1*uT/Myr),0.0,1.7,facecolor='cornsilk')
ax_5.axvspan(max(time_2*uT/Myr),max(time_1*uT/Myr),0.0,1.7,facecolor='cornsilk')


############################################################
# Plot arrow

x_ini = max(time_2*uT/Myr)
y_ini = 1.0
dx    = max(time_1*uT/Myr)-max(time_2*uT/Myr)
dy    = 0.0

ax_3.annotate(s='', xy=(x_ini+dx,y_ini), xytext=(x_ini,y_ini), arrowprops=dict(arrowstyle='<->'))
ax_5.annotate(s='', xy=(x_ini+dx,0.6), xytext=(x_ini,0.6), arrowprops=dict(arrowstyle='<->'))


############################################################
# Add text over arrow

ax_3.text(1.4,1.06, r'$0.47$ Myr')
ax_5.text(1.4,0.64,r'$0.47$ Myr')

dt = max(time_1*uT/Myr) - max(time_2*uT/Myr)
print dt

fig.savefig("prueba.png")


