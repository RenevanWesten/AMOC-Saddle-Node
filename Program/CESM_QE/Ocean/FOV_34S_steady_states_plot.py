#Program plots the FovS and the various steady states

from pylab import *
import numpy
import time
import glob, os
import netCDF4 as netcdf
import matplotlib.colors as colors


#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	#Fresh water

	fh.close()

	return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

time, FOV		    = ReadinData(directory+'CESM_QE/Ocean/FOV_section_34S.nc')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

for time_i in [600, 1500, 1550, 1600, 1650, 1700]:
	#Loop over each branch

	time_branch_QE, FOV_branch_QE	= ReadinData(directory+'CESM_'+str(time_i).zfill(4)+'/Ocean/FOV_section_34S_branch_QE.nc')


	fig, ax	= subplots()

	ax.plot(time, FOV, '-k', linewidth = 0.75, zorder = 6)
	ax.plot(time_branch_QE, FOV_branch_QE, '-', linewidth = 0.75, color = 'royalblue', zorder = 7)


	if time_i == 1600 or time_i == 1650 or time_i == 1700:
        #Read in the simulations initiated from the previous steady state
		time_branch_prev, FOV_branch_prev	= ReadinData(directory+'CESM_'+str(time_i).zfill(4)+'/Ocean/FOV_section_34S_branch_'+str(time_i-50).zfill(4)+'.nc')

		ax.plot(time_branch_prev, FOV_branch_prev, '-', linewidth = 0.75, color = 'firebrick', zorder = 7)

	ax.set_xlim(1480, 2220)
	ax.set_ylim(-0.25, 0.25)
	ax.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
	ax.set_xlabel('Model year')
	ax.set_ylabel('Freshwater transport (Sv)')
	ax.grid()

	if time_i == 600:
		ax.set_xlim(580, 1320)	
		
	if time_i == 600 or time_i == 1500 or time_i == 1550 or time_i == 1750:
		graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 1.5, label = 'Quasi-equilibrium (QE)')
		graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 1.5, label = 'Branched from QE')


	else:
		graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 1.5, label = 'Quasi-equilibrium (QE)')
		graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 1.5, label = 'Branched from QE')
		graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'Branched from $\overline{F_H} = $'+str(np.round((time_i-50)*0.0003, 3))+' Sv')

	legend_1	= ax.legend(loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(5)

	if time_i == 600: subplot_letter='i)'
	if time_i == 1500: subplot_letter='j)'
	if time_i == 1550: subplot_letter='k)'
	if time_i == 1600: subplot_letter='l)'
	if time_i == 1650: subplot_letter='m)'
	if time_i == 1700: subplot_letter='n)'
	
	ax.set_title(subplot_letter+' $F_{\mathrm{ovS}}$, $\overline{F_H} = $'+str(np.round(time_i*0.0003, 3))+' Sv (model year '+str(time_i)+')')
	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

time_0600, FOV_0600	= ReadinData(directory+'CESM_0600/Ocean/FOV_section_34S_branch_QE.nc')	
time_1500, FOV_1500	= ReadinData(directory+'CESM_1500/Ocean/FOV_section_34S_branch_QE.nc')
time_1550, FOV_1550	= ReadinData(directory+'CESM_1550/Ocean/FOV_section_34S_branch_QE.nc')	
time_1600, FOV_1600	= ReadinData(directory+'CESM_1600/Ocean/FOV_section_34S_branch_1550.nc')	
time_1650, FOV_1650	= ReadinData(directory+'CESM_1650/Ocean/FOV_section_34S_branch_1600.nc')	

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-100, 2500], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')


ax.plot((4400 - time[2199:])*0.0003, FOV[2199:], '-r', linewidth = 0.5)
ax.plot(time[:2200]*0.0003, FOV[:2200], '-k', linewidth = 0.5)

y_error_0600		= np.zeros((2, 1)) + 1
y_error_0600[:, 0]	= np.mean(FOV_0600[-50:]) - np.min(FOV_0600[-50:]), np.max(FOV_0600[-50:]) - np.mean(FOV_0600[-50:])
ax.errorbar(600*0.0003, np.mean(FOV_0600[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_0600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(FOV_1500[-50:]) - np.min(FOV_1500[-50:]), np.max(FOV_1500[-50:]) - np.mean(FOV_1500[-50:])
graph_3 		= ax.errorbar(1500*0.0003, np.mean(FOV_1500[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')

y_error_1550		= np.zeros((2, 1)) + 1
y_error_1550[:, 0]	= np.mean(FOV_1550[-50:]) - np.min(FOV_1550[-50:]), np.max(FOV_1550[-50:]) - np.mean(FOV_1550[-50:])
graph_3 		= ax.errorbar(1550*0.0003, np.mean(FOV_1550[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1550, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1600		= np.zeros((2, 1)) + 1
y_error_1600[:, 0]	= np.mean(FOV_1600[-50:]) - np.min(FOV_1600[-50:]), np.max(FOV_1600[-50:]) - np.mean(FOV_1600[-50:])
graph_3 		= ax.errorbar(1600*0.0003, np.mean(FOV_1600[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1650		= np.zeros((2, 1)) + 1
y_error_1650[:, 0]	= np.mean(FOV_1650[-50:]) - np.min(FOV_1650[-50:]), np.max(FOV_1650[-50:]) - np.mean(FOV_1650[-50:])
graph_3 		= ax.errorbar(1650*0.0003, np.mean(FOV_1650[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1650, linewidth = 0.0, elinewidth = 2.0, capsize=4)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(0, 0.66)
ax.set_ylim(-0.35, 0.35)
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Forward quasi-equilibrium')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = 'Backward quasi-equilibrium')

legend_1	= ax.legend(loc=(0.18, 0.015), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) Freshwater transport at 34$^{\circ}$S, $F_{\mathrm{ovS}}$')


show()
