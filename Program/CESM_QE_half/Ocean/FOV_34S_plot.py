#Program plots the FovS and determines early warning signals

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
import statsmodels.api as sm

#Making pathway to folder with all data
directory	      = '../../../Data/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV	= fh.variables['F_OV'][:]	#AMOC strength (Sv)

	fh.close()

	return time, FOV

def EWS_Classical(time_all, data_all, window = 50):
	"""Find the classical early warning signals for variance and (lag-1) autocorrelation"""


	auto_corr	= ma.masked_all(len(time_all))
	var		= ma.masked_all(len(time_all))

	for time_i in range(len(time_all) - window + 1):
		#Loop over each time window
		time_data	= time_all[time_i:time_i+window]
		data		= data_all[time_i:time_i+window]
		
		trend, base	= np.polyfit(time_data, data, 1)
		data		= data - ((trend * time_data) + base)

		#Now determine the auto-correlation and variance
		auto_corr[time_i + int(window / 2)]	= np.corrcoef(data[:-1], data[1:])[0,1]
		var[time_i + int(window / 2)]		= np.var(data)

	return auto_corr, var
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time_1, FOV_1	= ReadinData(directory+'CESM_QE/Ocean/FOV_section_34S.nc')	
time_2, FOV_2	= ReadinData(directory+'CESM_QE_half/Ocean/FOV_section_34S.nc')

time_1500, FOV_1500	= ReadinData(directory+'CESM_1500/Ocean/FOV_section_34S_branch_QE.nc')
time_1550, FOV_1550	= ReadinData(directory+'CESM_1550/Ocean/FOV_section_34S_branch_QE.nc')	
time_1600, FOV_1600	= ReadinData(directory+'CESM_1600/Ocean/FOV_section_34S_branch_1550.nc')	
time_1650, FOV_1650	= ReadinData(directory+'CESM_1650/Ocean/FOV_section_34S_branch_1600.nc')				

time_1, FOV_1		= time_1[:2200], FOV_1[:2200]
#-----------------------------------------------------------------------------------------

#Classical Early Warning Signals
auto_corr_1, var_1		    = EWS_Classical(time_1, FOV_1)
auto_corr_2, var_2		    = EWS_Classical(time_2, FOV_2)
auto_corr_1500, var_1500	= EWS_Classical(time_1500[-50:], FOV_1500[-50:])
auto_corr_1550, var_1550	= EWS_Classical(time_1550[-50:], FOV_1550[-50:])
auto_corr_1600, var_1600	= EWS_Classical(time_1600[-50:], FOV_1600[-50:])
auto_corr_1650, var_1650	= EWS_Classical(time_1650[-50:], FOV_1650[-50:])

#Convert to freshwater flux forcing
time_1				= time_1 * 0.0003
time_2				= (time_2 - 1999) * 0.00015 + 0.45

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

graph_1	= plot(time_1, FOV_1, '-k', linewidth = 0.5, alpha = 0.5)
graph_2	= plot(time_2, FOV_2, '-r', linewidth = 0.5)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(0.40, 0.66)
ax.set_ylim(-0.25, 0.25)
ax.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
ax.grid()

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

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', alpha = 0.5, linewidth = 1.5, label = r'$\partial_t F_H = 3.0 \times 10^{-4}$ Sv yr$^{-1}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = r'$\partial_t F_H = 1.5 \times 10^{-4}$ Sv yr$^{-1}$')

legend_1	= ax.legend(loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)


ax.set_title('b) Freshwater transport at 34$^{\circ}$S, $F_{\mathrm{ovS}}$')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1	= plot(time_1, var_1 * 10**4.0, '-k', linewidth = 0.75, alpha = 0.5)
graph_2	= plot(time_2, var_2 * 10**4.0, '-r', linewidth = 0.75)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', alpha = 0.5, linewidth = 1.5, label = r'$\partial_t F_H = 3.0 \times 10^{-4}$ Sv yr$^{-1}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = r'$\partial_t F_H = 1.5 \times 10^{-4}$ Sv yr$^{-1}$')

graph_3	= ax.scatter(1500 * 0.0003, np.mean(var_1500) * 10**4.0, s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10, label = 'Steady state')
ax.scatter(1550 * 0.0003, np.mean(var_1550) * 10**4.0, s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)
ax.scatter(1600 * 0.0003, np.mean(var_1600) * 10**4.0, s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)
ax.scatter(1650 * 0.0003, np.mean(var_1650) * 10**4.0, s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel(r'Freshwater transport variance ($\times 10^{-4}$ Sv$^2$)')
ax.set_xlim(0.40, 0.66)
ax.set_ylim(0, 10)
ax.grid()

legend_1	= ax.legend(loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) $F_{\mathrm{ovS}}$, variance')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

graph_1	= plot(time_1, auto_corr_1, '-k', linewidth = 0.75, alpha = 0.5)
graph_2	= plot(time_2, auto_corr_2, '-r', linewidth = 0.75)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', alpha = 0.5, linewidth = 1.5, label = r'$\partial_t F_H = 3.0 \times 10^{-4}$ Sv yr$^{-1}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = r'$\partial_t F_H = 1.5 \times 10^{-4}$ Sv yr$^{-1}$')

graph_3	= ax.scatter(1500 * 0.0003, np.mean(auto_corr_1500), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10, label = 'Steady state')
ax.scatter(1550 * 0.0003, np.mean(auto_corr_1550), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)
ax.scatter(1600 * 0.0003, np.mean(auto_corr_1600), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)
ax.scatter(1650 * 0.0003, np.mean(auto_corr_1650), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Autocorrelation')
ax.set_xlim(0.40, 0.66)
ax.set_ylim(-1, 1)
ax.grid()


legend_1	= ax.legend(loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('f) $F_{\mathrm{ovS}}$, autocorrelation')

show()


