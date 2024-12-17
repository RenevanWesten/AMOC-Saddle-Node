#Program plots the AMOC strength and determines early warning signals

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

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport

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

depth_min 	= 0
depth_max	= 1000

time_1, transport_1		    = ReadinData(directory+'CESM_QE/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	
time_2, transport_2		    = ReadinData(directory+'CESM_QE_half/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	

time_1500, transport_1500	= ReadinData(directory+'CESM_1500/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')
time_1550, transport_1550	= ReadinData(directory+'CESM_1550/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')	
time_1600, transport_1600	= ReadinData(directory+'CESM_1600/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_1550.nc')	
time_1650, transport_1650	= ReadinData(directory+'CESM_1650/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_1600.nc')

time_1, transport_1		    = time_1[:2200], transport_1[:2200]
#-----------------------------------------------------------------------------------------

#Classical Early Warning Signals
auto_corr_1, var_1		    = EWS_Classical(time_1, transport_1)
auto_corr_2, var_2		    = EWS_Classical(time_2, transport_2)
auto_corr_1500, var_1500	= EWS_Classical(time_1500[-50:], transport_1500[-50:])
auto_corr_1550, var_1550	= EWS_Classical(time_1550[-50:], transport_1550[-50:])
auto_corr_1600, var_1600	= EWS_Classical(time_1600[-50:], transport_1600[-50:])
auto_corr_1650, var_1650	= EWS_Classical(time_1650[-50:], transport_1650[-50:])

#Convert to freshwater flux forcing
time_1				        = time_1 * 0.0003
time_2				        = (time_2 - 1999) * 0.00015 + 0.45

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

graph_1	= plot(time_1, transport_1, '-k', linewidth = 0.5, alpha = 0.5)
graph_2	= plot(time_2, transport_2, '-', color = 'r', linewidth = 0.5)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(0.40, 0.66)
ax.set_ylim(-2, 18)
ax.set_yticks([0, 5, 10, 15])
ax.grid()

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(transport_1500[-50:]) - np.min(transport_1500[-50:]), np.max(transport_1500[-50:]) - np.mean(transport_1500[-50:])
graph_3 		= ax.errorbar(1500*0.0003, np.mean(transport_1500[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')

y_error_1550		= np.zeros((2, 1)) + 1
y_error_1550[:, 0]	= np.mean(transport_1550[-50:]) - np.min(transport_1550[-50:]), np.max(transport_1550[-50:]) - np.mean(transport_1550[-50:])
graph_3 		= ax.errorbar(1550*0.0003, np.mean(transport_1550[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1550, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1600		= np.zeros((2, 1)) + 1
y_error_1600[:, 0]	= np.mean(transport_1600[-50:]) - np.min(transport_1600[-50:]), np.max(transport_1600[-50:]) - np.mean(transport_1600[-50:])
graph_3 		= ax.errorbar(1600*0.0003, np.mean(transport_1600[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1650		= np.zeros((2, 1)) + 1
y_error_1650[:, 0]	= np.mean(transport_1650[-50:]) - np.min(transport_1650[-50:]), np.max(transport_1650[-50:]) - np.mean(transport_1650[-50:])
graph_3 		= ax.errorbar(1650*0.0003, np.mean(transport_1650[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1650, linewidth = 0.0, elinewidth = 2.0, capsize=4)


ax.plot([0.405, 0.495], [2.5, 2.5], '-k')
ax.plot([0.405, 0.405], [2.1, 2.9], '-k')
ax.plot([0.495, 0.495], [2.1, 2.9], '-k')
ax.text(0.450, 3.1, '300 years', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 11)

ax.plot([0.4275, 0.4725], [7.5, 7.5], '-', color = 'r')
ax.plot([0.4275, 0.4275], [7.1, 7.9], '-', color = 'r')
ax.plot([0.4725, 0.4725], [7.1, 7.9], '-', color = 'r')
ax.text(0.450, 8.1, '300 years', verticalalignment='center', horizontalalignment='center', color = 'r', fontsize = 11)


graph_1		= ax.plot([-100, -100], [-100, -100], '-k', alpha = 0.5, linewidth = 1.5, label = r'$\partial_t F_H = 3.0 \times 10^{-4}$ Sv yr$^{-1}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = r'$\partial_t F_H = 1.5 \times 10^{-4}$ Sv yr$^{-1}$')

legend_1	= ax.legend(loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) AMOC strength at 26$^{\circ}$N')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1	= plot(time_1, var_1, '-k', linewidth = 0.75, alpha = 0.5)
graph_2	= plot(time_2, var_2, '-r', linewidth = 0.75)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', alpha = 0.5, linewidth = 1.5, label = r'$\partial_t F_H = 3.0 \times 10^{-4}$ Sv yr$^{-1}$')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = r'$\partial_t F_H = 1.5 \times 10^{-4}$ Sv yr$^{-1}$')

graph_3	= ax.scatter(1500 * 0.0003, np.mean(var_1500), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10, label = 'Steady state')
ax.scatter(1550 * 0.0003, np.mean(var_1550), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)
ax.scatter(1600 * 0.0003, np.mean(var_1600), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)
ax.scatter(1650 * 0.0003, np.mean(var_1650), s = 40, color = 'darkorchid', marker = 's', edgecolor = 'mediumpurple', zorder = 10)

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel(r'Volume transport variance (Sv$^2$)')
ax.set_xlim(0.40, 0.66)
ax.set_ylim(0, 0.60)
ax.grid()

legend_1	= ax.legend(loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) AMOC strength, variance')

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

ax.set_title('e) AMOC strength, autocorrelation')

show()


