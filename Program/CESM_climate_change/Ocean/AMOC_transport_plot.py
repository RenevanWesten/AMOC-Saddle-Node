#Program plots the AMOC strength under climate change

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
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

time_hist_rcp45_branch_0600, transport_hist_rcp45_branch_0600		= ReadinData(directory+'CESM_0600_climate_change/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_RCP45.nc')
time_hist_rcp85_branch_0600, transport_hist_rcp85_branch_0600		= ReadinData(directory+'CESM_0600_climate_change/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_RCP85.nc')

time_hist_rcp45_branch_1500, transport_hist_rcp45_branch_1500		= ReadinData(directory+'CESM_1500_climate_change/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_RCP45.nc')
time_hist_rcp85_branch_1500, transport_hist_rcp85_branch_1500		= ReadinData(directory+'CESM_1500_climate_change/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_RCP85.nc')

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------	


fig, ax	= subplots()

ax.fill_between([-100, 5000], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_hist_rcp45_branch_0600[155:], transport_hist_rcp45_branch_0600[155:], '-', color = 'dodgerblue', linewidth = 1)
ax.plot(time_hist_rcp85_branch_0600[155:], transport_hist_rcp85_branch_0600[155:], '-', color = 'firebrick', linewidth = 1)
ax.plot(time_hist_rcp85_branch_0600[:156], transport_hist_rcp85_branch_0600[:156], '-', color = 'k', linewidth = 1)

ax.set_xlim(1850, 2500)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xticks([1900, 2000, 2100, 2200, 2300, 2400, 2500])
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = 'Historical')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'dodgerblue', linewidth = 2, label = 'RCP4.5')
graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = 'RCP8.5')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) AMOC strength at 26$^{\circ}$N ($\overline{F_H} = 0.18$ Sv)')

#-----------------------------------------------------------------------------------------	


fig, ax	= subplots()

ax.fill_between([-100, 5000], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(time_hist_rcp45_branch_1500[155:], transport_hist_rcp45_branch_1500[155:], '-', color = 'dodgerblue', linewidth = 1)
ax.plot(time_hist_rcp85_branch_1500[155:], transport_hist_rcp85_branch_1500[155:], '-', color = 'firebrick', linewidth = 1)
ax.plot(time_hist_rcp85_branch_1500[:156], transport_hist_rcp85_branch_1500[:156], '-', color = 'k', linewidth = 1)

ax.set_xlim(1850, 2500)
ax.set_ylim(-2, 22)
ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xticks([1900, 2000, 2100, 2200, 2300, 2400, 2500])
ax.grid()

ax.set_title('b) AMOC strength at 26$^{\circ}$N ($\overline{F_H} = 0.45$ Sv)')

graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 2, label = 'Historical')
graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'dodgerblue', linewidth = 2, label = 'RCP4.5')
graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 2, label = 'RCP8.5')

graphs	      	= graph_1 + graph_2 + graph_3
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

show()

