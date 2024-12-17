#Program plots the 2-meter surface temperature trend

from pylab import *
import numpy
import time
import glob, os
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinDataGMST(filename, year_start, year_end):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	temp	= fh.variables['TEMP_global'][:]	#Global mean surface temperature

	fh.close()
	
	time_start	= (np.abs(time - year_start)).argmin()
	time_end	= (np.abs(time - (year_end+1))).argmin()

	time		= time[time_start:time_end]
	temp		= temp[time_start:time_end]
	
	#Get the temperature trend
	a, b		= np.polyfit(time, temp, 1)

	return a * 100
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

year_start	= 2000
year_end	= 2100

#-----------------------------------------------------------------------------------------

temp_trend_global_rcp45_branch_0600	= ReadinDataGMST(directory+'CESM_0600_climate_change/Atmosphere/TEMP_2m_RCP45.nc', year_start, year_end)
temp_trend_global_rcp85_branch_0600	= ReadinDataGMST(directory+'CESM_0600_climate_change/Atmosphere/TEMP_2m_RCP85.nc', year_start, year_end)

temp_trend_global_rcp45_branch_1500	= ReadinDataGMST(directory+'CESM_1500_climate_change/Atmosphere/TEMP_2m_RCP45.nc', year_start, year_end)
temp_trend_global_rcp85_branch_1500	= ReadinDataGMST(directory+'CESM_1500_climate_change/Atmosphere/TEMP_2m_RCP85.nc', year_start, year_end)

#-----------------------------------------------------------------------------------------

fh 				= netcdf.Dataset(directory+'CESM_0600_climate_change/Atmosphere/TEMP_2m_trend_year_'+str(year_start)+'-'+str(year_end)+'_RCP45.nc', 'r')

lat				= fh.variables['lat'][:] 			
temp_trend_rcp45_branch_0600	= fh.variables['TEMP_trend'][:]

fh.close()

#-----------------------------------------------------------------------------------------

fh 				= netcdf.Dataset(directory+'CESM_0600_climate_change/Atmosphere/TEMP_2m_trend_year_'+str(year_start)+'-'+str(year_end)+'_RCP85.nc', 'r')

temp_trend_rcp85_branch_0600	= fh.variables['TEMP_trend'][:] 

fh.close()


#-----------------------------------------------------------------------------------------

fh 				= netcdf.Dataset(directory+'CESM_1500_climate_change/Atmosphere/TEMP_2m_trend_year_'+str(year_start)+'-'+str(year_end)+'_RCP45.nc', 'r')

temp_trend_rcp45_branch_1500	= fh.variables['TEMP_trend'][:]

fh.close()

#-----------------------------------------------------------------------------------------


fh 				= netcdf.Dataset(directory+'CESM_1500_climate_change/Atmosphere/TEMP_2m_trend_year_'+str(year_start)+'-'+str(year_end)+'_RCP85.nc', 'r')

temp_trend_rcp85_branch_1500	= fh.variables['TEMP_trend'][:] 

fh.close()

#-----------------------------------------------------------------------------------------
#Take the zonal means of the trends
temp_trend_rcp45_branch_0600	= np.mean(temp_trend_rcp45_branch_0600, axis = 1) 
temp_trend_rcp45_branch_1500	= np.mean(temp_trend_rcp45_branch_1500, axis = 1)
temp_trend_rcp85_branch_0600	= np.mean(temp_trend_rcp85_branch_0600, axis = 1) 
temp_trend_rcp85_branch_1500	= np.mean(temp_trend_rcp85_branch_1500, axis = 1) 
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1 = ax.plot(lat, temp_trend_rcp45_branch_0600, '-', linewidth = 2.0, color = 'dodgerblue', label = 'Hist/RCP4.5')
graph_2	= ax.plot(lat, temp_trend_rcp85_branch_0600, '-', linewidth = 2.0, color = 'firebrick', label = 'Hist/RCP8.5')
ax.plot(lat, np.zeros(len(lat))+temp_trend_global_rcp45_branch_0600, '--', linewidth = 2.0, color = 'dodgerblue')
ax.plot(lat, np.zeros(len(lat))+temp_trend_global_rcp85_branch_0600, '--', linewidth = 2.0, color = 'firebrick')

ax.set_ylabel('Temperature trend ($^{\circ}$C per century)')
ax.set_xlim(-90, 90)
ax.set_ylim(0, 9)
ax.grid()

ax.set_xticks(np.arange(-80, 80.1, 20))
ax.set_xticklabels(['80$^{\circ}$S', '60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N'])

graphs	      	= graph_1 + graph_2
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) Zonally-averaged temperature trend, 2000 - 2100 ($\overline{F_H} = 0.18$ Sv)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1 = ax.plot(lat, temp_trend_rcp45_branch_1500, '-', linewidth = 2.0, color = 'dodgerblue', label = 'Hist/RCP4.5')
graph_2	= ax.plot(lat, temp_trend_rcp85_branch_1500, '-', linewidth = 2.0, color = 'firebrick', label = 'Hist/RCP8.5')
ax.plot(lat, np.zeros(len(lat))+temp_trend_global_rcp45_branch_1500, '--', linewidth = 2.0, color = 'dodgerblue')
ax.plot(lat, np.zeros(len(lat))+temp_trend_global_rcp85_branch_1500, '--', linewidth = 2.0, color = 'firebrick')

ax.set_ylabel('Temperature trend ($^{\circ}$C per century)')
ax.set_xlim(-90, 90)
ax.set_ylim(0, 9)
ax.grid()

ax.set_xticks(np.arange(-80, 80.1, 20))
ax.set_xticklabels(['80$^{\circ}$S', '60$^{\circ}$S', '40$^{\circ}$S', '20$^{\circ}$S', 'Eq', '20$^{\circ}$N', '40$^{\circ}$N', '60$^{\circ}$N', '80$^{\circ}$N'])

graphs	      	= graph_1 + graph_2
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) Zonally-averaged temperature trend, 2000 - 2100 ($\overline{F_H} = 0.45$ Sv)')

show()


