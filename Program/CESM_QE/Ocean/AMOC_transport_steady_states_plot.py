#Program plots the AMOC strength and the various steady states

from pylab import *
import numpy
import time
import glob, os
import netCDF4 as netcdf
import matplotlib.colors as colors
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

time, transport			    = ReadinData(directory+'CESM_QE/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')	

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

for time_i in [600, 1500, 1550, 1600, 1650, 1700]:
	#Loop over each branch

	time_branch_QE, transport_branch_QE	= ReadinData(directory+'CESM_'+str(time_i).zfill(4)+'/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')

	fig, ax	= subplots()

	ax.plot(time, transport, '-k', linewidth = 0.75, zorder = 6)
	ax.plot(time_branch_QE, transport_branch_QE, '-', linewidth = 0.75, color = 'royalblue', zorder = 7)


	if time_i == 1600 or time_i == 1650 or time_i == 1700:
        #Read in the simulations initiated from the previous steady state
		time_branch_prev, transport_branch_prev	= ReadinData(directory+'CESM_'+str(time_i).zfill(4)+'/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_'+str(time_i-50).zfill(4)+'.nc')

		ax.plot(time_branch_prev, transport_branch_prev, '-', linewidth = 0.75, color = 'firebrick', zorder = 7)
		
		
	ax.set_xlim(1480, 2220)
	ax.set_ylim(-2, 18)
	ax.set_yticks([0, 5, 10, 15])
	ax.set_xlabel('Model year')
	ax.set_ylabel('Volume transport (Sv)')
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

	legend_1	= ax.legend(loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(5)

	if time_i == 600: subplot_letter='c)'
	if time_i == 1500: subplot_letter='d)'
	if time_i == 1550: subplot_letter='e)'
	if time_i == 1600: subplot_letter='f)'
	if time_i == 1650: subplot_letter='g)'
	if time_i == 1700: subplot_letter='h)'
	
	ax.set_title(subplot_letter+' AMOC strength, $\overline{F_H} = $'+str(np.round(time_i*0.0003, 3))+' Sv (model year '+str(time_i)+')')
	

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

time_0600, transport_0600	= ReadinData(directory+'CESM_0600/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')	
time_1500, transport_1500	= ReadinData(directory+'CESM_1500/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')
time_1550, transport_1550	= ReadinData(directory+'CESM_1550/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')	
time_1600, transport_1600	= ReadinData(directory+'CESM_1600/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_1550.nc')	
time_1650, transport_1650	= ReadinData(directory+'CESM_1650/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_1600.nc')	

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

ax.fill_between([-100, 2500], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')


ax.plot((4400 - time[2199:])*0.0003, transport[2199:], '-r', linewidth = 0.5)
ax.plot(time[:2200]*0.0003, transport[:2200], '-k', linewidth = 0.5)

y_error_0600		= np.zeros((2, 1)) + 1
y_error_0600[:, 0]	= np.mean(transport_0600[-50:]) - np.min(transport_0600[-50:]), np.max(transport_0600[-50:]) - np.mean(transport_0600[-50:])
ax.errorbar(600*0.0003, np.mean(transport_0600[-50:]), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_0600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

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


ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(0, 0.66)
ax.set_ylim(-2, 35)
ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
ax.grid()

ax.plot([0.205, 0.295], [22.5, 22.5], '-k')
ax.plot([0.205, 0.205], [22, 23], '-k')
ax.plot([0.295, 0.295], [22, 23], '-k')
ax.text(0.250, 24, '300 years', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10)

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Forward quasi-equilibrium')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = 'Backward quasi-equilibrium')

legend_1	= ax.legend(loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) AMOC strength at 26$^{\circ}$N')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.63, 0.40, 0.35, 0.35], projection = ccrs.Orthographic(-30, 10))

ax2.coastlines(resolution='50m')
ax2.gridlines()
ax2.add_feature(cfeature.LAND, zorder=10)
ax2.set_global()


lon     = np.arange(-1, 360)
lat     = np.arange(-90, 91)
field   = np.ones((len(lat), len(lon))) * -0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())


lon     = np.arange(-100, -5)
lat     = np.arange(20, 43)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

lon     = np.arange(-100, 3)
lat     = np.arange(42, 51)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax2.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

ax2.text(320, 38, '$+F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())
ax2.text(340, -10, '$-F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=11, transform=ccrs.PlateCarree())

x_1	= np.arange(-81, -9.99, 0.1)
y_1	= np.zeros(len(x_1)) + 26.0
y_2	= np.arange(24, 28.01, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(24, 28.01, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

x_1	= np.arange(-60, 20.01, 0.1)
y_1	= np.zeros(len(x_1)) - 34
y_2	= np.arange(-37, -30.99, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(-37, -30.99, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)


show()
