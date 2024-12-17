#Program plots the different gyre and FovN responses under the varying freshwater flux forcing

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/CESM_QE/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------
	
	
fh = netcdf.Dataset(directory+'Ocean/Freshwater_responses_Fgyre_FovN.nc', 'r')

time		= fh.variables['time'][:1700]   
F_ovN		= fh.variables['FovN'][:1700]    		
F_gyre		= fh.variables['Fgyre'][:1700]    		
salt_limb	= fh.variables['SALT_limb'][:1700]  

fh.close()

#-----------------------------------------------------------------------------------------

year_average 		= 20

time_average    	= int(len(time) / year_average)
time_average		= np.zeros(time_average)

F_ovN_average     	= np.zeros(len(time_average))
F_gyre_average     	= np.zeros(len(time_average))
salt_limb_average     	= np.zeros(len(time_average))

for time_i in range(len(time_average)):
	#Loop over each time slice
	time_average[time_i]		= np.mean(time[time_i*year_average:(time_i+1)*year_average])
	F_ovN_average[time_i] 		= np.mean(F_ovN[time_i*year_average:(time_i+1)*year_average])
	F_gyre_average[time_i] 		= np.mean(F_gyre[time_i*year_average:(time_i+1)*year_average])
	salt_limb_average[time_i] 	= np.mean(salt_limb[time_i*year_average:(time_i+1)*year_average])

#-----------------------------------------------------------------------------------------	

n_1, n_2	= np.polyfit(salt_limb_average, F_ovN_average, 1)
salt_fit	= np.arange(-0.6, 0.61, 0.1)

fig, ax	= subplots()

ax.plot(salt_fit, salt_fit * n_1 + n_2, '--k')

CS	= ax.scatter(salt_limb_average, F_ovN_average, c=time_average*0.0003, s = 40, vmin=0,vmax=np.max(time_average*0.0003),cmap = 'Spectral_r',zorder=10, edgecolor='k')
cbar	= colorbar(CS, ticks = np.arange(0, 0.51, 0.1))
cbar.set_label('Freshwater flux forcing, $F_H$ (Sv)')

ax.set_xlabel(r'$S_{\rightleftarrows}$ (g kg$^{-1}$)')
ax.set_ylabel('Freshwater transport by $F_{\mathrm{ovN}}$ (Sv)')
ax.set_xlim(-0.5, 0.5)
ax.set_ylim(-0.05, 0.0025)
ax.grid()

ax.set_title('b) Freshwater transport by $F_{\mathrm{ovN}}$')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

g_1, g_2	= np.polyfit(salt_limb_average, F_gyre_average, 1)

fig, ax	= subplots()

ax.plot(salt_fit, salt_fit * g_1 + g_2, '--k')

CS	= ax.scatter(salt_limb_average, F_gyre_average, c=time_average*0.0003, s = 40, vmin=0,vmax=np.max(time_average*0.0003),cmap = 'Spectral_r',zorder=10, edgecolor='k')
cbar	= colorbar(CS, ticks = np.arange(0, 0.51, 0.1))
cbar.set_label('Freshwater flux forcing, $F_H$ (Sv)')

ax.set_xlabel(r'$S_{\rightleftarrows}$ (g kg$^{-1}$)')
ax.set_ylabel('Freshwater transport by $F_{\mathrm{gyre}}$ (Sv)')
ax.set_xlim(-0.5, 0.5)
ax.set_ylim(0.4, 0.6)
ax.set_yticks([0.4, 0.45, 0.5, 0.55, 0.6])
ax.grid()

ax.set_title('a) Freshwater transport by $F_{\mathrm{gyre}}$')

show()


  	
