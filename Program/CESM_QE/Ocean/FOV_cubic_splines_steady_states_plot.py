#Program determines the FovS minimum using cubic splines

from pylab import *
import numpy
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from scipy.interpolate import CubicSpline

#Making pathway to folder with all data
directory	= '../../../Data/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	#Fresh water

	fh.close()

	return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time_1500, FOV_1500	    = ReadinData(directory+'CESM_1500/Ocean/FOV_section_34S_branch_QE.nc')
time_1550, FOV_1550	    = ReadinData(directory+'CESM_1550/Ocean/FOV_section_34S_branch_QE.nc')	
time_1600, FOV_1600	    = ReadinData(directory+'CESM_1600/Ocean/FOV_section_34S_branch_1550.nc')	
time_1650, FOV_1650	    = ReadinData(directory+'CESM_1650/Ocean/FOV_section_34S_branch_1600.nc')	

#Only consider the last 50 model years
FOV_1500	= FOV_1500[-50:]
FOV_1550	= FOV_1550[-50:]
FOV_1600	= FOV_1600[-50:]
FOV_1650	= FOV_1650[-50:]	

print('FovS (1500): ', np.mean(FOV_1500))
print('FovS (1550): ', np.mean(FOV_1550))
print('FovS (1600): ', np.mean(FOV_1600))
print('FovS (1650): ', np.mean(FOV_1550))
print()
#-----------------------------------------------------------------------------------------	

F_H_knot	= np.asarray([0.45, 0.465, 0.480, 0.495])
F_H		    = np.arange(0.45, 0.49501, 0.001)

realisation	= 100000

FOV_mean_1     	= ma.masked_all((realisation, len(F_H)))
FOV_min_1	= np.zeros(len(F_H))
FOV_min_F_H_1	= ma.masked_all(realisation)
FOV_mean_2     	= ma.masked_all((realisation, len(F_H)))
FOV_min_2	= np.zeros(len(F_H))
FOV_min_F_H_2	= ma.masked_all(realisation)

FOV_min_1_knot	= np.zeros(len(F_H_knot))
FOV_min_2_knot	= np.zeros(len(F_H_knot))
FOV_min_3_knot	= np.zeros(len(F_H_knot))

for real_i in range(realisation):
	#Get the different knots
	
	for random_i in range(3):
		#Draw a random ensemble
		FOV_knot_1500	= np.random.choice(FOV_1500)
		FOV_knot_1550	= np.random.choice(FOV_1550)
		FOV_knot_1600	= np.random.choice(FOV_1600)
		FOV_knot_1650	= np.random.choice(FOV_1650)

		#Get the knotes, based on a random sample from the steady states
		if random_i == 0:
			FOV_knot_1        			= np.asarray([FOV_knot_1500, FOV_knot_1550, FOV_knot_1600, FOV_knot_1650])
		if random_i == 1:
			FOV_knot_2        			= np.asarray([FOV_knot_1500, FOV_knot_1550, FOV_knot_1600, FOV_knot_1650])
		if random_i == 2:
			FOV_knot_3        			= np.asarray([FOV_knot_1500, FOV_knot_1550, FOV_knot_1600, FOV_knot_1650])
			
	FOV_min_3_knot[np.argmin(FOV_knot_3)]	+= 1

	#Fit the cubic splines
	cs_1             = CubicSpline(F_H_knot, FOV_knot_1, bc_type = 'not-a-knot')
	cs_2             = CubicSpline(F_H_knot, FOV_knot_2, bc_type = 'natural')
	
	#Save the values
	FOV_mean_1[real_i]				= cs_1(F_H)
	FOV_min_1[np.argmin(FOV_mean_1[real_i])]	+= 1
	FOV_min_F_H_1[real_i]				= F_H[np.argmin(FOV_mean_1[real_i])]
	FOV_mean_2[real_i]				= cs_2(F_H)
	FOV_min_2[np.argmin(FOV_mean_2[real_i])]	+= 1
	FOV_min_F_H_2[real_i]				= F_H[np.argmin(FOV_mean_2[real_i])]
	
	F_H_1	= F_H[np.argmin(FOV_mean_1[real_i])]
	F_H_2	= F_H[np.argmin(FOV_mean_2[real_i])]
	
	#Get the accumulated sums
	FOV_min_1_knot[(np.abs(F_H_knot - F_H_1)).argmin()]	+= 1
	FOV_min_2_knot[(np.abs(F_H_knot - F_H_2)).argmin()]	+= 1
	

print('FovS minimum (not-a-knot): ', np.mean(FOV_min_F_H_1))
print('FovS minimum (natural): ', np.mean(FOV_min_F_H_2))

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, (ax1, ax2)		= subplots(2,1, sharex=True)

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(FOV_1500) - np.min(FOV_1500), np.max(FOV_1500) - np.mean(FOV_1500)
graph_3 		= ax1.errorbar(1500*0.0003, np.mean(FOV_1500), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')
graph_3 		= ax2.errorbar(1500*0.0003, np.mean(FOV_1500), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')

y_error_1550		= np.zeros((2, 1)) + 1
y_error_1550[:, 0]	= np.mean(FOV_1550) - np.min(FOV_1550), np.max(FOV_1550) - np.mean(FOV_1550)
graph_3 		= ax1.errorbar(1550*0.0003, np.mean(FOV_1550), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1550, linewidth = 0.0, elinewidth = 2.0, capsize=4)
graph_3 		= ax2.errorbar(1550*0.0003, np.mean(FOV_1550), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1550, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1600		= np.zeros((2, 1)) + 1
y_error_1600[:, 0]	= np.mean(FOV_1600) - np.min(FOV_1600), np.max(FOV_1600) - np.mean(FOV_1600[-50:])
graph_3 		= ax1.errorbar(1600*0.0003, np.mean(FOV_1600), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1600, linewidth = 0.0, elinewidth = 2.0, capsize=4)
graph_3 		= ax2.errorbar(1600*0.0003, np.mean(FOV_1600), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1650		= np.zeros((2, 1)) + 1
y_error_1650[:, 0]	= np.mean(FOV_1650) - np.min(FOV_1650), np.max(FOV_1650) - np.mean(FOV_1650)
graph_3 		= ax1.errorbar(1650*0.0003, np.mean(FOV_1650), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1650, linewidth = 0.0, elinewidth = 2.0, capsize=4)
graph_3 		= ax2.errorbar(1650*0.0003, np.mean(FOV_1650), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1650, linewidth = 0.0, elinewidth = 2.0, capsize=4)


ax1.set_ylabel('Freshwater transport (Sv)')
ax1.set_xlim(0.44, 0.505)
ax1.set_ylim(-0.25, 0)
ax1.set_yticks([-0.2, -0.1, 0])
ax1.grid()

ax2.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax2.set_ylabel('Freshwater transport (Sv)')
ax2.set_xlim(0.44, 0.505)
ax2.set_ylim(-0.25, 0)
ax2.set_yticks([-0.2, -0.1, 0])
ax2.grid()

for real_i in range(10):
	#Plot a few random realisations
	ax1.plot(F_H, FOV_mean_1[real_i], '-', color = 'gray', alpha = 0.75, linewidth = 0.5)
	ax2.plot(F_H, FOV_mean_2[real_i], '-', color = 'firebrick', alpha = 0.5, linewidth = 0.5)

#Plot the means	
ax1.plot(F_H, np.mean(FOV_mean_1, axis = 0), '-k', linewidth = 2.0)
ax2.plot(F_H, np.mean(FOV_mean_2, axis = 0), '-r', linewidth = 2.0)

#For the legend
ax2.plot(F_H, np.mean(FOV_mean_1, axis = 0)-100, '-k', linewidth = 2.0, label = 'Not-a-knot')
ax2.plot(F_H, np.mean(FOV_mean_2, axis = 0)-100, '-r', linewidth = 2.0, label = 'Natural')

legend_1	= ax2.legend(loc=(0.015, 0.90), ncol=1, framealpha = 1.0, numpoints = 1)

ax1.set_title('a) $F_{\mathrm{ovS}}$ steady states, cubic splines fit')

#-----------------------------------------------------------------------------------------	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.bar(F_H_knot-0.0035, FOV_min_1_knot / realisation, width = 0.0028, facecolor='white', edgecolor = 'k',linewidth = 2, zorder =5)
ax.bar(F_H_knot, FOV_min_2_knot / realisation, width = 0.0028, facecolor='white', edgecolor = 'firebrick', linewidth = 2, zorder =5)
ax.bar(F_H_knot+0.0035, FOV_min_3_knot / realisation, width = 0.0028, facecolor='white', edgecolor = 'darkorchid', linewidth = 2, zorder = 5)

ax.bar(F_H_knot-0.0035, FOV_min_1_knot / realisation, width = 0.0028, facecolor='black', linewidth = 2, alpha = 0.30, zorder =5, label = 'Not-a-knot')
ax.bar(F_H_knot, FOV_min_2_knot / realisation, width = 0.0028, facecolor='red', linewidth = 2, alpha = 0.30, zorder = 5, label = 'Natural')
ax.bar(F_H_knot+0.0035, FOV_min_3_knot / realisation, width = 0.0028, facecolor='darkorchid', linewidth = 2, alpha = 0.30, zorder =5, label = 'Steady state')

ax.plot(F_H, FOV_min_1 / realisation, '-k', zorder = 10)
ax.plot(F_H, FOV_min_2 / realisation, '-r', zorder = 10)

ax.plot([0.4425, 0.4425], [-0.1, 0.5], '--', color = 'gray')
ax.plot([0.4575, 0.4575], [-0.1, 0.5], '--', color = 'gray')
ax.plot([0.4725, 0.4725], [-0.1, 0.5], '--', color = 'gray')
ax.plot([0.4875, 0.4875], [-0.1, 0.5], '--', color = 'gray')
ax.plot([0.5025, 0.5025], [-0.1, 0.5], '--', color = 'gray')

ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Frequency')
ax.set_xlim(0.44, 0.505)
ax.set_ylim(0, 0.45)
ax.grid()

legend_1	= ax.legend(loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) Probability distribution function of the $F_{\mathrm{ovS}}$ minimum')

show()


