#The continuation for the E-CCM (temperature only) and for varying hosing location

from pylab import *
import numpy as np
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory 	= '../../Data/E-CCM/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Continuation/Hosing_location/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_*.nc')
files.sort()

p_hosing    = np.zeros(len(files))
tipping_1   = np.zeros(len(files))
tipping_2   = np.zeros(len(files))
min_diff    = np.zeros(len(files))

for file_i in range(len(files)):
    
	#Get the p hosing value
	p_hosing[file_i]   = float(files[file_i][-6:-3]) / 100.

	fh = netcdf.Dataset(files[file_i], 'r')

	#Writing data to correct variable	
	E_A	    = fh.variables['E_A'][:]     	  	
	q_N	    = fh.variables['q_n'][:]     	  	
	F_ov    = fh.variables['F_ovS'][:]   
	T_t     = fh.variables['T_t'][:]     	  
	T_n     = fh.variables['T_n'][:]     	  
	T_ts    = fh.variables['T_ts'][:]     	  
	T_s     = fh.variables['T_s'][:]     	

	fh.close()

	E_A_1, q_N_1, F_ov_1	= E_A[0], q_N[0], F_ov[0]
	E_A_2, q_N_2, F_ov_2	= E_A[1], q_N[1], F_ov[1]
	E_A_3, q_N_3, F_ov_3	= E_A[2], q_N[2], F_ov[2]

	tipping_1[file_i]       = np.max(E_A_1)
	tipping_2[file_i]       = np.min(E_A_2)
	
	if p_hosing[file_i] < 0.84:
		min_diff[file_i]       	= np.max(E_A_1) - E_A_1[np.argmin(F_ov_1)]
		
	else:
		#The minimum is so close between two points (accuracy set-up), use the next one to obtain smooth curve
		min_diff[file_i]       	= np.max(E_A_1) - E_A_1[np.argmin(F_ov_1)+1]
	
#-----------------------------------------------------------------------------------------
    
fig, ax = subplots()

graph_1 = ax.plot(1.0-p_hosing, tipping_1, '-', color = 'darkorchid', linewidth = 2.0, label = '$E_{A}^1$ (AMOC on)')
graph_2  = ax.plot(1.0-p_hosing, tipping_2, '-', color = 'firebrick', linewidth = 2.0, label = '$E_{A}^2$ (AMOC off)')

ax.set_xlim(0, 1)
ax.set_ylim(0, 1.0)
ax.grid()

ax.set_xlabel(r'Hosing distribution, $\xi$')
ax.set_ylabel('Freshwater flux forcing, $E_A$ (Sv)')

ax2 	= ax.twinx()

graph_3 = ax2.plot(1.0-p_hosing, (tipping_1 - tipping_2)*10, '--', color = 'royalblue', linewidth = 2.0, label = r'$E_{A}^1 - E_{A}^2$ $(\times 10)$')
graph_4 = ax2.plot(1.0-p_hosing, min_diff * 100, '--', color = 'k', linewidth = 2.0, label = '$E_{A}^1 - E_A$ at $F_{\mathrm{ovS}}$ min.')

ax2.set_ylim(0, 5)
ax2.set_yticks([0, 1, 2, 3, 4, 5])
ax2.set_ylabel(r'Freshwater flux forcing difference ($\times 10^{-2}$ Sv)')

graphs	      	= graph_1 + graph_2
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1, fontsize = 9)

graphs	      	= graph_3 + graph_4
legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax2.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1, fontsize = 9)

ax.add_artist(legend_1)

ax.set_title('d) Saddle-node bifurcations and difference to $F_{\mathrm{ovS}}$ minimum')

show()