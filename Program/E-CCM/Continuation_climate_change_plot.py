#The continuation for the E-CCM (temperature only) and under climate change

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

p_hosing	= 1.0

fh = netcdf.Dataset(directory+'Continuation/Climate_change/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'_T_t_0_0_T_n_0_0.nc', 'r')

#Writing data to correct variable	
E_A		= fh.variables['E_A'][:]     	  	
q_N		= fh.variables['q_n'][:]  

fh.close()

E_A_1, q_N_1	= E_A[0], q_N[0]
E_A_2, q_N_2	= E_A[1], q_N[1]
E_A_3, q_N_3	= E_A[2], q_N[2]

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Continuation/Climate_change/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'_T_t_0_0_T_n_5_0.nc', 'r')

#Writing data to correct variable	
E_A		= fh.variables['E_A'][:]     	  	
q_N		= fh.variables['q_n'][:]  
F_ovS		= fh.variables['F_ovS'][:]          	  	        	  	

fh.close()

E_A_1_CC, q_N_1_CC, F_ovS_1_CC	= E_A[0], q_N[0], F_ovS[0]
E_A_2_CC, q_N_2_CC, F_ovS_2_CC	= E_A[1], q_N[1], F_ovS[1]
E_A_3_CC, q_N_3_CC, F_ovS_3_CC	= E_A[2], q_N[2], F_ovS[2]

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(E_A_3, q_N_3, ':', color = 'gray', linewidth = 1.5)
ax.plot(E_A_1, q_N_1, '-', color = 'darkviolet', linewidth = 2.0)
ax.plot(E_A_2, q_N_2, '-', color = 'firebrick', linewidth = 2.0)

ax.plot(E_A_3_CC, q_N_3_CC, linestyle= 'dashdot', color = 'gray', linewidth = 1.5)
ax.plot(E_A_1_CC, q_N_1_CC, '--', color = 'darkviolet', linewidth = 2.0, alpha = 0.9)
ax.plot(E_A_2_CC, q_N_2_CC, '--', color = 'firebrick', linewidth = 2.0, alpha = 0.9)

ax.set_xlim(0, 0.60)
ax.set_ylim(-2, 18)
ax.set_xlabel('Freshwater flux forcing, $E_A$ (Sv)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_yticks([0, 5, 10, 15])
ax.grid()


graph_tipping_1 = ax.plot(E_A_1, q_N_1-100, '-k', linewidth = 2.0, label = '$\Delta T_{\mathrm{n}}^a = 0^{\circ}$C')
graph_tipping_2 = ax.plot(E_A_1, q_N_1-100, '--k', linewidth = 2.0, label = '$\Delta T_{\mathrm{n}}^a = 5^{\circ}$C')

graphs_1	= graph_tipping_1 + graph_tipping_2
legend_labels_1 = [l.get_label() for l in graphs_1]
legend_1	= ax.legend(graphs_1, legend_labels_1, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

graph_AMOC_on  = ax.plot(E_A_1, q_N_1-100, '-', color = 'darkviolet', linewidth = 2.0, label = 'AMOC on')
graph_AMOC_un  = ax.plot(E_A_1, q_N_1-100, ':', color = 'gray', linewidth = 2.0, label = 'Unstable')
graph_AMOC_off = ax.plot(E_A_1, q_N_1-100, '-', color = 'firebrick', linewidth = 2.0, label = 'AMOC off')

graphs_2	= graph_AMOC_on + graph_AMOC_un + graph_AMOC_off

legend_labels_2 = [l.get_label() for l in graphs_2]
legend_2	= ax.legend(graphs_2, legend_labels_2, loc = 'center left', ncol=1, framealpha = 1.0)
ax.add_artist(legend_1)

ax.set_title('a) AMOC strength steady states')

show()

