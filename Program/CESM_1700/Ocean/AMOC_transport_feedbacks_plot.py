#Program plots the different AMOC feedbacks

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/CESM_1700/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	


fh = netcdf.Dataset(directory+'Ocean/AMOC_feedbacks.nc', 'r')

time		= fh.variables['time'][:]		
gyre		= fh.variables['Gyre_F'][:]	#Gyre responses (Sv)
residual	= fh.variables['Resi_F'][:]	#Residual responses (Sv)
FovN		= fh.variables['FovN_F'][:]	#FovN responses (Sv)
SAF		    = fh.variables['SAF_F'][:]	#Salt-advection feedback responses (Sv)
surf		= fh.variables['Surf_F'][:]	#Surface responses (Sv)
melt		= fh.variables['Melt_F'][:]	#Melt responses (Sv)
prec_evap	= fh.variables['PmE_F'][:]	#Precipitation minus Evaporation responses (Sv)

fh.close()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

graph_1	= ax.plot(time, gyre, '-', color = 'r', linewidth = 2, label = 'Gyres')
graph_2	= ax.plot(time, FovN, '-', color = 'b', linewidth = 2, label = 'Overturning at 65$^{\circ}$N ($F_{\mathrm{ovN}}$)')
graph_3	= ax.plot(time, SAF, '-', color = 'k', linewidth = 2, label = 'Salt-advection feedback ($\sim$$F_{\mathrm{ovS}}$)')
graph_4	= ax.plot(time, surf, '-', color = 'firebrick', linewidth = 2, label = 'Surface')


ax.set_xlabel('Model year after branching')
ax.set_ylabel('Volume transport response (Sv)')
ax.set_xlim(0, 400)
ax.set_ylim(-40, 40)
ax.grid()

ax2 	= fig.add_axes([0.20, 0.59, 0.37, 0.20])

graph_5	= ax2.plot(time, melt, '-', color = 'darkturquoise', linewidth = 2, label = 'Melt')
graph_6	= ax2.plot(time, prec_evap, '-', color = 'gray', linewidth = 2, label = 'P-E')
ax2.set_xlim(0, 400)
ax2.set_ylim(-30, 5)
ax2.grid()
ax2.set_title('Surface decomposition')

graphs	      = graph_1 + graph_2 + graph_3 + graph_4

legend_labels = [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

graphs	      	= graph_6 + graph_5

legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax2.legend(graphs, legend_labels, loc=(1.02, 0.30), ncol=1, framealpha = 1.0, numpoints = 1)


ax.set_title('b) AMOC response decomposition, $\overline{F_H} = 0.51$ Sv')

show()