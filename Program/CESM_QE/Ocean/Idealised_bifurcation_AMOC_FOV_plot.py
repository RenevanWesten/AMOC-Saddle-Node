#Program plots the idealised AMOC responses from the reduced model

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/'
 
def ReadinDataFOV(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	    #Freshwater (Sv)

	fh.close()

	return time, FOV

def ReadinDataAMOC(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	AMOC		= fh.variables['Transport'][:]	#AMOC strength

	fh.close()

	return time, AMOC

def AMOCvsHosing(Psi, F_ovS, F_H, c_1 = 0.52, c_2 = 20, n_1 = 0.025, g_1 = 0.032, S_0 = 35):
	"""The initial AMOC and FovS values, against the surface forcing"""
	
	#Repeated constant of salinity and temperaure AMOC responses
	con_salt_temp	= c_2 * (1 - c_1) * S_0	
	
	#The base part for both solutions
	AMOC_base	= (Psi / 2.0) - (con_salt_temp * F_ovS / (2 * Psi)) - ((n_1 + g_1) * S_0 / 2)
	
	#The square root part
	AMOC_sqrt	= ((Psi**2.0 + con_salt_temp * F_ovS + (n_1 + g_1) * S_0 * Psi) / (2 * Psi))**2.0 - con_salt_temp * F_H
	
	AMOC_on 	= AMOC_base + np.sqrt(AMOC_sqrt)
	AMOC_un 	= AMOC_base - np.sqrt(AMOC_sqrt)
	
	#Remove the masked elements
	AMOC_on		= ma.masked_where(AMOC_on != AMOC_on, AMOC_on)
	AMOC_un		= ma.masked_where(AMOC_un != AMOC_un, AMOC_un)
	
	#Determine the initial vertical salinity difference
	salt_vert_0	= - S_0 * F_ovS / Psi
	
	#Determine the vertical salinity difference against the hosing
	salt_vert_on	= salt_vert_0 + (Psi - AMOC_on) / (c_2 * (1 - c_1))
	salt_vert_un	= salt_vert_0 + (Psi - AMOC_un) / (c_2 * (1 - c_1))
	
	#Determine the FovS
	F_ovS_on	= - salt_vert_on * AMOC_on / S_0
	F_ovS_un	= - salt_vert_un * AMOC_un / S_0
	
	return AMOC_on, AMOC_un, F_ovS_on, F_ovS_un
	
def ForcingMax(Psi, F_ovS, c_1 = 0.52, c_2 = 20, n_1 = 0.025, g_1 = 0.032, S_0 = 35):
	"""The initial AMOC and FovS values, to find the maximum surface forcing"""

	#Repeated constant of salinity and temperaure AMOC responses
	con_salt_temp	= c_2 * (1 - c_1) * S_0	
	
	F_H_max		= 1 / con_salt_temp * ((Psi**2.0 + con_salt_temp * F_ovS + (n_1 + g_1) * S_0 * Psi) / (2 * Psi))**2.0
	
	return F_H_max
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000


time, AMOC	= ReadinDataAMOC(directory+'CESM_QE/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time, FOV	= ReadinDataFOV(directory+'CESM_QE/Ocean/FOV_section_34S.nc')

#Get the steady states
time_0600, AMOC_0600	= ReadinDataAMOC(directory+'CESM_0600/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')	
time_1500, AMOC_1500	= ReadinDataAMOC(directory+'CESM_1500/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')
time_1550, AMOC_1550	= ReadinDataAMOC(directory+'CESM_1550/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_QE.nc')	
time_1600, AMOC_1600	= ReadinDataAMOC(directory+'CESM_1600/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_1550.nc')	
time_1650, AMOC_1650	= ReadinDataAMOC(directory+'CESM_1650/Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_branch_1600.nc')	

time_0600, FOV_0600	    = ReadinDataFOV(directory+'CESM_0600/Ocean/FOV_section_34S_branch_QE.nc')	
time_1500, FOV_1500	    = ReadinDataFOV(directory+'CESM_1500/Ocean/FOV_section_34S_branch_QE.nc')
time_1550, FOV_1550	    = ReadinDataFOV(directory+'CESM_1550/Ocean/FOV_section_34S_branch_QE.nc')	
time_1600, FOV_1600	    = ReadinDataFOV(directory+'CESM_1600/Ocean/FOV_section_34S_branch_1550.nc')	
time_1650, FOV_1650	    = ReadinDataFOV(directory+'CESM_1650/Ocean/FOV_section_34S_branch_1600.nc')	

#Now only get the last 50 years
AMOC_0600, FOV_0600	= AMOC_0600[-50:], FOV_0600[-50:]
AMOC_1500, FOV_1500	= AMOC_1500[-50:], FOV_1500[-50:]
AMOC_1550, FOV_1550	= AMOC_1550[-50:], FOV_1550[-50:]
AMOC_1600, FOV_1600	= AMOC_1600[-50:], FOV_1600[-50:]
AMOC_1650, FOV_1650	= AMOC_1650[-50:], FOV_1650[-50:]

F_H_steady_states		= np.asarray([0.000, 0.180, 0.450, 0.465, 0.480, 0.495])
AMOC_steady_states		= np.asarray([np.mean(AMOC[:50]), np.mean(AMOC_0600), np.mean(AMOC_1500), np.mean(AMOC_1550), np.mean(AMOC_1600), np.mean(AMOC_1650)])


#-----------------------------------------------------------------------------------------

print('Mean (1-50): ', np.mean(AMOC[:50]), np.mean(FOV[:50]))
print('Max (1-50) : ', np.max(AMOC[:50]), np.max(FOV[:50]))
print('Min (1-50) : ', np.min(AMOC[:50]), np.min(FOV[:50]))

#-----------------------------------------------------------------------------------------

AMOC_0		= np.arange(5, 25.1, 0.1)
FovS_0		= np.arange(-0.3, 0.301, 0.01)

x, y		= np.meshgrid(AMOC_0, FovS_0)

F_H_max		= ForcingMax(x, y)

fig, ax	= subplots()

ax.fill_between([0, 25], y1 = np.zeros(2) + -0.5, y2 = np.zeros(2) + 0.5, color = 'gray', alpha = 0.25)

#The 86% is needed to take the total surface forcing, so +F_H and -F_H
CS	= ax.contourf(x, y, F_H_max / 0.86, levels = np.arange(0, 0.601, 0.025), extend = 'max', cmap = 'Reds')
cbar	= colorbar(CS, ticks = np.arange(0, 0.601, 0.1))
cbar.set_label('Critical freshwater flux forcing, $F_H^{\mathrm{c}}$ (Sv)')

#CESM first 50 years
ax.plot([np.min(AMOC[:50]), np.min(AMOC[:50]), np.max(AMOC[:50]), np.max(AMOC[:50]), np.min(AMOC[:50])], [np.min(FOV[:50]), np.max(FOV[:50]), np.max(FOV[:50]), np.min(FOV[:50]), np.min(FOV[:50])], '-k', linewidth = 2.0)

#Observations
ax.plot([16, 16, 19, 19, 16], [-0.28, -0.05, -0.05, -0.28, -0.28], '-b', linewidth = 2.0)

ax.text(np.mean(AMOC[:50]), np.min(FOV[:50])-0.01, 'CESM', verticalalignment='top', horizontalalignment='center', color = 'k', fontsize = 13)
ax.text(17.5, -0.05+0.01, 'Observations', verticalalignment='bottom', horizontalalignment='center', color = 'b', fontsize = 13)

ax.set_xlim(0, 25)
ax.set_ylim(-0.3, 0.3)
ax.set_xlabel('Initial AMOC strength (Sv)')
ax.set_ylabel('Initial $F_{\mathrm{ovS}}$ (Sv)')
ax.grid()

ax.set_title('c) Critical freshwater flux forcing, $F_H^{\mathrm{c}}$')

#-----------------------------------------------------------------------------------------

AMOC_0		= np.mean(AMOC[:50])
FovS_0		= np.mean(FOV[:50])
F_H		   = np.arange(0, 0.6000001, 0.00001)

AMOC_1, AMOC_un_1, FovS_1, FovS_un_1	= AMOCvsHosing(AMOC_0, FovS_0, F_H * 0.86)
AMOC_2, AMOC_un_2, FovS_2, FovS_un_2 	= AMOCvsHosing(np.max(AMOC[:50]), np.max(FOV[:50]), F_H * 0.86)
AMOC_3, AMOC_un_3, FovS_3, FovS_un_3 	= AMOCvsHosing(np.min(AMOC[:50]), np.min(FOV[:50]), F_H * 0.86)
AMOC_4, AMOC_un_4, FovS_4, FovS_un_4	= AMOCvsHosing(17, -0.15, F_H * 0.86)	#Observations

#-----------------------------------------------------------------------------------------

fig, ax		= subplots()

graph_1		= ax.plot(F_H, AMOC_1, '-', color = 'k', linewidth = 1.5, label = 'CESM (mean)')
graph_2		= ax.plot(F_H, AMOC_2, '-', color = 'r', linewidth = 1.5, label = 'CESM (max)')
graph_3		= ax.plot(F_H, AMOC_3, '-', color = 'b', linewidth = 1.5, label = 'CESM (min)')
graph_4		= ax.plot(F_H, AMOC_4, '-', color = 'firebrick', linewidth = 1.5, label = 'Observations')

ax.plot(F_H, AMOC_un_1, ':', color = 'k', linewidth = 1.5)
ax.plot(F_H, AMOC_un_2, ':', color = 'r', linewidth = 1.5)
ax.plot(F_H, AMOC_un_3, ':', color = 'b', linewidth = 1.5)
ax.plot(F_H, AMOC_un_4, ':', color = 'firebrick', linewidth = 1.5)

y_error_0600		= np.zeros((2, 1)) + 1
y_error_0600[:, 0]	= np.mean(AMOC_0600) - np.min(AMOC_0600), np.max(AMOC_0600) - np.mean(AMOC_0600)
ax.errorbar(600*0.0003, np.mean(AMOC_0600), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_0600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(AMOC_1500) - np.min(AMOC_1500), np.max(AMOC_1500) - np.mean(AMOC_1500)
graph_3 		= ax.errorbar(1500*0.0003, np.mean(AMOC_1500), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')

y_error_1550		= np.zeros((2, 1)) + 1
y_error_1550[:, 0]	= np.mean(AMOC_1550) - np.min(AMOC_1550), np.max(AMOC_1550) - np.mean(AMOC_1550)
graph_3 		= ax.errorbar(1550*0.0003, np.mean(AMOC_1550), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1550, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1600		= np.zeros((2, 1)) + 1
y_error_1600[:, 0]	= np.mean(AMOC_1600) - np.min(AMOC_1600), np.max(AMOC_1600) - np.mean(AMOC_1600)
graph_3 		= ax.errorbar(1600*0.0003, np.mean(AMOC_1600), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1650		= np.zeros((2, 1)) + 1
y_error_1650[:, 0]	= np.mean(AMOC_1650) - np.min(AMOC_1650), np.max(AMOC_1650) - np.mean(AMOC_1650)
graph_3 		= ax.errorbar(1650*0.0003, np.mean(AMOC_1650), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1650, linewidth = 0.0, elinewidth = 2.0, capsize=4)


ax.set_xlim(0, 0.6)
ax.set_ylim(0, 20)
ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('AMOC strength (Sv)')
ax.grid()
ax.set_yticks([0, 5, 10, 15, 20])

legend_1	= ax.legend(loc=(0.82, 0.69), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) Idealised AMOC responses under hosing')

#-----------------------------------------------------------------------------------------

fig, ax		= subplots()

graph_1		= ax.plot(F_H, FovS_1, '-', color = 'k', linewidth = 1.5, label = 'CESM (mean)')
graph_2		= ax.plot(F_H, FovS_2, '-', color = 'r', linewidth = 1.5, label = 'CESM (max)')
graph_3		= ax.plot(F_H, FovS_3, '-', color = 'b', linewidth = 1.5, label = 'CESM (min)')
graph_4		= ax.plot(F_H, FovS_4, '-', color = 'firebrick', linewidth = 1.5, label = 'Observations')

ax.plot(F_H, FovS_un_1, ':', color = 'k', linewidth = 1.5)
ax.plot(F_H, FovS_un_2, ':', color = 'r', linewidth = 1.5)
ax.plot(F_H, FovS_un_3, ':', color = 'b', linewidth = 1.5)
ax.plot(F_H, FovS_un_4, ':', color = 'firebrick', linewidth = 1.5)

y_error_0600		= np.zeros((2, 1)) + 1
y_error_0600[:, 0]	= np.mean(FOV_0600) - np.min(FOV_0600), np.max(FOV_0600) - np.mean(FOV_0600)
ax.errorbar(600*0.0003, np.mean(FOV_0600), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_0600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1500		= np.zeros((2, 1)) + 1
y_error_1500[:, 0]	= np.mean(FOV_1500) - np.min(FOV_1500), np.max(FOV_1500) - np.mean(FOV_1500)
graph_3 		= ax.errorbar(1500*0.0003, np.mean(FOV_1500), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1500, linewidth = 0.0, elinewidth = 2.0, capsize=4, label = 'Steady state')

y_error_1550		= np.zeros((2, 1)) + 1
y_error_1550[:, 0]	= np.mean(FOV_1550) - np.min(FOV_1550), np.max(FOV_1550) - np.mean(FOV_1550)
graph_3 		= ax.errorbar(1550*0.0003, np.mean(FOV_1550), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1550, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1600		= np.zeros((2, 1)) + 1
y_error_1600[:, 0]	= np.mean(FOV_1600) - np.min(FOV_1600), np.max(FOV_1600) - np.mean(FOV_1600)
graph_3 		= ax.errorbar(1600*0.0003, np.mean(FOV_1600), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1600, linewidth = 0.0, elinewidth = 2.0, capsize=4)

y_error_1650		= np.zeros((2, 1)) + 1
y_error_1650[:, 0]	= np.mean(FOV_1650) - np.min(FOV_1650), np.max(FOV_1650) - np.mean(FOV_1650)
graph_3 		= ax.errorbar(1650*0.0003, np.mean(FOV_1650), color = 'darkorchid', marker = 's', markerfacecolor = 'darkorchid', markeredgecolor = 'mediumpurple', yerr = y_error_1650, linewidth = 0.0, elinewidth = 2.0, capsize=4)

ax.set_xlim(0, 0.6)
ax.set_ylim(-0.35, 0.35)
ax.set_xlabel('Freshwater flux forcing, $F_H$ (Sv)')
ax.set_ylabel('Freshwater transport (Sv)')
ax.grid()

legend_1	= ax.legend(loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) Idealised $F_{\mathrm{ovS}}$ responses under hosing')

#-----------------------------------------------------------------------------------------

AMOC_0		= np.mean(AMOC[:50])
FovS_0		= np.mean(FOV[:50])
F_H		= np.arange(0, 0.6000001, 0.00001)

gyre_coef		= np.arange(0, 0.101, 0.001)
F_H_max_gyre_1		= ma.masked_all(len(gyre_coef))
F_H_max_gyre_2		= ma.masked_all(len(gyre_coef))
F_H_min_diff_gyre_1	= ma.masked_all(len(gyre_coef))
F_H_min_diff_gyre_2	= ma.masked_all(len(gyre_coef))

for gyre_i in range(len(gyre_coef)):
	#Determine the idealised bifurcation diagram
	AMOC_1, AMOC_un_1, FovS_1, FovS_un_1	= AMOCvsHosing(AMOC_0, FovS_0, F_H * 0.86, n_1 = 0, g_1 = gyre_coef[gyre_i])
	F_H_max_gyre_1[gyre_i]			= ForcingMax(AMOC_0, FovS_0, n_1 = 0, g_1 = gyre_coef[gyre_i]) / 0.86
	F_H_min					= F_H[np.argmin(FovS_1)]
	F_H_min_diff_gyre_1[gyre_i] 		= F_H_max_gyre_1[gyre_i] - F_H_min

	AMOC_2, AMOC_un_2, FovS_2, FovS_un_2	= AMOCvsHosing(AMOC_0, FovS_0, F_H * 0.86, g_1 = gyre_coef[gyre_i])
	F_H_max_gyre_2[gyre_i]			= ForcingMax(AMOC_0, FovS_0, g_1 = gyre_coef[gyre_i]) / 0.86
	F_H_min					= F_H[np.argmin(FovS_2)]
	F_H_min_diff_gyre_2[gyre_i] 		= F_H_max_gyre_2[gyre_i] - F_H_min
	
fig, ax		= subplots()

ax.plot([0.032, 0.032], [-5, 5], ':', linewidth = 2.0, color = 'royalblue')
graph_1		= ax.plot(gyre_coef, F_H_max_gyre_2, '-', color = 'k', linewidth = 2, label = '$n_1$ = 0.025 Sv kg$^{-1}$ g$^{-1}$')
graph_2		= ax.plot(gyre_coef, F_H_max_gyre_1, '-', color = 'r', linewidth = 2, label = '$n_1$ = 0 Sv kg$^{-1}$ g$^{-1}$')

legend_1	= ax.legend(loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_xlim(0, 0.1)
ax.set_ylim(0.30, 0.60)
ax.set_xlabel('Gyre sensitivity, $g_1$ (Sv kg$^{-1}$ g$^{-1}$)')
ax.set_ylabel('Critical freshwater flux forcing, $F_H^{\mathrm{c}}$ (Sv)')

ax.grid()

ax2 	= ax.twinx()

graph_1	= ax2.plot(gyre_coef, F_H_min_diff_gyre_2 * 10**2.0, '--', color = 'k', linewidth = 2, label = '$n_1$ = 0.025 Sv kg$^{-1}$ g$^{-1}$')
graph_2	= ax2.plot(gyre_coef, F_H_min_diff_gyre_1 * 10**2.0, '--', color = 'r', linewidth = 2, label = '$n_1$ = 0 Sv kg$^{-1}$ g$^{-1}$')



graph_3  	= ax.plot(0, -100, '-k', linewidth = 2.0, label = '$F_H^{\mathrm{c}}$')
graph_4 	= ax.plot(0, -100, '--k', linewidth = 2.0, label = '$F_H^{\mathrm{c}} - F_H$ at $F_{\mathrm{ovS}}$ min. ')

graphs_2	= graph_3 + graph_4

legend_labels_2 = [l.get_label() for l in graphs_2]
legend_2	= ax.legend(graphs_2, legend_labels_2, loc = 'lower right', ncol=1, framealpha = 1.0)
ax.add_artist(legend_1)


ax2.set_ylim(-0.1, 2)
ax2.set_yticks([0, 0.5, 1, 1.5, 2])
ax2.set_ylabel(r'Freshwater flux forcing difference ($\times 10^{-2}$ Sv)')

ax.set_title('d) $F_H^{\mathrm{c}}$ and difference to $F_{\mathrm{ovS}}$ minimum')

show()






