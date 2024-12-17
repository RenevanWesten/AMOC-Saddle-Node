#The E-CCM (temperature only) and for varying hosing location

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

def Heaviside(x):
	"""Returns 1 for positive values, 0 otherwise"""
	if x > 0:
		return 1.0
	return 0

	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#Fixed parameters
V_0	    = 3.0 * 10**17.0 	#Total volume of basin
V_n	    = 3.0 * 10**15.0	#Volume of northern box
V_s	    = 9.0 * 10**15.0	#Volume of southern box
A	    = 1.0 * 10**14.0	#Area of Atlantic pycnocline
L_xA	= 1.0 * 10**7.0		#Zonal extent of Atlantic Ocean at southern end
L_y	    = 1.0 * 10**6.0		#Meridional extent of frontal region
L_xS	= 3.0 * 10**7.0		#Zonal extent of Southern Ocean
tau	    = 0.1			#Wind stress
A_GM	= 1700.0		#Eddy diffusivity
f_s	    = -10**(-4.0)		#Coriolis parameter
rho_0	= 1027.5		#Reference density
kappa	= 10**(-5.0)		#Vertical diffusivity
S_0	    = 35.0			#Reference salinity
T_0	    = 5.0			#Reference temperature
eta	    = 3.0 * 10**4.0		#Hydraulic constant
alpha	= 2.0 * 10**(-4.0)	#Thermal expansion coefficient
beta	= 8.0 * 10**(-4.0)	#Haline contraction coefficient
r_S	    = 1.0 * 10**7.0		#Transport by the southern subtropical gyre
r_N	    = 5.0 * 10**6.0		#Transport by the northern subtropical gyre
E_S	    = 0.17 * 10**6.0	#Symmetric freshwater flux
#-----------------------------------------------------------------------------------------

p_hosing		    = 1.0             #xi = 1 - p
lambda_a		    = 0.0000035		  #Atmospheric exchange coefficient (per time)
filename_parameter	= 'MOC_box_model_lambda_a_0_0000035.txt'
filename_output		= 'MOC_box_model_temp_equilibrium_states_lambda_a_0_0000035.nc'

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
T_ts_a, T_n_a, S_t, S_n, S_ts, S_s, T_t, T_n, T_ts, T_s, T_d, D = np.loadtxt(directory+'Initialisation/'+filename_parameter)

#Assume that the tropical and southern box have the same temperatures as the ts-box and n-box, respectively
T_t_a	= np.copy(T_ts_a)
T_s_a	= np.copy(T_n_a)

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
fh = netcdf.Dataset(directory+'Continuation/Hosing_location/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'.nc', 'r')

S_t = fh.variables['S_t'][0,0]     	  	
S_n = fh.variables['S_n'][0,0]     	  	
S_ts= fh.variables['S_ts'][0,0]   	  	
S_s = fh.variables['S_s'][0,0]     	  
S_d = fh.variables['S_d'][0,0]     	 
T_t = fh.variables['T_t'][0,0]     	  
T_n = fh.variables['T_n'][0,0]     	  
T_ts= fh.variables['T_ts'][0,0]     	  
T_s = fh.variables['T_s'][0,0]     	
T_d = fh.variables['T_d'][0,0]     	  	  	  
D   = fh.variables['D'][0,0]     	  	  	  

fh.close()

#Constant depending on the standard parameters above
rho_n	= rho_0 * (1.0 - alpha * (T_n - T_0) + beta * (S_n - S_0))
rho_ts	= rho_0 * (1.0 - alpha * (T_ts - T_0) + beta * (S_ts - S_0))
q_Ek	= (tau * L_xS) / (rho_0 * np.fabs(f_s))
q_e	= A_GM * L_xA * D / L_y
q_U	= kappa * A / D
q_N	= (eta * (rho_n - rho_ts) * D**2.0 / rho_0)
q_S	= q_Ek - q_e

V_t	= A * D
V_ts	= L_xA * L_y * D / 2.0
V_d	= V_0 - V_n - V_s - V_t - V_ts

#Salinity conservation
VS_0	= V_0 * S_0
VS_t	= V_t * S_t
VS_ts	= V_ts * S_ts
VS_n	= V_n * S_n
VS_s	= V_s * S_s
VS_d	= VS_0 - VS_n - VS_t - VS_ts - VS_s
S_d	= VS_d / V_d

#Ocean Heat content (assume constant density)
OHC_t	= V_t * (T_t+273.15)
OHC_ts	= V_ts * (T_ts+273.15)
OHC_n	= V_n * (T_n+273.15)
OHC_s	= V_s * (T_s+273.15)
OHC_d	= V_d * (T_d+273.15)

#MOV-index
M_ov	= -q_S / S_0 * (S_ts - S_d)

#Horizontal areas for the exchange with the surface
A_t	= np.copy(A)
A_ts	= L_y * L_xA
A_s	= L_y * L_xS
A_n	= np.copy(A_ts)

#-----------------------------------------------------------------------------------------

delta_t			= 86400 * 5

if p_hosing	== 1.0:
	E_A_max		= 0.474
	
if p_hosing	== 0.5:
	E_A_max		= 0.609
	
if p_hosing	== 0.0:
	E_A_max		= 0.8348

num_year		= int(E_A_max / 0.0003)

if p_hosing == 1.0:
    time			= np.arange((num_year+1000) * 365 * 86400 / delta_t) / (365 * 86400 / delta_t)

if p_hosing == 0.5:
    time			= np.arange((num_year+3000) * 365 * 86400 / delta_t) / (365 * 86400 / delta_t)
    
if p_hosing == 0.0:
    time			= np.arange((num_year+10000) * 365 * 86400 / delta_t) / (365 * 86400 / delta_t)
    
#Empty arrays for each time step
S_t_all		= np.ma.masked_all(len(time))
S_n_all		= np.ma.masked_all(len(time))
S_ts_all	= np.ma.masked_all(len(time))	#Salinity
S_s_all		= np.ma.masked_all(len(time))
S_d_all		= np.ma.masked_all(len(time))

T_t_all		= np.ma.masked_all(len(time))
T_n_all		= np.ma.masked_all(len(time))
T_ts_all	= np.ma.masked_all(len(time))	#Temperature
T_s_all		= np.ma.masked_all(len(time))
T_d_all		= np.ma.masked_all(len(time))

D_all		= np.ma.masked_all(len(time))
q_N_all		= np.ma.masked_all(len(time))	#Remaining quantities
q_S_all		= np.ma.masked_all(len(time))
M_ov_all	= np.ma.masked_all(len(time))
ice_n_all	= np.ma.masked_all(len(time))

E_A_all						                = np.zeros(len(time))
E_A_all[:int(num_year*365*86400/delta_t)]	= np.linspace(0, E_A_max, int(num_year*365*86400/delta_t), endpoint = True)
E_A_all[int(num_year*365*86400/delta_t):]	= E_A_max

E_A_all		= E_A_all * 10**6.0

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

for time_i in range(len(time)):
	#Add stochastic noise
	E_A	= E_A_all[time_i]

	#T-box
	VS_t_delta 	= q_S * (Heaviside(q_S) * S_ts + Heaviside(-q_S) * S_t) + q_U * S_d - Heaviside(q_N) * q_N * S_t + r_S * (S_ts - S_t) + r_N * (S_n - S_t) + (2.0 * E_S - (1.0 - p_hosing) * E_A) * S_0
	OHC_t_delta 	= q_S * (Heaviside(q_S) * (T_ts+273.15) + Heaviside(-q_S) * (T_t+273.15)) + q_U * (T_d+273.15) - Heaviside(q_N) * q_N * (T_t+273.15) + r_S * (T_ts - T_t) + r_N * (T_n - T_t) - A_t * lambda_a * (T_t - T_t_a)
	VS_t		= VS_t + VS_t_delta * delta_t
	OHC_t		= OHC_t + OHC_t_delta * delta_t

	#TS-box
	VS_ts_delta	= q_Ek * S_s - q_e * S_ts - q_S * (Heaviside(q_S) * S_ts + Heaviside(-q_S) * S_t) + r_S * (S_t - S_ts)
	OHC_ts_delta	= q_Ek * (T_s+273.15) - q_e * (T_ts+273.15) - q_S * (Heaviside(q_S) * (T_ts+273.15) + Heaviside(-q_S) * (T_t+273.15)) + r_S * (T_t - T_ts) - A_ts * lambda_a * (T_ts - T_ts_a)
	VS_ts		= VS_ts + VS_ts_delta * delta_t
	OHC_ts		= OHC_ts + OHC_ts_delta * delta_t

	#N-box
	VS_n_delta	= Heaviside(q_N) * q_N * (S_t - S_n) + r_N * (S_t - S_n) - (E_S + p_hosing * E_A) * S_0
	OHC_n_delta	= Heaviside(q_N) * q_N * (T_t - T_n) + r_N * (T_t - T_n) - A_n * lambda_a * (T_n - T_n_a)
	VS_n		= VS_n + VS_n_delta * delta_t
	OHC_n		= OHC_n + OHC_n_delta * delta_t

	#S-box
	VS_s_delta	= q_S * (Heaviside(q_S) * S_d + Heaviside(-q_S) * S_s) + q_e * S_ts - q_Ek * S_s - (E_S - E_A) * S_0
	OHC_s_delta	= q_S * (Heaviside(q_S) * (T_d+273.15) + Heaviside(-q_S) * (T_s+273.15)) + q_e * (T_ts+273.15) - q_Ek * (T_s+273.15) - A_s * lambda_a * (T_s - T_s_a)
	VS_s		= VS_s + VS_s_delta * delta_t
	OHC_s		= OHC_s + OHC_s_delta * delta_t

	#D-box (only for temperature)
	OHC_d_delta	= Heaviside(q_N) * q_N * (T_n+273.15) - q_S * (Heaviside(q_S) * (T_d+273.15) + Heaviside(-q_S) * (T_s+273.15)) - q_U * (T_d+273.15)
	OHC_d		= OHC_d + OHC_d_delta * delta_t

	#Depth pycnocline
	D_delta		= q_U + q_Ek - q_e - Heaviside(q_N) * q_N
	D_delta		= D_delta / (A + 0.5 * L_xA * L_y)
	D		= D + D_delta * delta_t

	#D-box, conservation of salt
	VS_d	= VS_0 - VS_n - VS_t - VS_ts - VS_s

	#Update volume of boxes
	V_t	= A * D
	V_ts	= L_xA * L_y * D / 2.0
	V_d	= V_0 - V_n - V_s - V_t - V_ts

	#Update the salinity
	S_t	= VS_t / V_t
	S_ts	= VS_ts / V_ts
	S_n	= VS_n / V_n
	S_s	= VS_s / V_s
	S_d	= VS_d / V_d

	#Update the temperature
	T_t	   = OHC_t / V_t - 273.15
	T_ts   	   = OHC_ts / V_ts - 273.15
	T_n	   = OHC_n / V_n - 273.15
	T_s	   = OHC_s / V_s - 273.15
	T_d	   = OHC_d / V_d - 273.15

	#Update the fluxes
	rho_n	= rho_0 * (1.0 - alpha * (T_n - T_0) + beta * (S_n - S_0))
	rho_ts	= rho_0 * (1.0 - alpha * (T_ts - T_0) + beta * (S_ts - S_0))
	q_Ek	= (tau * L_xS) / (rho_0 * np.fabs(f_s))
	q_e	= A_GM * L_xA * D / L_y
	q_U	= kappa * A / D
	q_N	= eta * (rho_n - rho_ts) * D**2.0 / rho_0
	q_S	= q_Ek - q_e

	#MOV-index
	M_ov	= -q_S / S_0 * (S_ts - S_d)

	if time_i % (365 * 50 * 86400.0 / delta_t) == 0:
		print(time_i * (delta_t / 86400.0) / 365.0, E_A / 10**6.0, q_N / 10**6.0)

	#Save data
	S_t_all[time_i]		= S_t
	S_n_all[time_i]		= S_n
	S_ts_all[time_i]	= S_ts	#Salinity
	S_s_all[time_i]		= S_s
	S_d_all[time_i]		= S_d
	T_t_all[time_i]		= T_t
	T_n_all[time_i]		= T_n
	T_ts_all[time_i]	= T_ts	#Temperature
	T_s_all[time_i]		= T_s
	T_d_all[time_i]		= T_d
	D_all[time_i]		= D
	q_N_all[time_i]		= Heaviside(q_N) * q_N / 10**6.0
	q_S_all[time_i]		= q_S / 10**6.0
	M_ov_all[time_i]	= M_ov / 10**6.0	

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Continuation/Hosing_location/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'.nc', 'r')

E_A_equil	    = fh.variables['E_A'][:]     	  	
q_N_equil	    = fh.variables['q_n'][:]     	  	  	
FOV_equil	    = fh.variables['F_ovS'][:]  
S_n_equil	    = fh.variables['S_n'][:]    	  	  	   	  	  	

fh.close()
   

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()


graph_1 = ax.plot(E_A_equil[0], q_N_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
graph_2 = ax.plot(E_A_equil[1], q_N_equil[1], '-', color = 'firebrick', label='AMOC off')
graph_3 = ax.plot(E_A_equil[2], q_N_equil[2], '--', color = 'gray', label='Unstable')

graph_4 = ax.plot(E_A_all / 10**6.0, q_N_all, '-k', linewidth = 2, label = 'Quasi-equilibrium')
ax.plot(E_A_all[int(num_year*365*86400/delta_t)] / 10**6.0, q_N_all[int(num_year*365*86400/delta_t)], 'ok')


ax.set_xlabel('Freshwater flux forcing, $E_A$ (Sv)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(0, 1)
ax.set_ylim(-2, 18)
ax.set_yticks([0, 5, 10, 15])
ax.grid()

graphs	      	= graph_1 + graph_2 + graph_3 + graph_4
legend_labels 	= [l.get_label() for l in graphs]

if p_hosing == 1.0:

	legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)
	
	ax2 	= fig.add_axes([0.175, 0.40, 0.25, 0.20])
	
	graph_1 = ax2.plot(E_A_equil[0], q_N_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
	graph_3 = ax2.plot(E_A_equil[2], q_N_equil[2], '--', color = 'gray', label='Unstable')
	graph_4 = ax2.plot(E_A_all / 10**6.0, q_N_all, '-k', linewidth = 2, label = 'Quasi-equilibrium')
	ax2.plot(E_A_all[int(num_year*365*86400/delta_t)] / 10**6.0, q_N_all[int(num_year*365*86400/delta_t)], 'ok')

	ax2.set_xlim(0.47, 0.49)
	ax2.set_ylim(8, 12)
	ax2.grid()
	ax2.set_title('AMOC strength (Sv)', fontsize = 10)
	
	ax3	= fig.add_axes([0.68, 0.40, 0.25, 0.20])
	
	graph_1 = ax3.plot(E_A_equil[0], S_n_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
	graph_3 = ax3.plot(E_A_equil[2], S_n_equil[2], '--', color = 'gray', label='Unstable')
	graph_4 = ax3.plot(E_A_all / 10**6.0, S_n_all, '-k', linewidth = 2, label = 'Quasi-equilibrium')
	ax3.plot(E_A_all[int(num_year*365*86400/delta_t)] / 10**6.0, S_n_all[int(num_year*365*86400/delta_t)], 'ok')

	ax3.set_xlim(0.47, 0.49)
	ax3.set_ylim(34.1, 34.5)
	ax3.grid()
	ax3.set_title('Salinity box n, $S_{\mathrm{n}}$ (g kg$^{-1}$)', fontsize = 10)

	ax.set_title(r'a) AMOC strength, $\xi = 0$')

	
if p_hosing == 0.5:

	legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)
	
	ax2 	= fig.add_axes([0.175, 0.40, 0.25, 0.20])
	
	graph_1 = ax2.plot(E_A_equil[0], q_N_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
	graph_3 = ax2.plot(E_A_equil[2], q_N_equil[2], '--', color = 'gray', label='Unstable')
	graph_4 = ax2.plot(E_A_all / 10**6.0, q_N_all, '-k', linewidth = 2, label = 'Quasi-equilibrium')
	ax2.plot(E_A_all[int(num_year*365*86400/delta_t)] / 10**6.0, q_N_all[int(num_year*365*86400/delta_t)], 'ok')

	ax2.set_xlim(0.60, 0.62)
	ax2.set_ylim(7, 11)
	ax2.grid()
	ax2.set_title('AMOC strength (Sv)', fontsize = 10)
	
	ax3	= fig.add_axes([0.68, 0.40, 0.25, 0.20])
	
	graph_1 = ax3.plot(E_A_equil[0], S_n_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
	graph_3 = ax3.plot(E_A_equil[2], S_n_equil[2], '--', color = 'gray', label='Unstable')
	graph_4 = ax3.plot(E_A_all / 10**6.0, S_n_all, '-k', linewidth = 2, label = 'Quasi-equilibrium')
	ax3.plot(E_A_all[int(num_year*365*86400/delta_t)] / 10**6.0, S_n_all[int(num_year*365*86400/delta_t)], 'ok')

	ax3.set_xlim(0.60, 0.62)
	ax3.set_ylim(34.2, 34.6)
	ax3.grid()
	ax3.set_title('Salinity box n, $S_{\mathrm{n}}$ (g kg$^{-1}$)', fontsize = 10)


	ax.set_title(r'b) AMOC strength, $\xi = 0.5$')


if p_hosing == 0:
	legend_1	= ax.legend(graphs, legend_labels, loc=(0.72, 0.745), ncol=1, framealpha = 1.0, numpoints = 1)
	
	ax2 	= fig.add_axes([0.21, 0.51, 0.25, 0.20])
	
	graph_1 = ax2.plot(E_A_equil[0], q_N_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
	graph_3 = ax2.plot(E_A_equil[2], q_N_equil[2], '--', color = 'gray', label='Unstable')
	graph_4 = ax2.plot(E_A_all / 10**6.0, q_N_all, '-k', linewidth = 2, label = 'Quasi-equilibrium')
	ax2.plot(E_A_all[int(num_year*365*86400/delta_t)] / 10**6.0, q_N_all[int(num_year*365*86400/delta_t)], 'ok')

	ax2.set_xlim(0.83, 0.84)
	ax2.set_ylim(7, 11)
	ax2.grid()
	ax2.set_title('AMOC strength (Sv)', fontsize = 10)
	
	ax3	= fig.add_axes([0.21, 0.17, 0.25, 0.20])
	
	graph_1 = ax3.plot(E_A_equil[0], S_n_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
	graph_3 = ax3.plot(E_A_equil[2], S_n_equil[2], '--', color = 'gray', label='Unstable')
	graph_4 = ax3.plot(E_A_all / 10**6.0, S_n_all, '-k', linewidth = 2, label = 'Quasi-equilibrium')
	ax3.plot(E_A_all[int(num_year*365*86400/delta_t)] / 10**6.0, S_n_all[int(num_year*365*86400/delta_t)], 'ok')

	ax3.set_xlim(0.83, 0.84)
	ax3.set_ylim(34.6, 34.8)
	ax3.grid()
	ax3.set_title('Salinity box n, $S_{\mathrm{n}}$ (g kg$^{-1}$)', fontsize = 10)
	
	ax.set_title(r'c) AMOC strength, $\xi = 1.0$')


show()