#The E-CCM (temperature only) and under climate change

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

def SteadyStates(filename_continuation, E_A_value):
	"""Get the AMOC on, AMOC off and the unstable branch
	for the given value of E_A"""

	variable_names	= ['S_t', 'S_n', 'S_ts', 'S_s', 'S_d', 'T_t', 'T_n', 'T_ts', 'T_s', 'T_d', 'q_n', 'F_ovS', 'D']

	variable_branches	= ma.masked_all((len(variable_names), 3))

	for variable_i in range(len(variable_names)):
		#Loop over each variable and determine the equilibrium value for the given E_A
		fh = netcdf.Dataset(filename_continuation, 'r')

		E_A		= fh.variables['E_A'][:]     	  	
		variable	= fh.variables[variable_names[variable_i]][:]     	  	  	

		fh.close()

		for branch_i in range(3):
			#Loop over the on branch, off branch, and unstable branch
			E_A_branch, variable_branch	= E_A[branch_i], variable[branch_i]
			E_A_branch, variable_branch	= E_A_branch[E_A_branch.mask == False], variable_branch[E_A_branch.mask == False]
			sort_index			= np.argsort(E_A_branch)
			E_A_branch, variable_branch	= E_A_branch[sort_index], variable_branch[sort_index]

			#Interpolate
			variable_int			= np.interp(E_A_value, E_A_branch, variable_branch, left = -999, right = -999)

			if variable_int == -999:
				#Variable is outside the relevant range
				continue

			#Save the variable
			variable_branches[variable_i, branch_i]	= variable_int
  	  
	return variable_branches

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
E_A	    = 0.0644 * 10**6.0	#Asymmetric freshwater flux
#-----------------------------------------------------------------------------------------

p_hosing		    = 1.0             #xi = 1 - p
lambda_a		    = 0.0000035		#Atmospheric exchange coefficient (per time)
filename_parameter	= 'MOC_box_model_lambda_a_0_0000035.txt'
filename_output		= 'MOC_box_model_temp_equilibrium_states_lambda_a_0_0000035.nc'

E_A 			    = 0.335	#Freshwater flux forcing

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
T_ts_a, T_n_a, S_t, S_n, S_ts, S_s, T_t, T_n, T_ts, T_s, T_d, D = np.loadtxt(directory+'Initialisation/'+filename_parameter)

#Assume that the tropical and southern box have the same temperatures as the ts-box and n-box, respectively
T_t_a	= np.copy(T_ts_a)
T_s_a	= np.copy(T_n_a)

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
fh = netcdf.Dataset(directory+'Steady_states/'+filename_output, 'r')

E_A_all 	= fh.variables['E_A'][:120] 	#Always start from the on state
E_A_index	= (np.abs(E_A_all - E_A)).argmin()

E_A		    = E_A_all[E_A_index] * 10**6.0	  	
S_t 		= fh.variables['S_t'][E_A_index]     	  	
S_n 		= fh.variables['S_n'][E_A_index]     	  	
S_ts		= fh.variables['S_ts'][E_A_index]     	  	
S_s 		= fh.variables['S_s'][E_A_index]     	  
S_d 		= fh.variables['S_d'][E_A_index]     	 
T_t 		= fh.variables['T_t'][E_A_index]     	  
T_n 		= fh.variables['T_n'][E_A_index]     	  
T_ts		= fh.variables['T_ts'][E_A_index]     	  
T_s 		= fh.variables['T_s'][E_A_index]     	
T_d 		= fh.variables['T_d'][E_A_index]     	  	  	  
D   		= fh.variables['D'][E_A_index]     	  	  	  

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

num_year		= 500
time			= np.arange((1000) * 365 * 86400 / delta_t) / (365 * 86400 / delta_t)

#Empty arrays for each time step
S_t_all		= np.ma.masked_all((len(time), 4))
S_n_all		= np.ma.masked_all((len(time), 4))
S_ts_all	= np.ma.masked_all((len(time), 4))	#Salinity
S_s_all		= np.ma.masked_all((len(time), 4))
S_d_all		= np.ma.masked_all((len(time), 4))

T_t_all		= np.ma.masked_all((len(time), 4))
T_n_all		= np.ma.masked_all((len(time), 4))
T_ts_all	= np.ma.masked_all((len(time), 4))	#Temperature
T_s_all		= np.ma.masked_all((len(time), 4))
T_d_all		= np.ma.masked_all((len(time), 4))

D_all		= np.ma.masked_all((len(time), 4))
q_N_all		= np.ma.masked_all((len(time), 4))	#Remaining quantities
q_S_all		= np.ma.masked_all((len(time), 4))
M_ov_all	= np.ma.masked_all((len(time), 4))


temp_n_anom					= np.zeros(len(time))
temp_n_anom[:int(num_year*365*86400/delta_t)]	= np.linspace(0, 5.0, int(num_year*365*86400/delta_t), endpoint = True)
temp_n_anom[int(num_year*365*86400/delta_t):]	= 5.0
temp_t_anom					= np.zeros(len(time))

for time_i in range(len(time)):
	#Add stochastic noise

	#T-box
	VS_t_delta 	= q_S * (Heaviside(q_S) * S_ts + Heaviside(-q_S) * S_t) + q_U * S_d - Heaviside(q_N) * q_N * S_t + r_S * (S_ts - S_t) + r_N * (S_n - S_t) + (2.0 * E_S - (1.0 - p_hosing) * E_A) * S_0
	OHC_t_delta 	= q_S * (Heaviside(q_S) * (T_ts+273.15) + Heaviside(-q_S) * (T_t+273.15)) + q_U * (T_d+273.15) - Heaviside(q_N) * q_N * (T_t+273.15) + r_S * (T_ts - T_t) + r_N * (T_n - T_t) - A_t * lambda_a * (T_t - (T_t_a+temp_t_anom[time_i]))
	VS_t		= VS_t + VS_t_delta * delta_t
	OHC_t		= OHC_t + OHC_t_delta * delta_t

	#TS-box
	VS_ts_delta	= q_Ek * S_s - q_e * S_ts - q_S * (Heaviside(q_S) * S_ts + Heaviside(-q_S) * S_t) + r_S * (S_t - S_ts)
	OHC_ts_delta	= q_Ek * (T_s+273.15) - q_e * (T_ts+273.15) - q_S * (Heaviside(q_S) * (T_ts+273.15) + Heaviside(-q_S) * (T_t+273.15)) + r_S * (T_t - T_ts) - A_ts * lambda_a * (T_ts - (T_ts_a+temp_t_anom[time_i]))
	VS_ts		= VS_ts + VS_ts_delta * delta_t
	OHC_ts		= OHC_ts + OHC_ts_delta * delta_t

	#N-box
	VS_n_delta	= Heaviside(q_N) * q_N * (S_t - S_n) + r_N * (S_t - S_n) - (E_S + p_hosing * E_A) * S_0
	OHC_n_delta	= Heaviside(q_N) * q_N * (T_t - T_n) + r_N * (T_t - T_n) - A_n * lambda_a * (T_n - (T_n_a+temp_n_anom[time_i]))
	VS_n		= VS_n + VS_n_delta * delta_t
	OHC_n		= OHC_n + OHC_n_delta * delta_t

	#S-box
	VS_s_delta	= q_S * (Heaviside(q_S) * S_d + Heaviside(-q_S) * S_s) + q_e * S_ts - q_Ek * S_s - (E_S - E_A) * S_0
	OHC_s_delta	= q_S * (Heaviside(q_S) * (T_d+273.15) + Heaviside(-q_S) * (T_s+273.15)) + q_e * (T_ts+273.15) - q_Ek * (T_s+273.15) - A_s * lambda_a * (T_s - (T_s_a+temp_n_anom[time_i]))
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
	S_t_all[time_i, 0]		= S_t
	S_n_all[time_i, 0]		= S_n
	S_ts_all[time_i, 0]		= S_ts	#Salinity
	S_s_all[time_i, 0]		= S_s
	S_d_all[time_i, 0]		= S_d
	T_t_all[time_i, 0]		= T_t
	T_n_all[time_i, 0]		= T_n
	T_ts_all[time_i, 0]		= T_ts	#Temperature
	T_s_all[time_i, 0]		= T_s
	T_d_all[time_i, 0]		= T_d
	D_all[time_i, 0]		= D
	q_N_all[time_i, 0]		= Heaviside(q_N) * q_N / 10**6.0
	q_S_all[time_i, 0]		= q_S / 10**6.0
	M_ov_all[time_i, 0]		= M_ov / 10**6.0

	#Get the relevant steady states for each parameter
	T_n_str		= str(round(temp_n_anom[time_i], 1))

	if time_i == 0:
		filename_con	= directory+'Continuation/Climate_change/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'_T_t_0_0_T_n_'+str(T_n_str[0])+'_'+str(T_n_str[2])+'.nc'

		steady_states	= SteadyStates(filename_con, E_A / 10**6.0)

	if filename_con	!= directory+'Continuation/Climate_change/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'_T_t_0_0_T_n_'+str(T_n_str[0])+'_'+str(T_n_str[2])+'.nc':
		#Different steady states
		steady_states	= SteadyStates(filename_con, E_A / 10**6.0)
	
		#Update continuation file
		filename_con	= directory+'Continuation/Climate_change/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'_T_t_0_0_T_n_'+str(T_n_str[0])+'_'+str(T_n_str[2])+'.nc'


	S_t_all[time_i, 1:]	= steady_states[0]
	S_n_all[time_i, 1:]	= steady_states[1]
	S_ts_all[time_i, 1:]	= steady_states[2]
	S_s_all[time_i, 1:]	= steady_states[3]
	S_d_all[time_i, 1:]	= steady_states[4]
	T_t_all[time_i, 1:]	= steady_states[5]
	T_n_all[time_i, 1:]	= steady_states[6]
	T_ts_all[time_i, 1:]	= steady_states[7]
	T_s_all[time_i, 1:]	= steady_states[8]
	T_d_all[time_i, 1:]	= steady_states[9]
	q_N_all[time_i, 1:]	= steady_states[10]
	M_ov_all[time_i, 1:]	= steady_states[11]
	D_all[time_i, 1:]	= steady_states[12]


#-----------------------------------------------------------------------------------------


fig, ax = subplots()

ax.plot(time, q_N_all[:, 0], '-k')
ax.plot(time, q_N_all[:, 1], '-', color = 'darkviolet')
ax.plot(time, q_N_all[:, 2], '-', color = 'firebrick')
ax.plot(time, q_N_all[:, 3], ':', color = 'gray')

ax.set_xlim(-20, 1000)
ax.set_ylim(-2, 18)
ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_yticks([0, 5, 10, 15])
ax.grid()

graph_1  = ax.plot(0, -100, '-', color = 'k', linewidth = 2.0, label = 'Transient')
graph_2  = ax.plot(0, -100, '-', color = 'darkviolet', linewidth = 2.0, label = 'AMOC on')
graph_3  = ax.plot(0, -100, ':', color = 'gray', linewidth = 2.0, label = 'Unstable')
graph_4  = ax.plot(0, -100, '-', color = 'firebrick', linewidth = 2.0, label = 'AMOC off')

graphs_1	= graph_1 + graph_2 + graph_3 + graph_4

legend_labels_2 = [l.get_label() for l in graphs_1]
legend_1	= ax.legend(graphs_1, legend_labels_2, loc = 'upper right', ncol=1, framealpha = 1.0)

ax.set_title('b) AMOC strength under climate change ($\overline{E_A} = 0.335$ Sv)')

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Continuation/Climate_change/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'_T_t_0_0_T_n_2_0.nc', 'r')

#Writing data to correct variable	
E_A_equil		= fh.variables['E_A'][:]     	  	
q_N_equil		= fh.variables['q_n'][:]  

fh.close()

ax2 	= fig.add_axes([0.175, 0.45, 0.20, 0.15])

ax2.plot(E_A_equil[2], q_N_equil[2], ':', color = 'gray', label='Unstable')
ax2.plot(E_A_equil[0], q_N_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
ax2.plot(E_A_equil[1], q_N_equil[1], '-', color = 'firebrick', label = 'AMOC off')
ax2.plot(E_A / 10**6.0, q_N_all[int(200*365*86400/delta_t), 0], 'ok')

ax2.set_xlim(0, 0.60)
ax2.set_ylim(-2, 18)
ax2.grid()
ax2.set_title('$\Delta T_{\mathrm{n}}^a = 2^{\circ}$C', fontsize = 10)

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Continuation/Climate_change/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_'+str(int(p_hosing*100)).zfill(3)+'_T_t_0_0_T_n_4_0.nc', 'r')

#Writing data to correct variable	
E_A_equil		= fh.variables['E_A'][:]     	  	
q_N_equil		= fh.variables['q_n'][:]  

fh.close()

ax3 	= fig.add_axes([0.50, 0.25, 0.20, 0.15])

ax3.plot(E_A_equil[2], q_N_equil[2], ':', color = 'gray', label='Unstable')
ax3.plot(E_A_equil[0], q_N_equil[0], '-', color = 'darkorchid', label = 'AMOC on')
ax3.plot(E_A_equil[1], q_N_equil[1], '-', color = 'firebrick', label = 'AMOC off')
ax3.plot(E_A / 10**6.0, q_N_all[int(400*365*86400/delta_t), 0], 'ok')

ax3.set_xlim(0, 0.60)
ax3.set_ylim(-2, 18)
ax3.grid()
ax3.set_title('$\Delta T_{\mathrm{n}}^a = 4^{\circ}$C', fontsize = 10)
	

show()