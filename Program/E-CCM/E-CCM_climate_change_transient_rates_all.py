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
	return x > 0

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

E_A 			= 0.33	#Freshwater flux forcing

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
T_ts_a, T_n_a, S_t, S_n, S_ts, S_s, T_t, T_n, T_ts, T_s, T_d, D = np.loadtxt(directory+'Initialisation/'+filename_parameter)

#Assume that the tropical and southern box have the same temperatures as the ts-box and n-box, respectively
T_t_a	= np.copy(T_ts_a)
T_s_a	= np.copy(T_n_a)


delta_t			= 86400 * 5

trend_temp		= np.arange(1, 21)	#Trend per century
trend_temp[-1]		= 11.85			#Maximum forcing for AMOC recovery

time			= np.arange((1000) * 365 * 86400 / delta_t) / (365 * 86400 / delta_t)

#Save only for AMOC strength
q_N_all			= ma.masked_all((len(time), len(trend_temp)))

temp_n_anom		= ma.masked_all((len(time), len(trend_temp)))

for trend_i in range(len(trend_temp)):
	#Loop over each trend to get the trend
	temp_n_anom[:, trend_i]	= time * trend_temp[trend_i] * 0.01

temp_n_anom[temp_n_anom >= 5.0] 	= temp_n_anom[temp_n_anom >= 5.0] * 0.0 + 5.0
temp_t_anom				= np.zeros((len(time), len(trend_temp)))

#The instant forcing
temp_n_anom[:, -2] = 5.0
temp_n_anom[0, -2] = 0.0

#-----------------------------------------------------------------------------------------

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
fh = netcdf.Dataset(directory+'Steady_states/'+filename_output, 'r')

E_A_all 	= fh.variables['E_A'][:120] 	#Always start from the on state
E_A_index	= (np.abs(E_A_all - E_A)).argmin()

E_A		    = E_A_all[E_A_index] * 10**6.0	  	
S_t 		= fh.variables['S_t'][E_A_index]    	+ np.zeros(len(trend_temp))	  	
S_n 		= fh.variables['S_n'][E_A_index]    	+ np.zeros(len(trend_temp))	
S_ts		= fh.variables['S_ts'][E_A_index]   	+ np.zeros(len(trend_temp))	  	
S_s 		= fh.variables['S_s'][E_A_index]  	+ np.zeros(len(trend_temp))	  
S_d 		= fh.variables['S_d'][E_A_index] 	+ np.zeros(len(trend_temp))    	 
T_t 		= fh.variables['T_t'][E_A_index] 	+ np.zeros(len(trend_temp)) 	  
T_n 		= fh.variables['T_n'][E_A_index] 	+ np.zeros(len(trend_temp))    	  
T_ts		= fh.variables['T_ts'][E_A_index] 	+ np.zeros(len(trend_temp))     	  
T_s 		= fh.variables['T_s'][E_A_index] 	+ np.zeros(len(trend_temp))    	
T_d 		= fh.variables['T_d'][E_A_index] 	+ np.zeros(len(trend_temp))     	  	  	  
D   		= fh.variables['D'][E_A_index]  	+ np.zeros(len(trend_temp))   	  	  	  

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
	q_N_all[time_i]		= Heaviside(q_N) * q_N / 10**6.0

#-----------------------------------------------------------------------------------------
#Remove the critical trend
q_N_crit	= q_N_all[:, -1]
q_N_all		= q_N_all[:, :-1]
trend_temp	= trend_temp[:-1]


#-----------------------------------------------------------------------------------------

fig, ax = subplots()

cmap 		= mpl.cm.get_cmap('viridis', len(trend_temp))
CS	 	= ax.scatter(np.arange(len(trend_temp)), np.arange(len(trend_temp))-100, c=np.arange(len(trend_temp)), cmap=cmap)
cbar		= colorbar(CS, ticks=np.arange((len(trend_temp)-1) / len(trend_temp) / 2.0, len(trend_temp), (len(trend_temp)-1) / len(trend_temp)))
cbar.set_label('Temperature trend $T_{\mathrm{n}}^a$ ($^{\circ}$C per century)')
cbar.ax.set_yticklabels(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', 'instant', '20'])

for trend_i in range(len(trend_temp)):
	ax.plot(time, q_N_all[:, trend_i], c=cmap(trend_i))

ax.set_xlim(-20, 1000)
ax.set_ylim(-2, 18)
ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_yticks([0, 5, 10, 15])
ax.grid()

ax2 	= fig.add_axes([0.515, 0.67, 0.20, 0.15])

ax2.plot(time, q_N_crit, '-k')
ax2.set_xlim(-20, 1000)
ax2.set_ylim(-2, 18)
ax2.grid()
ax2.set_title('11.85$^{\circ}$C per century', fontsize = 10)

ax.set_title('d) AMOC strength under varying climate change ($\overline{E_A} = 0.33$ Sv)')

show()
