#The E-CCM (temperature only)

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


def ECCM(time_all, E_A_all, S_t, S_n, S_ts, S_s, S_d, T_t, T_n, T_ts, T_s, T_d, D):

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

	#Horizontal areas for the exchange with the surface
	A_t	= np.copy(A)
	A_ts	= L_y * L_xA
	A_s	= L_y * L_xS
	A_n	= np.copy(A_ts)

	#Empty arrays for each time step
	S_t_all		= np.ma.masked_all(len(time_all))
	S_n_all		= np.ma.masked_all(len(time_all))
	S_ts_all	= np.ma.masked_all(len(time_all))	#Salinity
	S_s_all		= np.ma.masked_all(len(time_all))
	S_d_all		= np.ma.masked_all(len(time_all))

	T_t_all		= np.ma.masked_all(len(time_all))
	T_n_all		= np.ma.masked_all(len(time_all))
	T_ts_all	= np.ma.masked_all(len(time_all))	#Temperature
	T_s_all		= np.ma.masked_all(len(time_all))
	T_d_all		= np.ma.masked_all(len(time_all))

	D_all		= np.ma.masked_all(len(time_all))
	q_N_all		= np.ma.masked_all(len(time_all))	#Remaining quantities
	q_S_all		= np.ma.masked_all(len(time_all))
	M_ov_all	= np.ma.masked_all(len(time_all))

	for time_i in range(len(time_all)):
		#Add stochastic noise
		E_A	= E_A_all[time_i] #+ np.random.normal(0, 1.5) * 10**6.0

		#T-box
		VS_t_delta 	= q_S * (Heaviside(q_S) * S_ts + Heaviside(-q_S) * S_t) + q_U * S_d - Heaviside(q_N) * q_N * S_t + r_S * (S_ts - S_t) + r_N * (S_n - S_t) + 2.0 * E_S * S_0
		OHC_t_delta 	= q_S * (Heaviside(q_S) * (T_ts+273.15) + Heaviside(-q_S) * (T_t+273.15)) + q_U * (T_d+273.15) - Heaviside(q_N) * q_N * (T_t+273.15) + r_S * (T_ts - T_t) + r_N * (T_n - T_t) - A_t * lambda_a * (T_t - T_t_a)
		VS_t		= VS_t + VS_t_delta * delta_t
		OHC_t		= OHC_t + OHC_t_delta * delta_t

		#TS-box
		VS_ts_delta	= q_Ek * S_s - q_e * S_ts - q_S * (Heaviside(q_S) * S_ts + Heaviside(-q_S) * S_t) + r_S * (S_t - S_ts)
		OHC_ts_delta	= q_Ek * (T_s+273.15) - q_e * (T_ts+273.15) - q_S * (Heaviside(q_S) * (T_ts+273.15) + Heaviside(-q_S) * (T_t+273.15)) + r_S * (T_t - T_ts) - A_ts * lambda_a * (T_ts - T_ts_a)
		VS_ts		= VS_ts + VS_ts_delta * delta_t
		OHC_ts		= OHC_ts + OHC_ts_delta * delta_t

		#N-box
		VS_n_delta	= Heaviside(q_N) * q_N * (S_t - S_n) + r_N * (S_t - S_n) - (E_S + E_A) * S_0
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
	
	return S_t_all, S_n_all, S_ts_all, S_s_all, S_d_all, T_t_all, T_n_all, T_ts_all, T_s_all, T_d_all, D_all, q_N_all, M_ov_all
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
lambda_a= 0.000004		#Atmospheric exchange coefficient (per time)
eta	    = 3.0 * 10**4.0		#Hydraulic constant
alpha	= 2.0 * 10**(-4.0)	#Thermal expansion coefficient
beta	= 8.0 * 10**(-4.0)	#Haline contraction coefficient
r_S	    = 1.0 * 10**7.0		#Transport by the southern subtropical gyre
r_N	    = 5.0 * 10**6.0		#Transport by the northern subtropical gyre
E_S	    = 0.17 * 10**6.0	#Symmetric freshwater flux
E_A	    = 0.0644 * 10**6.0	#Asymmetric freshwater flux
#-----------------------------------------------------------------------------------------

lambda_a		    = 0.0000035		#Atmospheric exchange coefficient (per time)
filename_parameter	= 'MOC_box_model_lambda_a_0_0000035.txt'
filename_output		= 'MOC_box_model_temp_equilibrium_states_lambda_a_0_0000035.nc'

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
T_ts_a, T_n_a, S_t, S_n, S_ts, S_s, T_t, T_n, T_ts, T_s, T_d, D = np.loadtxt(directory+'Initialisation/'+filename_parameter)

#Assume that the tropical and southern box have the same temperatures as the ts-box and n-box, respectively
T_t_a	= np.copy(T_ts_a)
T_s_a	= np.copy(T_n_a)

#Get the steady state and use these as input parameters, but these are freely allowed to evolve
fh = netcdf.Dataset(directory+'Steady_states/'+filename_output, 'r')

S_t = fh.variables['S_t'][0]     	  	
S_n = fh.variables['S_n'][0]     	  	
S_ts= fh.variables['S_ts'][0]     	  	
S_s = fh.variables['S_s'][0]     	  
S_d = fh.variables['S_d'][0]     	 
T_t = fh.variables['T_t'][0]     	  
T_n = fh.variables['T_n'][0]     	  
T_ts= fh.variables['T_ts'][0]     	  
T_s = fh.variables['T_s'][0]     	
T_d = fh.variables['T_d'][0]     	  	  	  
D   = fh.variables['D'][0]     	  	  	  

fh.close()


#-----------------------------------------------------------------------------------------

delta_t		= 86400 * 5
time_hys	= np.arange(4400 * 365 * 86400 / delta_t) / (365 * 86400 / delta_t)

E_A_hys		= 0.66-np.abs((time_hys - 2200)) * 0.0003
E_A_hys		= E_A_hys * 10**6.0

S_t_hys, S_n_hys, S_ts_hys, S_s_hys, S_d_hys, T_t_hys, T_n_hys, T_ts_hys, T_s_hys, T_d_hys, D_hys, q_N_hys, F_ov_hys = ECCM(time_hys, E_A_hys, S_t, S_n, S_ts, S_s, S_d, T_t, T_n, T_ts, T_s, T_d, D)

for time_branch_i in [600, 1500, 1530, 1560, 1590, 1620]:

	time_index	= (np.abs(time_hys-time_branch_i)).argmin()

	time_branch	= np.arange(500 * 365 * 86400 / delta_t) / (365 * 86400 / delta_t)
	E_A_branch 	= np.zeros(len(time_branch))+E_A_hys[time_index]

	if time_branch_i == 1530:

		S_t_branch_2, S_n_branch_2, S_ts_branch_2, S_s_branch_2, S_d_branch_2, T_t_branch_2, T_n_branch_2, T_ts_branch_2, T_s_branch_2, T_d_branch_2, D_branch_2, q_N_branch_2, F_ov_branch_2 = ECCM(time_branch, E_A_branch, S_t_branch_1[-1], S_n_branch_1[-1], S_ts_branch_1[-1], S_s_branch_1[-1], S_d_branch_1[-1], T_t_branch_1[-1], T_n_branch_1[-1], T_ts_branch_1[-1], T_s_branch_1[-1], T_d_branch_1[-1], D_branch_1[-1])

	if time_branch_i > 1530:

		S_t_branch_2, S_n_branch_2, S_ts_branch_2, S_s_branch_2, S_d_branch_2, T_t_branch_2, T_n_branch_2, T_ts_branch_2, T_s_branch_2, T_d_branch_2, D_branch_2, q_N_branch_2, F_ov_branch_2 = ECCM(time_branch, E_A_branch, S_t_branch_2[-1], S_n_branch_2[-1], S_ts_branch_2[-1], S_s_branch_2[-1], S_d_branch_2[-1], T_t_branch_2[-1], T_n_branch_2[-1], T_ts_branch_2[-1], T_s_branch_2[-1], T_d_branch_2[-1], D_branch_2[-1])

	S_t_branch_1, S_n_branch_1, S_ts_branch_1, S_s_branch_1, S_d_branch_1, T_t_branch_1, T_n_branch_1, T_ts_branch_1, T_s_branch_1, T_d_branch_1, D_branch_1, q_N_branch_1, F_ov_branch_1 = ECCM(time_branch, E_A_branch, S_t_hys[time_index], S_n_hys[time_index], S_ts_hys[time_index], S_s_hys[time_index], S_d_hys[time_index], T_t_hys[time_index], T_n_hys[time_index], T_ts_hys[time_index], T_s_hys[time_index], T_d_hys[time_index], D_hys[time_index])

	#-----------------------------------------------------------------------------------------

	fig, ax	= subplots()

	ax.plot(time_hys, q_N_hys, '-k', linewidth = 2, zorder = 6)

	if time_branch_i >= 1560:

		ax.plot(time_branch+time_branch_i, q_N_branch_1, '-', linewidth = 2, color = 'royalblue', zorder = 7)
		ax.plot(time_branch+time_branch_i, q_N_branch_2, '-', linewidth = 2, color = 'firebrick', zorder = 7)

	else:
		ax.plot(time_branch+time_branch_i, q_N_branch_1, '-', linewidth = 2, color = 'royalblue', zorder = 7)
		
	ax.set_xlim(1480, 2220)
	ax.set_ylim(-2, 18)
	ax.set_yticks([0, 5, 10, 15])
	ax.set_xlabel('Model year')
	ax.set_ylabel('Volume transport (Sv)')
	ax.grid()
	
	if time_branch_i <= 1530:
		graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 1.5, label = 'Quasi-equilibrium (QE)')
		graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 1.5, label = 'Branched from QE')

	else:
		graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 1.5, label = 'Quasi-equilibrium (QE)')
		graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 1.5, label = 'Branched from QE')
		graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'Branched from $\overline{E_A} = $'+str(np.round((time_branch_i-10)*0.0003, 3))+' Sv')

	if time_branch_i == 600:
		ax.set_xlim(580, 1320)	
		legend_1	= ax.legend(loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(5)
	
	else:
		legend_1	= ax.legend(loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(5)

	if time_branch_i == 600: subplot_letter='c)'
	if time_branch_i == 1500: subplot_letter='d)'
	if time_branch_i == 1530: subplot_letter='e)'
	if time_branch_i == 1560: subplot_letter='f)'
	if time_branch_i == 1590: subplot_letter='g)'
	if time_branch_i == 1620: subplot_letter='h)'
	
	ax.set_title(subplot_letter+' AMOC strength, $\overline{E_A} = $'+str(np.round(time_branch_i*0.0003, 3))+' Sv (model year '+str(time_branch_i)+')')
	
	#-----------------------------------------------------------------------------------------

	fig, ax	= subplots()

	ax.plot(time_hys, F_ov_hys, '-k', linewidth = 2, zorder = 6)

	if time_branch_i >= 1560:

		ax.plot(time_branch+time_branch_i, F_ov_branch_1, '-', linewidth = 2, color = 'royalblue', zorder = 7)
		ax.plot(time_branch+time_branch_i, F_ov_branch_2, '-', linewidth = 2, color = 'firebrick', zorder = 7)

	else:
		ax.plot(time_branch+time_branch_i, F_ov_branch_1, '-', linewidth = 2, color = 'royalblue', zorder = 7)
		
	ax.set_xlim(1480, 2220)
	ax.set_ylim(-0.4, 0.4)
	#ax.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
	ax.set_xlabel('Model year')
	ax.set_ylabel('Freshwater transport (Sv)')
	ax.grid()
	
	if time_branch_i <= 1530:
		graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 1.5, label = 'Quasi-equilibrium (QE)')
		graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 1.5, label = 'Branched from QE')

	else:
		graph_1		= ax.plot([-100, -100], [-100, -100], '-', color = 'k', linewidth = 1.5, label = 'Quasi-equilibrium (QE)')
		graph_2		= ax.plot([-100, -100], [-100, -100], '-', color = 'royalblue', linewidth = 1.5, label = 'Branched from QE')
		graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'Branched from $\overline{E_A} = $'+str(np.round((time_branch_i-10)*0.0003, 3))+' Sv')

	if time_branch_i == 600:
		ax.set_xlim(580, 1320)	
		
	legend_1	= ax.legend(loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1).set_zorder(5)

	if time_branch_i == 600: subplot_letter='i)'
	if time_branch_i == 1500: subplot_letter='j)'
	if time_branch_i == 1530: subplot_letter='k)'
	if time_branch_i == 1560: subplot_letter='l)'
	if time_branch_i == 1590: subplot_letter='m)'
	if time_branch_i == 1620: subplot_letter='n)'
	
	ax.set_title(subplot_letter+' $F_{\mathrm{ovS}}$, $\overline{E_A} = $'+str(np.round(time_branch_i*0.0003, 3))+' Sv (model year '+str(time_branch_i)+')')
	
#-----------------------------------------------------------------------------------------

filename_con	= directory+'Continuation/Hosing_location/Continuation_MOC_box_model_temp_lambda_a_0_000035_p_hosing_100.nc'

fh = netcdf.Dataset(filename_con, 'r')

#Writing data to correct variable	
E_A_equil	= fh.variables['E_A'][:]     	  	
q_N_equil	= fh.variables['q_n'][:]          	  	
F_ov_equil	= fh.variables['F_ovS'][:]     

fh.close()

print('E_A max:', np.max(E_A_equil[0]), 'Sv')
print('FovS min:', E_A_equil[0, np.argmin(F_ov_equil[0])])
print('diff E_A:', np.max(E_A_equil[0]) - E_A_equil[0, np.argmin(F_ov_equil[0])])

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1 = ax.plot(E_A_equil[0], q_N_equil[0], '-', linewidth = 2.0, color = 'darkorchid')
graph_2 = ax.plot(E_A_equil[2], q_N_equil[2], ':', linewidth = 2.0, color = 'gray')
graph_3 = ax.plot(E_A_equil[1], q_N_equil[1], '-', linewidth = 2.0,  color = 'firebrick')

ax.plot((4400 - time_hys)*0.0003, q_N_hys, '-r', linewidth = 1.5)
ax.plot(time_hys*0.0003, q_N_hys, '-k', linewidth = 1.5)

ax.set_xlabel('Freshwater flux forcing, $E_A$ (Sv)')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(0, 0.66)
ax.set_ylim(-2, 35)
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Forward quasi-equilibrium')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = 'Backward quasi-equilibrium')
graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'darkorchid', linewidth = 1.5, label = 'AMOC on (steady state)')
graph_4		= ax.plot([-100, -100], [-100, -100], ':', color = 'gray', linewidth = 1.5, label = 'Unstable (steady state)')
graph_5		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'AMOC off (steady state)')

graphs	      	= graph_1 + graph_2 + graph_3 + graph_4 + graph_5
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.plot([0.205, 0.295], [22.5, 22.5], '-k')
ax.plot([0.205, 0.205], [22, 23], '-k')
ax.plot([0.295, 0.295], [22, 23], '-k')
ax.text(0.250, 24, '300 years', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize = 10)

ax.set_title('a) AMOC strength')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_1 = ax.plot(E_A_equil[0], F_ov_equil[0], '-', linewidth = 2.0, color = 'darkorchid')
graph_2 = ax.plot(E_A_equil[2], F_ov_equil[2], ':', linewidth = 2.0, color = 'gray')
graph_3 = ax.plot(E_A_equil[1], F_ov_equil[1], '-', linewidth = 2.0,  color = 'firebrick')

ax.plot((4400 - time_hys)*0.0003, F_ov_hys, '-r', linewidth = 1.5)
ax.plot(time_hys*0.0003, F_ov_hys, '-k', linewidth = 1.5)

ax.set_xlabel('Freshwater flux forcing, $E_A$ (Sv)')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(0, 0.66)
ax.set_ylim(-0.4, 0.4)
ax.grid()

graph_1		= ax.plot([-100, -100], [-100, -100], '-k', linewidth = 1.5, label = 'Forward quasi-equilibrium')
graph_2		= ax.plot([-100, -100], [-100, -100], '-r', linewidth = 1.5, label = 'Backward quasi-equilibrium')
graph_3		= ax.plot([-100, -100], [-100, -100], '-', color = 'darkorchid', linewidth = 1.5, label = 'AMOC on (steady state)')
graph_4		= ax.plot([-100, -100], [-100, -100], ':', color = 'gray', linewidth = 1.5, label = 'Unstable (steady state)')
graph_5		= ax.plot([-100, -100], [-100, -100], '-', color = 'firebrick', linewidth = 1.5, label = 'AMOC off (steady state)')

graphs	      	= graph_1 + graph_2 + graph_3 + graph_4 + graph_5
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) Freshwater transport at 30$^{\circ}$S, $F_{\mathrm{ovS}}$')

show()


