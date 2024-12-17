#Program plots the analytical solutions of the Stommel box model

from pylab import *
import numpy
import glob, os
import math
import matplotlib.colors as colors

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

alpha	= 2 * 10**(-4.0)
beta	= 8 * 10**(-4.0)
k	= 2 * 1000

#-----------------------------------------------------------------------------------------

delta_T	= 5

point	= k * alpha**2.0 * delta_T**2.0 / (4 * beta)

fresh_1	= np.linspace(-0.1, point, 200, endpoint = True)
fresh_2	= np.linspace(0, point, 200, endpoint = True)
fresh_4	= np.linspace(0, 0.7, 200, endpoint = True)

psi_1	= (k * alpha * delta_T / 2.0) + np.sqrt( (k * alpha * delta_T / 2.0)**2.0 - k * beta * fresh_1)
psi_2	= (k * alpha * delta_T / 2.0) - np.sqrt( (k * alpha * delta_T / 2.0)**2.0 - k * beta * fresh_2)
psi_4	= (k * alpha * delta_T / 2.0) - np.sqrt( (k * alpha * delta_T / 2.0)**2.0 + k * beta * fresh_4)

fig, ax	= subplots()

graph_1 = plot(fresh_1, psi_1, '-', color = 'darkorchid', linewidth = 1.5, label = 'AMOC on ($\psi_1$)')
graph_2 = plot(fresh_2, psi_2, ':', color = 'gray', linewidth = 1.5, label = 'Unstable ($\psi_2$)')
graph_4 = plot(fresh_4, psi_4, '-', color = 'firebrick', linewidth = 1.5, label = 'AMOC off ($\psi_4$)')
plot(fresh_1[-1], psi_1[-1], 'ok')
plot(fresh_4[0], psi_4[0], 'ok')

graphs	      	= graph_1 + graph_2 + graph_4
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_xlim(-0.05, 0.7)
ax.set_ylim(-0.5, 2.5)
ax.set_xlabel('Freshwater flux forcing, $\eta$')
ax.set_ylabel('Circulation strength, $\psi$')
ax.grid()
ax.set_title('a) Stommel model, $\Delta T^a = 5$')

#-----------------------------------------------------------------------------------------

delta_T	= 3

point	= k * alpha**2.0 * delta_T**2.0 / (4 * beta)

fresh_1	= np.linspace(-0.1, point, 200, endpoint = True)
fresh_2	= np.linspace(0, point, 200, endpoint = True)
fresh_4	= np.linspace(0, 0.7, 200, endpoint = True)

psi_1	= (k * alpha * delta_T / 2.0) + np.sqrt( (k * alpha * delta_T / 2.0)**2.0 - k * beta * fresh_1)
psi_2	= (k * alpha * delta_T / 2.0) - np.sqrt( (k * alpha * delta_T / 2.0)**2.0 - k * beta * fresh_2)
psi_4	= (k * alpha * delta_T / 2.0) - np.sqrt( (k * alpha * delta_T / 2.0)**2.0 + k * beta * fresh_4)

fig, ax	= subplots()

graph_1 = plot(fresh_1, psi_1, '-', color = 'darkorchid', linewidth = 1.5, label = 'AMOC on ($\psi_1$)')
graph_2 = plot(fresh_2, psi_2, ':', color = 'gray', linewidth = 1.5, label = 'Unstable ($\psi_2$)')
graph_4 = plot(fresh_4, psi_4, '-', color = 'firebrick', linewidth = 1.5, label = 'AMOC off ($\psi_4$)')
plot(fresh_1[-1], psi_1[-1], 'ok')
plot(fresh_4[0], psi_4[0], 'ok')

graphs	      	= graph_1 + graph_2 + graph_4
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_xlim(-0.05, 0.7)
ax.set_ylim(-0.5, 2.5)
ax.set_xlabel('Freshwater flux forcing, $\eta$')
ax.set_ylabel('Circulation strength, $\psi$')
ax.grid()
ax.set_title('b) Stommel model, $\Delta T^a = 3$')

show()




