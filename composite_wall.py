"""
    Composite wall.
    
    Simulation made as part of a project delivery in 
    Transport Processes (PEF1110) at the University of South-Eastern Norway
"""

import numpy as np
import math as mt
import matplotlib.pyplot as plt

# Parameters
Ti = 20                                     # Initial time (x = 0 and t < 0) [°C]
Ts = 30                                     # Change of temperature at startup (x = 0 and t >= 0) [°C]
Tend = 10                                   # Temperature of end of the wall (x = L) [°C]

L = 0.4                                     # Length of wall [m]
L_cs = 0.1                                  # Length of carbon steel (first and last section of wall) [m]
L_al = 0.2                                  # Length of aluminium (middle section) [m]

interface_nodes = [L_cs, L_cs + L_al]       # Interface nodes

# Carbon steel
k_cs = 45                                   # Thermal conductivity of carbon steel [W/(mK)]
rho_cs = 7850                               # Density of carbon steel [kg/m³]
cp_cs = 500                                 # Specific heat caoacity of carcon steel [J/(kgK)]
alpha_cs = k_cs/(rho_cs*cp_cs)              # Thermal diffusivity of carbon steel [W/(mK)]

# Aluminium
k_al = 205                                  # Thermal conductivity of aluminium steel [W/(mK)]
rho_al = 2700                               # Density of aluminium [kg/m³]
cp_al = 897                                 # Specific heat caoacity of aluminium [J/(kgK)]
alpha_al = k_al/(rho_al*cp_al)              # Thermal diffusivity of aluminium [W/(mK)]

avg_rho_cp = (rho_cs*cp_cs + rho_al*cp_al)/2

# Simulation parameters
dt = 0.05                                   # Timestep [s]
end_time = 7200                             # End time of simulation [s]
N_t = mt.floor(end_time/dt) + 1             # Number of time instances

N_x = 101                                   # Number of equidistant points in x-direction
dx = L/(N_x - 1)                            # Distance between adjacent spatial points being evaluated

# Storage arrays
temp_init = np.zeros(N_x)
position = np.linspace(0, L, N_x)

# Calculate initial temperature profile: steady-state
temp_array = np.zeros((N_t, N_x))
heat_flux = -(Tend - Ti)/(L_cs/k_cs + L_al/k_al + L_cs/k_cs)
T1 = Ti - heat_flux*L_cs/k_cs
T2 = T1 - heat_flux*L_al/k_al

first_index = round(L_cs/dx)                # Index position of first wall interface
second_index = round((L_cs + L_al)/dx)      # Index position of second wall interface

for i in range(N_x):
    if (i == 0):
        temp_init[i] = Ti
    elif (i <= first_index):
        temp_init[i] = temp_init[i-1] + (T1 - Ti)/L_cs*dx
    elif (i <= second_index):
        temp_init[i] = temp_init[i-1] + (T2 - T1)/L_al*dx
    else:
        temp_init[i] = temp_init[i-1] + (T1 - Ti)/L_cs*dx

temp_array[0,:] = temp_init

# Solving equations
for i in range(1, N_t):
    for j in range (N_x):
        
        # Which heat diffusivity to use
        pos = j*dx
        if (pos >= L_cs and pos <= (L_al + L_cs)):
            alpha = alpha_al
        else:
            alpha = alpha_cs
        
        if (j == 0):
            temp_array[i, j] = Ts
        elif (j == N_x-1):
            temp_array[i, j] = Tend
        else:
            if (mt.isclose(pos, interface_nodes[0])):
                # Must handle carbon steel - aluminium interface node
                temp_array[i, j] = temp_array[i-1, j] + dt/(avg_rho_cp*dx**2)*(k_al*(temp_array[i-1, j+1] - temp_array[i-1, j]) - k_cs*(temp_array[i-1, j]-temp_array[i-1, j-1]))
            elif (mt.isclose(pos, interface_nodes[1])):
                # Must handle aluminium - carbon steel interface node
                temp_array[i, j] = temp_array[i-1, j] + dt/(avg_rho_cp*dx**2)*(k_cs*(temp_array[i-1, j+1] - temp_array[i-1, j]) - k_al*(temp_array[i-1, j]-temp_array[i-1, j-1]))
            else:
                # Solving non-interface nodes
                temp_array[i, j] = temp_array[i-1, j] + alpha*dt/(dx**2)*(temp_array[i-1, j+1] - 2*temp_array[i-1, j] + temp_array[i-1, j-1])

 # Add colormap
num_plots = 7
plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.PuRd(np.linspace(0.25, 1, num_plots))))
            
# Plot temperature distribution for different time instances
n = mt.floor(0/dt)
plt.plot(position, temp_array[n,:], label="t < 0 min")
n = mt.floor(1*60/dt)
plt.plot(position, temp_array[n,:], label="t = 1 min")
n = mt.floor(5*60/dt)
plt.plot(position, temp_array[n,:], label="t = 5 min")
n = mt.floor(10*60/dt)
plt.plot(position, temp_array[n,:], label="t = 10 min")
n = mt.floor(30*60/dt)
plt.plot(position, temp_array[n,:], label="t = 30 min")
n = mt.floor(60*60/dt)
plt.plot(position, temp_array[n,:], label="t = 60 min")
plt.plot(position, temp_array[-1,:], label="t = 120 min")

# Add vertical axes at the material interfaces
plt.axvline(x = 0, color='gray', linestyle='--', alpha = 0.7)
plt.axvline(x = L_cs, color='gray', linestyle='--', alpha = 0.7)
plt.axvline(x = L_cs + L_al, color = 'gray', linestyle = '--', alpha = 0.7)
plt.axvline(x = L_cs + L_al + L_cs, color = 'gray', linestyle = '--', alpha = 0.7)

plt.text(L_cs/2, 28, "Cs", 
         horizontalalignment='center', color='dimgray', fontsize=10, weight='bold')

plt.text(L_cs + L_al/2, 28, "Al", 
         horizontalalignment='center', color='dimgray', fontsize=10, weight='bold')

plt.text(L_cs + L_al + L_cs/2, 18, "Cs", 
         horizontalalignment='center', color='dimgray', fontsize=10, weight='bold')

plt.grid(True, linestyle=":", alpha=0.6)

plt.title("Composite wall")
plt.xlabel("Position [m]")
plt.ylabel("Temperature [°C]")
plt.legend()
plt.show()
