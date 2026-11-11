"""
    Single carbon steel wall
    
    Simulation made as part of a project delivery in 
    Transport Processes (PEF1110) at the University of South-Eastern Norway
"""

import numpy as np
import math as mt
import matplotlib.pyplot as plt

# Parameters
Ti = 20                             # Initial time (x = 0 and t < 0) [°C]
Ts = 30                             # Change of temperature at startup (x = 0 and t >= 0) [°C]
Tend = 10                           # Temperature of end of the wall (x = L) [°C]

# Carbon steel
k = 45                              # Thermal conductivity of carbon steel [W/(mK)]
rho = 7850                          # Density of carbon steel [kg/m³]
cp = 500                            # Specific heat caoacity of carcon steel [J/(kgK)]
L = 0.4                             # Length of wall [m]

alpha = k/(rho*cp)                  # Thermal diffusivity [m^2/s]

# Simulation parameters
dt = 0.05                           # Timestep [s]
end_time = 7200                     # End time of simulation [s]
N_t = mt.floor(end_time/dt) + 1     # Number of time instances

N_x = 101                           # Number of equidistant points in x-direction
dx = L/(N_x - 1)                    # Distance between adjacent spatial points being evaluated

# Storage arrays
temp_init = np.linspace(Ti, Tend, N_x)
position = np.linspace(0, L, N_x)

temp_array = np.zeros((N_t, N_x))
temp_array[0,:] = temp_init

# Solving equations
for i in range(1, N_t):
    for j in range (N_x):
        if (j == 0):
            temp_array[i, j] = Ts
        elif (j == N_x-1):
            temp_array[i, j] = Tend
        else:
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

plt.axvline(x = 0, color='gray', linestyle='--', alpha = 0.7)
plt.axvline(x = L, color='gray', linestyle='--', alpha = 0.7)

plt.grid(True, linestyle=":", alpha=0.6)

plt.title("Single carbon steel wall")
plt.xlabel("Position [m]")
plt.ylabel("Temperature [°C]")
plt.legend()
plt.show()
