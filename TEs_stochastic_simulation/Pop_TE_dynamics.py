import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

####### Constants for the ODE
alpha = 2.0 / 10
beta1 = 1.0 / 100
beta2 = 0.5 / 100
delta1 = 0.5 / 10
delta2 = 1.5 / 100
V = 500  # Volume for normalization

####### Define the hex color for Bittersweet
bittersweet = '#C04F17'

####### ODE system
def odes(t, z):
    x, y = z
    dxdt = x * ((beta1 * alpha) / (beta1 + delta2 + beta2 * y) - delta1)
    dydt = y * ((beta2 * alpha * x) / (beta1 + delta2 + beta2 * y) - delta1)
    return [dxdt, dydt]

####### Time span and initial conditions
t_span = (0, 800)
t_eval = np.linspace(0, 800, 500)

initial_conditions = [0.4, 0.4]

####### Solve the ODE
sol = solve_ivp(odes, t_span, initial_conditions, t_eval=t_eval)

####### ODE solution
x_values = sol.y[0]
y_values = sol.y[1]

####### Read the CSV file
file_path = "/paths/pop_TE_dynamics.csv"

data = pd.read_csv(file_path)

####### Normalize by V = 500
data['X_scaled'] = data['X'] / 500
data['Y_scaled'] = data['Y'] / 500

plt.figure(figsize=(10, 6.5))

######## Plot deterministic ODE results
plt.plot(sol.t, x_values, color='blue', linestyle='-', linewidth=2)
plt.plot(sol.t, y_values, color='green', linestyle='-', linewidth=2)

######## Plot scaled stochastic data
plt.plot(data['generation'], data['X_scaled'], color='blue', linestyle='-', linewidth=0.8)
plt.plot(data['generation'], data['Y_scaled'], color='green', linestyle='-', linewidth=0.8)

######## Set font size and make tick labels bold
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')

######## Add labels and title
# plt.xlabel('Time (generation)')
# plt.ylabel('Copy number (scaled)')

######## Gridline visibility
plt.grid(True, alpha=0.1)

######## Show the plot
plt.show()
