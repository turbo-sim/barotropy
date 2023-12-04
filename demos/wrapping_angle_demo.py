import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp

# Meridional length (obtained from FreeCAD)
M = 0.1
m = np.linspace(0, M, 50)

# Metal angle distribution (random example)
beta_1 = -70
beta_2 = -30
beta = lambda m: beta_1 + (beta_2 - beta_1)/M*m

# Radius distribution (random example)
r_1 = 0.2
r_2 = 0.5
r = lambda m: r_1 + (r_2 - r_1)*(m/M)**2

# Radius distribution (random example)
z_1 = 0.0
z_2 = 0.2
z = lambda m: z_1 + (z_2 - z_1)*(m/M)**2

# Slope of the wrap angle distribution
def get_wrap_angle_slope(m, theta):
    return np.tan(np.deg2rad(beta(m)))/r(m)

# Compute the wrap angle distribution
theta_0 = 0.00
sol = solve_ivp(get_wrap_angle_slope, [0, M], [theta_0], t_eval=m, method='RK45', rtol=1e6, atol=1e-6)
theta = sol.y.flatten()

# Plot the data
fig, ax = plt.subplots()
ax.set_xlabel('Meridional length (m)')
ax.set_ylabel('Wrap angle distribution (degrees)')
ax.plot(m, np.rad2deg(theta), color='black', marker="o", markerfacecolor="white", markersize=3.0)
ax.legend(fontsize=10)
fig.tight_layout(pad=1, w_pad=None, h_pad=None)

# Show figures
plt.show()