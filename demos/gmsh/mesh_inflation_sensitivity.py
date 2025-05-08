import numpy as np
import matplotlib.pyplot as plt
import barotropy as bpy

# Parameters
h_wall = 1e-6     # m
r = 1.1           # growth ratio
H_total = 0.1      # m
N_total = 100     # total cells

# Residual functio
def residual(N_bl):
    H_bl = h_wall * (1 - r**N_bl) / (1 - r)
    h_bl = h_wall * r ** (N_bl - 1)
    N_core = N_total - N_bl
    H_core = H_total - H_bl
    h_core = H_core / N_core
    return h_core - h_bl

# Evaluate residual for integer values
N_vals = np.arange(1, N_total)
residuals = [residual(N) for N in N_vals]



# Plot
plt.figure()
plt.plot(N_vals, residuals, marker='o')
plt.xlabel("Number of inflation cells (N_inflation)")
plt.ylabel("Residual |h_core - h_last| [m]")
plt.title("Residual vs. number of inflation cells")
plt.grid(True)
# plt.yscale("log")
# plt.ylim([-0.01, 0.01])
plt.tight_layout()


bpy.find_optimal_inflation_layer(h_wall, r, H_total, N_total)

bpy.plot_y_mesh(h_wall, r, H_total, N_total)

plt.show()