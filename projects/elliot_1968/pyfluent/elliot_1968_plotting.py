import os
import pandas as pd
import matplotlib.pyplot as plt

import barotropy as bpy

# TODO @Davoud: If I remember correctly, the agreement between the experimental and the simulated pressure distributions was better when you first showed it to me. Are we getting the same simulation results?

bpy.print_package_info()
bpy.set_plot_options()

# Create output directory
DIR_OUT = "output/validation"
case_name = "pressure_distribution_case_14"
os.makedirs(DIR_OUT, exist_ok=True)

# Create plot
fig, ax = plt.subplots(figsize=(6.4, 4.8))
# ax.axis("equal")
ax.set_xlabel('Distance [mm]')
ax.set_ylabel('Pressure [Pa]')

# Plot simulation data
data_sim = bpy.parse_fluent_xy("./output/case_14/latest/fluent_simulation/case_14_pressure.xy")
x0 = data_sim["centerline_conv"]["X-Coordinate"][0]
for segment in ["centerline_conv", "centerline_throat", "centerline_dive"]:
    x_sim = 1000*(data_sim[segment]["X-Coordinate"]-x0)
    y_sim = data_sim[segment]["Static Pressure"]
    ax.plot(x_sim, y_sim, linestyle="-",color=bpy.COLORS_MATLAB[0])

# Plot experimental data
data_exp = pd.read_excel("./pressure_distribution_case_14.xlsx")
data_exp["distance"] = data_exp["distance_in"] * 25.4  # Convert inches to mm
data_exp["pressure"] = data_exp["pressure_psi"] * 6894.76  # Convert psi to Pa
ax.plot(data_exp["distance"]-data_exp["distance"][0], data_exp["pressure"], marker="o", markerfacecolor="white", linestyle="none",color=bpy.COLORS_MATLAB[1])
plt.tight_layout(pad=1)

# Save figure
bpy.savefig_in_formats(fig, os.path.join(DIR_OUT, case_name))

# Show figure
plt.show()