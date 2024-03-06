import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

import barotropy as bpy
import barotropy.fluids.core as props
import barotropy.sCO2_utilities as sco2


# ---------------------------------------------------------------------------- #
# Load and postprocess data
# ---------------------------------------------------------------------------- #
# FILENAME = "sCO2_experimental_data.xlsx"
OUTPUT_DIR = "figures"
SAVE_FIGS = True
bpy.set_plot_options(minor_ticks=False)
Nx = 20
Ny = 20


# # Define plot settings
# bpy.set_plot_options(minor_ticks=False, grid=False)
# markers = ["P", "v", "s", "^", "o", "X", '<']
# marker_sizes = {'o': 3.5, 'v': 4.5, 's': 3.5, '^': 4.5, 'P': 4.5, 'X': 4.5, '<': 3.5}

# n = len(df["compressor"].unique())
# colors = plt.get_cmap("magma")(np.linspace(0.80, 0.3, n))

# Create fluid object
fluid = props.Fluid(name="CO2", backend="HEOS", exceptions=True)



T = 300
p = 75e5
state = fluid.set_state(props.PT_INPUTS, p, T)

a = props.set_state_Qs(fluid, 0.1, state.s)


fluid.triple_point_liquid
print(a)


# state = fluid.set_state(props.PT_INPUTS, p, T)

# state = fluid.set_state(props.QSmass_INPUTS, 0.00, state.s)

# print(state)

# for key in state.to_dict():
#     print(key)

# for key in state:
#     print(key)

# print(state)
# print(state.fluid_name)

# # print(state._properties)





# # fl

# # ---------------------------------------------------------------------------- #
# # Plot inlet states on p-s diagram
# # ---------------------------------------------------------------------------- #
# # Create figure
# fig, ax = plt.subplots(figsize=(6.8, 4.8))
# ax.set_xlabel("Entropy [kJ/kg/K]")
# ax.set_ylabel("Pressure [bar]")
# range_x = np.linspace(1000, 2000, Nx)
# range_y = np.linspace(59, 120, Ny)*1e5
# range_z = np.asarray([30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]) + 273.15
# ax.set_xlim(np.asarray([range_x[0], range_x[-1]]) / 1e3)
# ax.set_ylim(np.asarray([range_y[0], range_y[-1]]) / 1e5)

# # Plot diagram
# ax = sco2.plot_ps_diagram(fluid, ax, range_x, range_y, range_z)


# T = 300
# p = 75e5
# state = fluid.set_state(props.PT_INPUTS, p, T)
# ax.plot(state.s, state.p, 'ro')
# print(state)


# state = fluid.set_state(props.QSmass_INPUTS, 0.00, state.s)
# ax.plot(state.s, state.p, 'ro')
# print(state)


# # Place the legend outside of the figure
# # ax.legend(loc="upper right", fontsize=10)
# fig.tight_layout(pad=1)

# # Scale to nice units
# bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
# bpy.scale_graphics_y(fig, +1e-5, mode="multiply")

# # # Save figure
# # if SAVE_FIGS:
# #     base_filename = os.path.join(OUTPUT_DIR, f"testing")
# #     bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])



# # # ---------------------------------------------------------------------------- #
# # # Plot inlet states on T-s diagram
# # # ---------------------------------------------------------------------------- #
# # # Create figure
# # fig, ax = plt.subplots(figsize=(6.8, 4.8))
# # ax.set_xlabel("Entropy [kJ/kg/K]")
# # ax.set_ylabel("Temperature [$^\circ$C]")
# # range_x = np.linspace(1000, 2000, Nx)
# # range_y = np.linspace(5, 55, Ny) + 273.15
# # range_z = np.asarray([50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]) * 1e5
# # ax.set_xlim(np.asarray([range_x[0], range_x[-1]]) / 1e3)
# # ax.set_ylim(np.asarray([range_y[0], range_y[-1]]) - 273.15)

# # # Plot diagram
# # ax = sco2.plot_Ts_diagram(fluid, ax, range_x, range_y, range_z)

# # # Place the legend outside of the figure
# # ax.legend(loc="lower right", fontsize=10)

# # # Scale to nice units
# # bpy.scale_graphics_x(fig, +1e-3, mode="multiply")
# # bpy.scale_graphics_y(fig, -273.15, mode="add")

# # # # Save figure
# # # if SAVE_FIGS:
# # #     base_filename = os.path.join(OUTPUT_DIR, f"sCO2_experimental_data_Ts_diagram")
# # #     bpy.savefig_in_formats(fig, base_filename, formats=[".png", ".svg"])




# # Show figures
# plt.show()