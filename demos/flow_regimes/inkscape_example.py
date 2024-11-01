
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import barotropy as bpy
from svgpathtools import svg2paths


def apply_conformal_mapping(x, y, x1, y1, r1, r2, c_ax, theta_0=0):
    """Applies conformal mapping to transform the coordinates (x, y) into a new set of coordinates."""

    # Convert (x, y) to (r, theta) coordinates preserving the wrap angle
    r = r1 * np.exp(np.log(r2 / r1) * (x - x1) / c_ax)
    theta = theta_0 + np.log(r2 / r1) / c_ax * (y - y1)
    x_new = r * np.cos(theta)
    y_new = r * np.sin(theta)
    
    return x_new, y_new


def get_coordinates_from_svg(svg_file):
    paths, attributes = svg2paths(svg_file)
    xy_points = []
    for path in paths:
        for segment in path:
            if segment.length() > 0:
                for i in np.linspace(0, 1, num=25):  # 100 points per segment
                    point = segment.point(i)
                    xy_points.append((point.real, point.imag))
    return xy_points



# Create the folder to save figures
bpy.set_plot_options(grid=False)
colors = bpy.COLORS_MATLAB
fig_dir = "output"
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# Define cascade parameters
r_1 = 2.2
r_2 = 3.0
r_3 = 3.2
r_4 = 4.2
N_stator = 15
N_rotor = 33

# Get stator coordinates
xy_stator = get_coordinates_from_svg('stator.svg')

x_stator = np.asarray([y for x, y in xy_stator])
y_stator = np.asarray([x for x, y in xy_stator])
c_stator = (np.max(x_stator) - np.min(x_stator))
x_stator = (x_stator - np.min(x_stator)) / c_stator
y_stator = (y_stator - np.min(y_stator)) / c_stator

# Get rotor blade coordinates
xy_rotor = get_coordinates_from_svg('rotor.svg')
x_rotor = np.asarray([x for x, y in xy_rotor])
y_rotor = np.asarray([-y for x, y in xy_rotor])
c_rotor = (np.max(x_rotor) - np.min(x_rotor))
x_rotor = (x_rotor - np.min(x_rotor)) / c_rotor
y_rotor = (y_rotor - np.min(y_rotor)) / c_rotor

# Create a figure with two subplots side by side
fig, ax = plt.subplots(figsize=(6, 4.2))
ax.set_aspect('equal')

# Plot stator blades
for theta in np.linspace(0, 2*np.pi, N_stator):
    # print(x_stator)
    # a = b
    x, y = apply_conformal_mapping(x_stator, y_stator, x1=0, y1=0, r1=r_1, r2=r_2, c_ax=1, theta_0=theta)
    ax.plot(x, y, color="black", linewidth=1)
    # ax.plot(x_stator, y_stator)

# Plot rotor blades
for theta in np.linspace(0, 2*np.pi, N_rotor):
    x, y = apply_conformal_mapping(x_rotor, y_rotor, x1=0, y1=0, r1=r_3, r2=r_4, c_ax=1, theta_0=theta)
    ax.plot(x, y, color="black", linewidth=1)
    # ax.plot(x_rotor, y_rotor)



theta = np.linspace(0, 2*np.pi, 200)
a, b = 1.00, 1.00
ax.plot(a*r_1*np.cos(theta), a*r_1*np.sin(theta), linestyle="--", color='black', linewidth=0.75)
ax.plot(b*r_2*np.cos(theta), b*r_2*np.sin(theta), linestyle="--", color='black', linewidth=0.75)
ax.plot(a*r_3*np.cos(theta), a*r_3*np.sin(theta), linestyle="--", color='black', linewidth=0.75)
ax.plot(b*r_4*np.cos(theta), b*r_4*np.sin(theta), linestyle="--", color='black', linewidth=0.75)


bpy.savefig_in_formats(fig, os.path.join(fig_dir, "radial_turbine"), formats=[".svg"])
plt.show()