import matplotlib.pyplot as plt
import numpy as np
import nozzle_hyperbolic as hyp
import pandas as pd

hyp.set_plot_options(grid=False)

# Plot function1 for showing shock at position 25 50 75 100
# yaml file: euler_adapted_expansion


# Load the data from the files
data01 = np.load("Results/Case1_rusanov_N180.npz")
data= np.load("Results/trial90.npz")
# data02 = np.load("HLL_Case1_Nx115.npz")
data03 = np.load("Results/Case5_rusanov_N180.npz")
# data04 = np.load("HLL_Case5_Nx115.npz")

h, x1, u1, p1, rho1, Q1, a1 = data01["h"], data01["x"], data01["u"], data01["p"], data01["rho"], data01["Q"], data01["a"]
# x2, u2, p2, rho2, Q2, a2 = data02["x"], data02["u"], data02["p"], data02["rho"], data02["Q"], data02["a"]
x3, u3, p3, rho3, Q3, a3 = data03["x"], data03["u"], data03["p"], data03["rho"], data03["Q"], data03["a"]
# x4, u4, p4, rho4, Q4, a4 = data04["x"], data04["u"], data04["p"], data04["rho"], data04["Q"], data04["a"]
x, p = data["x"], data["p"]

# m1 = rho1[-1]*u1[-1]*h[-1]
# print(m1)
# m3 = rho3[-1]*u3[-1]*h[-1]
# print(m3)

exp_case1 = 'Case1.xlsx'
exp_case1 = pd.read_excel(exp_case1)

CFD_case1 = 'Case1_CFD.xlsx'
CFD_case1 = pd.read_excel(CFD_case1)

exp_case5 = 'Case5.xlsx'
exp_case5 = pd.read_excel(exp_case5)

CFD_case5 = 'Case5_CFD.xlsx'
CFD_case5 = pd.read_excel(CFD_case5)

# min_index = np.argmin(h)

## Case 1
x_nucleation_case1 = 3.748

plt.figure(figsize=(6, 5))
ax1 = plt.gca()
ax1.set_ylim([15, 60])
ax1.set_xlim([-35, 60])
plt.xlabel(r'$x$ (mm)')
plt.ylabel(r'$P$ (bar)')
plt.plot(x1*1e3, p1/1e5, '-', color='darkorange', linewidth=1.8, markerfacecolor='white', label="1D HEM")
# plt.plot(x*1e3, p/1e5, '-', color='purple', linewidth=1.8, markerfacecolor='white', label="Q_onset=0.992")
# plt.plot(x2, p2/1e5, '-', color='orange', linewidth=1.5, label="HLL")
plt.plot(CFD_case1['Position'], CFD_case1['Pressure'], '-', color='dimgray', linewidth=1.8, label='Romei et al. (2021)')
plt.plot(exp_case1['Position'], exp_case1['Pressure'], 's', color='black', markersize=6, markerfacecolor='k', label='Lettieri et al. (2017)')

# Add the second y-axis for the vapour quality area
plt.grid(False)
ax2 = ax1.twinx()
ax2.plot(x1*1e3, (1-Q1), '--', linewidth=1.8, color="darkcyan", label='Liquid mass fraction')
ax2.set_ylabel(r'Liquid mass fraction (-)')  # Adjust labelpad to move it away from the frame
ax2.set_ylim([0, 0.5])
y2_color = "darkcyan"
ax2.spines['right'].set_color(y2_color)  # Change the spine color
ax2.tick_params(axis='y', colors=y2_color)  # Change the tick color
ax2.yaxis.label.set_color(y2_color)  # Change the label color
plt.tight_layout()
plt.axvline(x=x_nucleation_case1, color='red', linestyle='-', linewidth=2, ymin=0, ymax=0.035, label="Experimental condensation onset")

# Get handles and labels from both axes
handles1, labels1 = ax1.get_legend_handles_labels()  # From primary axis
handles2, labels2 = ax2.get_legend_handles_labels()  # From secondary axis
# Combine handles and labels
handles = handles1 + handles2
labels = labels1 + labels2

# Create a combined legend with your custom settings
plt.legend(
    handles, labels,
    loc='upper right',
    bbox_to_anchor=(0.98, 0.98),  # Custom positioning
    fontsize='x-small',        # Custom font size
    handlelength=1,            # Shorten the line markers
    handletextpad=0.3,         # Reduce space between markers and text
    borderpad=0.3,             # Decrease padding inside the legend box
    frameon=True               # Keep the frame around the legend
)



## Case 5
x_nucleation_case5 = -4.675

plt.figure(figsize=(6, 5))
ax1 = plt.gca()
ax1.set_ylim([25, 85])
ax1.set_xlim([-40, 60])
plt.xlabel(r'$x$ (mm)')
plt.ylabel(r'$P$ (bar)')
plt.plot(x3*1e3, p3/1e5, '-', color='darkorange', linewidth=1.8, markerfacecolor='white', label="1D HEM")
# plt.plot(x4, p4/1e5, '-', color='orange', linewidth=1.5, label="HLL")
plt.plot(CFD_case5['Position'], CFD_case5['Pressure'], '-', color='dimgray', linewidth=1.8, label='Romei et al. (2021)')
plt.plot(exp_case5['Position'], exp_case5['Pressure'], 's', color='black', markersize=6, markerfacecolor='k', label='Lettieri et al. (2017)')

# Add the second y-axis for the vapour quality area
plt.grid(False)
ax2 = ax1.twinx()
ax2.plot(x3*1e3, (1-Q3), '--', linewidth=1.8, color="darkcyan", label='Liquid mass fraction ')
ax2.set_ylabel(r'Liquid mass fraction (-)')  # Adjust labelpad to move it away from the frame
ax2.set_ylim([0, 0.5])
y2_color = "darkcyan"
ax2.spines['right'].set_color(y2_color)  # Change the spine color
ax2.tick_params(axis='y', colors=y2_color)  # Change the tick color
ax2.yaxis.label.set_color(y2_color)  # Change the label color
plt.tight_layout()
plt.axvline(x=x_nucleation_case5, color='red', linestyle='-', linewidth=2, ymin=0, ymax=0.035, label="Experimental condensation onset")

# Get handles and labels from both axes
handles1, labels1 = ax1.get_legend_handles_labels()  # From primary axis
handles2, labels2 = ax2.get_legend_handles_labels()  # From secondary axis
# Combine handles and labels
handles = handles1 + handles2
labels = labels1 + labels2

# Create a combined legend with your custom settings
plt.legend(
    handles, labels,
    loc='upper right',
    bbox_to_anchor=(0.98, 0.98),  # Custom positioning
    fontsize='x-small',        # Custom font size
    handlelength=1,            # Shorten the line markers
    handletextpad=0.3,         # Reduce space between markers and text
    borderpad=0.3,             # Decrease padding inside the legend box
    frameon=True               # Keep the frame around the legend
)

plt.show()
