import pandas as pd
import matplotlib.pyplot as plt

# Compare data from simulation and experimental test


# Experimental Data
path = "C:/Users/alber/OneDrive - Danmarks Tekniske Universitet/Documents/CFD/barotropy/projects/Simoneau_Hendricks_1979/data_experiments/"
name = 'case_069_nitrogen_T124_10K_P39_20bar.csv'
case_number = 69
exp_res = path + name
data = pd.read_csv(exp_res, skiprows=1, header=None)
data.columns = ['Position', 'Pressure_Pa']
plt.figure(figsize=(8, 6))
plt.plot(data['Position'], data['Pressure_Pa'], marker='o', linestyle='-', color='k', label='Experimental Data')

#  Simulation data
path = "C:/Users/alber/OneDrive - Danmarks Tekniske Universitet/Documents/CFD/barotropy/projects/Simoneau_Hendricks_1979/automate_simulations/"
name = 'test_plot_pressure.csv'
exp_res = path + name
data = pd.read_csv(exp_res, skiprows=4, sep='\s+', header=None)
data.columns = ['Position', 'Pressure_Pa']
data['Position'] = pd.to_numeric(data['Position'], errors='coerce')
data['Pressure_Pa'] = pd.to_numeric(data['Pressure_Pa'], errors='coerce')
plt.plot(data['Position'], data['Pressure_Pa'], linestyle='-', color='b', label='CFD Data')


plt.xlabel('Position')
plt.ylabel('Pressure (Pa)')
plt.title(f"Case {case_number} with barotropic model")
plt.grid(True)
plt.legend()
plt.show()
