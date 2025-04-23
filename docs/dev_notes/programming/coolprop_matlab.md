
# CoolProp-MATLAB installations instructions

## Installation with Conda
In order to use the code you need a MATLAB installation. Also, make sure that the `your_script.m` script includes a code-line to add the folder [`dependencies`](./dependencies) to the MATLAB path:
``` matlab
addpath(genpath("dependencies"))
```
In addition, you need to install the CoolProp fluid property library and interface it with MATLAB through Python. Check the step-by-step instructions below to learn how to interface MATLAB with CoolProp in a Windows operating system. The installation for Linux or Mac operating systems would be similar.

#### Step 1 - Download and install Miniconda
[Download the installer](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html) and follow the instructions until the installation is completed.

In order to use the Conda from from the [Git Bash](https://gitforwindows.org/) teminal, you have to initialize it in your terminal session. The easiest way to achieve this in a permanent way is to modify your `.bashrc` file to initialize conda each time that bash is started. Open your bash terminal and type the following command:
  ``` bash
  echo ". \${HOME}/AppData/Local/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc
  ```
This will add a new line to your `.bashrc` file that will execute `conda.sh` each time that you open the bash terminal. To check if the command worked, you can open your `.bashrc` file:
  ```bash
  notepad ~/.bashrc 
  ```
  
#### Step 2 - Create a virtual environment and install CoolProp
Open a new Git Bash terminal and type the following command to create a new virtual environment named `coolprop_env`:
```shell
conda create --name coolprop_env python=3.8
```
Now that the environment is created, you can activate it. Just use the following command:
```shell
conda activate coolprop_env
```
Finally, type the following command to install CoolProp:
```shell
conda install CoolProp --channel conda-forge
```
That's it! Note that it was necessary to tell Miniconda that it should look for Coolprop in the `conda-forge` channel.

#### Step 3 - Interface MATLAB and Coolprop
Open MATLAB (or close and open it if it was already open) and type the following command to let MATLAB know where is the Python executable (python.exe) of the virtual environment that you have just created:
```matlab
pyversion(fullfile(char(java.lang.System.getProperty('user.home')), '\AppData\Local\miniconda3\envs\coolprop_env\python.exe'))
```
This command will find the Python executable of the `coolprop_env` environment installed in your home directory.

Good! You have installed CoolProp and interfaced it with MATLAB. Let's do a simple test to check if the installation was successful. We are going to use CoolProp to compute the saturation temperature of water at atmospheric pressure. Just type the following command in MATLAB:
```matlab
py.CoolProp.CoolProp.PropsSI('T','P',101325,'Q',0,'Water')
```
If this does not throw and error and returns 373.1243 K, the installation was successful.



## More advanced use of Conda
If you are using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html), you can use the following [Bash](https://gitforwindows.org/) command to create a new virtual environment with all the dependencies required to run the code in this repository:
``` bash
conda env create --file environment.yaml
```
This will create the `barotropic_env` virtual environment and install all the packages in the specified in the YAML file.

To activate the virtual environment use:
``` bash
conda activate barotropic_env
```
If you need to install additional packages you can use the following command:
``` bash
conda install <name of the package>
```
You can also install new packages by adding their names to the `environment.yaml` file and updating the environment (using `--prune` removes any dependencies that are no longer required):
``` bash
conda env update --file environment.yaml --prune
```


