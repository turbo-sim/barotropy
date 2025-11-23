
## Derivation of the barotropic model


The barotropic model for two components is based on the following assumptions:

- Mechanical equilibrium: 
  $$p_1=p_2=p$$
- Thermal equilibrium
  $$T_1=T_2=T$$
- No mass transfer between phases
  $$y_1=\mathrm{constant} \quad \text{and} \quad y_2=\mathrm{constant}$$

- The process is characterized by a polytropic efficiency defined as the ratio of the actual enthalpy change to the ideal (or isentropic) enthalpy change over an infinitesimal small segment of the process:

    $$\frac{\mathrm{d}h}{\mathrm{d}h_s} = \eta_p$$

    where:
    - $h$ is the enthalpy of the fluid for the actual process
    - $h_s$ is the enthalpy of the fluid for an ideal isentropic process
    - $\eta_p$ is the polytropic efficiency of the machine



From the fundamental $T\mathrm{d}s$ equation we know that for an isentropic process:

$$\mathrm{d}h_s = \frac{\mathrm{d}p}{\rho}$$

such that the polytropic relation can be expressed as:

$$\frac{\mathrm{d}h}{\mathrm{d}p} = \frac{\eta_p}{\rho}$$


This equation relates the enthalpy of the mixture with the pressure and densities of the mixture.

In order to evaluate the right hand side, we need the density of the mixture, which is given by:

$$\frac{1}{\rho} = \frac{y_1}{\rho_1(p,T)} +  \frac{y_2}{\rho_2(p,T)}$$

This means that we need the temperature and pressure to evaluate the density of the mixture and integrate the ODE posed by the polytropic relation.

We can derive an expression to compute the temperature from the definition of the enthalpy of the mixture

$$h = y_1 \, h_1(p,T) + y_2 \, h_2(p,T)$$

The differential of enthalpy can be expressed as:


$$ \mathrm{d}h = 
y_1 \, \left[\left(\frac{\partial h_1}{\partial T}\right)_T  \mathrm{d}T + \left(\frac{\partial h_1}{\partial p}\right)_T  \mathrm{d}p \right] +
y_2 \, \left[\left(\frac{\partial h_2}{\partial T}\right)_T  \mathrm{d}T + \left(\frac{\partial h_2}{\partial p}\right)_T  \mathrm{d}p \right] 
$$

$$ \mathrm{d}h = 
y_1 \, \left(c_{p,1} \mathrm{d}T + \mu_{T,1} \mathrm{d}p \right) +
y_2 \, \left(c_{p,2} \mathrm{d}T + \mu_{T,2} \mathrm{d}p \right)
$$

$$ \mathrm{d}h = 
\left( y_1 \, c_{p,1} +y_2 \, c_{p,2} \right) \mathrm{d}T + 
\left( y_1 \, \mu_{p,1} +y_2 \, \mu_{T,2} \right) \mathrm{d}p
$$

$$ \mathrm{d}h = c_p \, \mathrm{d}T + \mu_{T} \, \mathrm{d}p
$$

with:

$$
c_p = y_1 \, c_{p,1} +y_2 \, c_{p,2} \\
\mu_{T} = y_1 \, \mu_{p,1} +y_2 \, \mu_{T,2} 
$$

Therefore, the total derivative of the mixture enthalpy with respect to pressur ecan be expressed as:

$$\frac{\mathrm{d}h}{\mathrm{d}p} = c_p \frac{\mathrm{d}T}{\mathrm{d}p} + \mu_{T}$$

solving for the derivative of temperature we get:

$$
\frac{\mathrm{d}T}{\mathrm{d}p} = \frac{1}{c_p} \left(\frac{\mathrm{d}h}{\mathrm{d}p} - \mu_T \right)
$$


replacing the derivative of enthalpy with respect to pressure from the polytropic equation we find the following system of ordinary differential equations:


$$
\begin{align}
\frac{\mathrm{d}h}{\mathrm{d}p} &= \frac{\eta_p}{\rho} \\
\frac{\mathrm{d}T}{\mathrm{d}p} &= \frac{1}{c_p} \left(\frac{\eta_p}{\rho} - \mu_T \right)
\end{align}
$$

These equations can be integrated from the initial condition numerically.


The properties of the nitrogen-water mixture were plotted as a function of pressure for the different experimental mixture ratios.


### Mixture composition
The mixture ratio is defined as the quotient between the mass fraction of water to nitrogen:

$$R = y_\mathrm{1}/y_\mathrm{2}$$

where subscripts 1 and 2 denote water and nitrogen, respectively.

The mass fraction of each component can be expressed as a function of the mixture ratio
$$y_1 = \frac{R}{1+R}$$
$$y_2 = \frac{1}{1+R}$$


### Thermodynamic property calculation
The thermodynamic properties of the mixture as a function of pressure are computed solving the following system of ordinary differential equations

$$
\begin{align}
\frac{\mathrm{d}h}{\mathrm{d}p} &= \frac{\eta_p}{\rho} \\
\frac{\mathrm{d}T}{\mathrm{d}p} &= \frac{1}{c_p} \left(\frac{\eta_p}{\rho} - \mu_T \right)
\end{align}
$$

with initial conditions:
The initial conditions of the problem are given by:

$$
\begin{gather}
p_0 = p_\mathrm{in} \\
T_0 = T_\mathrm{in} \\
h_0 = y_1 \, h_1(p_\mathrm{in},T_\mathrm{in}) + y_2 \, h_2(p_\mathrm{in},T_\mathrm{in})  \\
\end{gather}
$$

The thermophysical properties appearing in the ODE system are:
- $h$ - specific enthalpy of the mixture
- $T$ - temperature (same for both phases)
- $p$ - pressure (same for both phases)
- $\rho$ - density of the mixture
- $c_p=\left(\frac{\mathrm{\partial}h}{\mathrm{\partial}T}\right)_p$ - isobaric heat capacity of the mixture
- $\mu_T=\left(\frac{\mathrm{\partial}h}{\mathrm{\partial}p}\right)_T$ - isothermal Joule-Thompson coefficient of the mixture

In addition, $\eta_p$ is the polytropic efficiency of the process, with:
  - $\eta_p = 1$ for an isentropic process
  - $\eta_p = 0$ for an isenthalpic process


In order to integrate the ODE system numerically we have to evaluate the right hand side of the system. This involves the evaluation of the mixture density, isobaric heat capacity, and derivative of enthalpy with respect to pressure at constant enthalpy.

#### Evaluation of density

The density of the mixture is given by the mass-weighted harmonic average of the density of each component:

$$\frac{1}{\rho} = \frac{y_1}{\rho_1(p,T)} +  \frac{y_2}{\rho_2(p,T)}$$



$$-\frac{1}{\rho^2} \left(\frac{\partial \rho}{\partial p}\right)_s = -\frac{y_1}{\rho_1^2} \left(\frac{\partial \rho_1}{\partial p}\right)_{s_1} -\frac{y_2}{\rho_2^2} \left(\frac{\partial \rho_2}{\partial p}\right)_{s_2} $$

$$
y_1 / \rho_1 = \alpha_1 / \rho \qquad y_2 / \rho_2 = \alpha_2 / \rho
$$

$$\frac{1}{\rho^2} \left(\frac{\partial \rho}{\partial p}\right)_s = \frac{\alpha_1}{\rho \rho_1} \left(\frac{\partial \rho_1}{\partial p}\right)_{s_1} + \frac{\alpha_2}{\rho \rho_2} \left(\frac{\partial \rho_2}{\partial p}\right)_{s_2} $$

$$\frac{1}{\rho} \left(\frac{\partial \rho}{\partial p}\right)_s = \frac{\alpha_1}{\rho_1} \left(\frac{\partial \rho_1}{\partial p}\right)_{s_1} + \frac{\alpha_2}{\rho_2} \left(\frac{\partial \rho_2}{\partial p}\right)_{s_2} $$

$$\frac{1}{\rho c^2} = \frac{\alpha_1}{\rho_1 c_1^2} + \frac{\alpha_2}{\rho_2 c_2^2} $$


$$
\rho = \rho (p, T(s, p)) ????
$$


#### Isobaric heat capacity

The isobaric heat capacity of the mixture is given by the mass-weighted arithmetic average of the heat capacity of each component:
$$c_p = y_1 \, c_{p,1}(p,T) + y_2 \, c_{p,2}(p,T)$$



#### Isothermal Joule-Thompson coefficient
The isothermal Joule-Thompson coefficient of the mixture is given by the mass-weighted arithmetic average of the coefficient of each component
$$\mu_T = y_1 \, \mu_{T,1}(p,T) + y_2 \, \mu_{T,2}(p,T)$$

The isothermal Joule-Thompson coefficient of each component can be calculated as:

$$\mu_{T,i} = \left(\frac{\mathrm{\partial}h_i}{\mathrm{\partial}p}\right)_T = \frac{1}{\rho_i}(1-\alpha_i\,T)$$

where $\alpha_i = -\frac{1}{\rho_i} \left(\frac{\partial \rho_i}{\partial T}\right)_p$ is the isobaric expansion coefficient of component $i$.


## Solving the ODE

Numerical integration with Runge-Kutta methods


Once the right hand side of the equations is evaluated, the system of equations can be integrated using a numerical algorithm such as ODE45




## Postprocessing

Once the system of ODEs is solved, the enthalpy and temperature of the mixture are known as a function of pressure. Other thermodynamic properties such as void fraction, speed of sound, and entropy are calculated from the ODE solution in a postprocessing step.

#### Entropy

The entropy of the mixture is evaluated as:
$$s = y_1 \, s_1(p,T) + y_2 \, s_2(p,T)$$

The entropy of the mixture should be constant for a process with $\eta_p=1$ this serves as a verification step.

#### Enthalpy
Similarly, the enthalpy of the mixture can be re-calculated as:

$$h = y_1 \, h_1(p,T) + y_2 \, h_2(p,T)$$

This evaluation of the enthalpy serves as a verification that the enthalpy obtained from the the numerical integration of the ODE system does not deviate from its definition due to, for example, numerical tolerances or discretization error of the numerical ODE solver.

The enthalpy of the mixture should be constant for a process with $\eta_p=0$ this serves as a verification step.

#### Volume fractions
The volume fractions of each component can be calculated from the density of the mixture and the density of each component as:

$$\phi_1 = \left(\frac{\rho}{\rho_1}\right) y_1$$

$$\phi_2 = \left(\frac{\rho}{\rho_1}\right) y_2$$


#### Viscosity

The viscosity of the mixture is given by the volume-weighted arithmetic average of the viscosity of each component:

$$\mu = \phi_1 \, \mu_1(p,T) + \phi_2 \, \mu_2(p,T)$$

#### Speed of sound

The speed of sound of the mixture is evaluated by differentiation of the density equation at constant entropy.

The final relation indicates that the bulk modulus of the mixture $\rho c^2$ is given by the volume-weighted average of the bulk-moduli of each component:

$$\frac{1}{\rho c^2} = \frac{\phi_1}{\rho_1 c_1^2} +  \frac{\phi_2}{\rho_2 c_2^2}$$


This expression is the well known Wood's formula for the speed of sound under the assumption of Frozen eqauilibrium, also referred to frozen equilibrium model (HFM)






## Verification of the model

Check 1: Verify that the enthalpy of the mixture obtained by numerical integration of the ODE system agrees with the enthalpy of the mixture calculated from the mixing rule

Check 2: verify that entropy is constant for eta=1

Check 3: verify that the enthalpy is constant for eta=0


It seems the equation sare not properly satisfied.
The error does not seem numerical.
Try for the barotropic model for a single component to see if I can identify the error in a simpler case 




## Minimal working example

The following script summarizes how to integrate the system of equations in MATLAB:

```matlab
%% Barotropic fluid model calculation minimal working example
% Clear the workspace
clear all
close all
clc

% Define inlet and outlet conditions
T_in = 22.00 + 273.15;
p_in = 2000.0e3;
p_out = 98.6e3;

% Define the mixture composition    
R = 50.00;
y_1 = R./(1 + R);
y_2 = 1./(1 + R);

% Define polytropic efficiency of the process
eta_poly = 1.00;  % Isentropic process

% Initialize the fluid objects
py.importlib.import_module('CoolProp.CoolProp');
fluid_1 = py.CoolProp.CoolProp.AbstractState('HEOS', 'Water');
fluid_2 = py.CoolProp.CoolProp.AbstractState('HEOS', 'Nitrogen');

% Calculate the initial condition for enthalpy
fluid_1.update(py.CoolProp.CoolProp.PT_INPUTS, p_in, T_in);
fluid_2.update(py.CoolProp.CoolProp.PT_INPUTS, p_in, T_in);
h_in = y_1*fluid_1.hmass + y_2*fluid_2.hmass;

% Solve the ODE system
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
ode_handle = @(p, grad_hT) polytropic_expansion_ode(p, grad_hT, y_1, y_2, fluid_1, fluid_2, eta_poly);
[p, hT] = ode45(ode_handle, [p_in, p_out], [h_in, T_in], options);
T = hT(:, 2);

% Plot evolution of temperature
figure(); hold on; box on;
xlabel('Pressure (bar)')
ylabel('Temperature ($^{\circ}$C)')
plot(p/1e5, T-273.15, 'k-')

% Other thermodynamic properties can be evaluated from pressure,
% temperature and mixture composition as a post-processing step after the
% ODE system has been solved...

function grad_hT = polytropic_expansion_ode(p, hT, y_1, y_2, fluid_1, fluid_2, eta_poly)
    
    % Rename variables
    h = hT(1);
    T = hT(2);

    % Update thermodynamic state
    fluid_1.update(py.CoolProp.CoolProp.PT_INPUTS, p, T);
    fluid_2.update(py.CoolProp.CoolProp.PT_INPUTS, p, T);
    
    % Compute additional thermodynamic properties
    cp = y_1*fluid_1.cpmass + y_2*fluid_2.cpmass;
    rho = 1/(y_1/fluid_1.rhomass + y_2/fluid_2.rhomass);
    dhdp_T_1 = (1 - T*fluid_1.isobaric_expansion_coefficient)/fluid_1.rhomass;
    dhdp_T_2 = (1 - T*fluid_2.isobaric_expansion_coefficient)/fluid_2.rhomass;
    dhdp_T = y_1*dhdp_T_1 - y_2*dhdp_T_2;

    % Compute the slope of the polytropic process
    dhdp = eta_poly/rho;
    dTdp = (dhdp - dhdp_T)/cp;
    grad_hT = [dhdp; dTdp];  % Right hand side of ODE system

end
```



### Discussion of the property calculation Elliot-Case

It was checked that the numerical value of the polytropic efficiency had a negligible effect on all thermodynamic properties (except enthalpy and entropy):
- Density
- Viscosity
- Void fraction
- Speed of sound
- Temperature
- Enthalpy
- Entropy

Since we are interested in the properties that affect the barotropic model in the context of CFD simulations, the properties were evaluated at a polytropic efficiency of 100% (i.e., isentropic process).





<img src="./property_plots/water_nitrogen_mixture_density.png" alt="" width="500"/>


<img src="./property_plots/water_nitrogen_mixture_viscosity.png" alt="" width="500"/>


<img src="./property_plots/water_nitrogen_mixture_temperature.png" alt="" width="500"/>


<img src="./property_plots/water_nitrogen_mixture_void_fraction.png" alt="" width="500"/>


<img src="./property_plots/water_nitrogen_mixture_speed_sound.png" alt="" width="500"/>


<img src="./property_plots/water_nitrogen_mixture_enthalpy.png" alt="" width="500"/>

<img src="./property_plots/water_nitrogen_mixture_entropy.png" alt="" width="500"/>






## Density problem
- Problem exists with inviscid simulations as well
- Problem is independent of the density relaxation factor

I have to try low order polynomial

I changed from 8th order to 4th order and the deviation changed from +4% to +2%. This is an improvement, but it is still not clear why I am getting a deviation
The deviation is most significant at the exit because there is a shock wave at the exit
The problem is not the order of the polynomial leading to numerical error due to bad conditioning. I tried with order 8 and 12 and the results were essentially the same.
The order of the polynomials is not introducing errors in the density calculation.



I have to try UDF




## Density variation Fluent

Important:
If you plan to define density using a UDF, note that the solution convergence will become
poor as the density variation becomes large. Specifying a compressible law (density as a
function of pressure) or multiphase behavior (spatially varying density) may lead to divergence.
It is recommended that you restrict the use of UDFs for density to weakly compressible
flows with mild density variations.

## Speed of osund Fluent

Liquid density is not a constant but is instead a function of the pressure field. In order to stabilize
the pressure solution for compressible flows in Ansys Fluent, an extra term related to the speed of
sound is needed in the pressure correction equation. Consequently, when you want to define a
custom density function for a compressible flow, your model must also include a speed of sound
function. Although you can direct Ansys Fluent to calculate a speed of sound function by choosing
one of the available methods (for example, piecewise-linear, polynomial) in the Create/Edit Mater-

ials dialog box, as a general guideline you should define a speed of sound function along with
your density UDF using the formulation:
For simplicity, it is recommended that you concatenate the density and speed of sound functions
into a single UDF source file.
The following UDF source code example contains two concatenated functions: a density function
named superfluid_density that is defined in terms of pressure and a custom speed of sound
function named sound_speed.


``` cpp
/********************************************************************
Density and speed of sound UDFs.
*********************************************************************/
#include "udf.h"
#define BMODULUS 2.2e9
#define rho_ref 1000.0
#define p_ref 101325

DEFINE_PROPERTY(superfluid_density, c, t)
{
  real rho;
  real p, dp;
  p = C_P(c,t) + op_pres;
  dp = p-p_ref;
  rho = rho_ref/(1.0-dp/BMODULUS);
  return rho;
}
DEFINE_PROPERTY(sound_speed, c,t)
{
  real a;
  real p, dp;
  p = C_P(c,t) + op_pres;
  dp = p-p_ref; a = sqrt(BMODULUS*(1.-dp/BMODULUS)/rho_ref);
  return a;
}
```