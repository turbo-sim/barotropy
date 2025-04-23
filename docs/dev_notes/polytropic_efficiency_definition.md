## General real-gas definition
The polytropic efficiency can be defied as the limit of the isentropic efficiency when the enthalpy difference tends to zero (i.e., very small pressure ratio):
$$\eta_\text{p} = \lim_{\Delta h \to 0} \frac{\Delta h_s}{\Delta h} = \frac{\text{d} h_s}{\text{d} h}$$
This expression can be reformulated using the Gibbs relation:
$$T\text{d}s = \text{d}h - \text{d}p/\rho \to \text{d}h_s =  \text{d}p/\rho = \text{d}h- T\text{d}s$$
to yield:
$$ \eta_\text{p} =  \frac{\text{d} h_s}{\text{d} h} = \frac{1}{\rho} \frac{\text{d}p}{\text{d}h} = 1 - T \frac{\text{d}s}{\text{d}h}$$
The states along the polytropic process can be calculated integrating either of the following two ordinary differential equations:
$$ \begin{gather}
\frac{\text{d}h}{\text{d}p} = \frac{1}{\rho \,\eta_\text{p}} \quad \text{Using pressure-enthalpy function calls}\\
\frac{\text{d}s}{\text{d}h} = \frac{1}{T}(1-\eta_\text{p})  \quad \text{Using enthalpy-entropy function calls} 
\end{gather}$$


## Approximate real-gas definition (Mallen and Saville)

When the enthalpy is proportional to temperature (as when the fluid behaves as an ideal gas), the polytropic process can be solved analytically:
$$ \begin{gather}
\text{d}s = (1-\eta_\text{p}) \frac{\text{d}h}{T} = (1-\eta_\text{p}) c_p\frac{\text{d}T}{T} \\
s_2 - s_1 = (1-\eta_\text{p}) c_p \ln{\left(T_2/T_1\right)} \\
s_2 - s_1 = (1-\eta_\text{p}) c_p (T_2-T_1)\frac{\ln{\left(T_2/T_1\right)}}{(T_2 - T_1)} \\
\hat{T}(s_2 - s_1) = (1-\eta_\text{p}) (h_2-h_1) \\
\hat{T}\Delta s = \Delta h - \Delta h_p \\
\end{gather}$$
where $\hat{T} = (T_2-T_1)/\ln{(T_2/T_1)}$ is the mean logarithmic temperature between inlet and outlet states. This expression can be reformulated rewritten as follows:
$$
\begin{gather}
\eta_\text{p} = 1- \frac{\Delta s}{\Delta h} \hat{T} = 1 - \left(\frac{s_2-s_1}{h_2-h_1}\right)  \left(\frac{T_2-T_1}{\ln{(T_2/T_1)}}\right)
\end{gather}
$$
When the inlet state and polytropic efficiency are known, the outlet state can be solved iteratively to satisfy the equation above. 

## Ideal gas polytropic process
In the case of a perfect gas, this equation can be solved explicitly using the definition for the enthalpy and entropy changes:
$$
\begin{gather}
h_2 - h_1 = c_p(T_2- T_1) \\
s_2 - s_1 = c_p \ln(T_2/T_1) - R \ln(p_2/p_1)
\end{gather}
$$
leading to the following expression:
$$
\eta_\text{p} = \frac{\gamma - 1 }{\gamma} \frac{\ln(p_2/p_1)}{\ln(T_2/T_1)}
$$
which can be re-arranged as:
$$
\left( \frac{T_2}{T_1} \right) = \left( \frac{p_2}{p_1} \right)^{\frac{1}{\eta_\text{p}} \left(\frac{\gamma-1}{\gamma}\right)}
$$
This same expression can also be derived by replacing the ideal gas equation into the polytropic efficiency definition:
$$
\begin{gather}
\eta_\text{p} =  \frac{\text{d} h_s}{\text{d} h} = \frac{1}{\rho\,c_p} \frac{\text{d}p}{\text{d}T} = \frac{RT}{p\,c_p} \frac{\text{d}p}{\text{d}T} \\
\eta_\text{p} c_p \frac{\text{d}T}{T} = R\frac{\text{d}p}{p} \\
\frac{\eta_\text{p} \, \gamma }{\gamma -1} \frac{\text{d}T}{T} = \frac{\text{d}p}{p} \\  
\left( \frac{T_2}{T_1} \right) = \left( \frac{p_2}{p_1} \right)^{\frac{1}{\eta_\text{p}} \left(\frac{\gamma-1}{\gamma}\right)}
\end{gather}
$$


## Comparison between isentropic and polytropic efficiency

**Actual enthalpy difference:**
$$
h_2 - h_1 = c_p(T_2-T_1) = \gamma R T_1/(\gamma-1) \left[ \left( \frac{p_2}{p_1} \right)^{\frac{1}{\eta_\text{p}} \left(\frac{\gamma-1}{\gamma}\right)} - 1 \right]
$$
**Isentropic enthalpy difference:**
$$
h_{2s} - h_1 = c_p(T_{2s}-T_1) = \gamma R T_1/(\gamma-1) \left[ \left( \frac{p_2}{p_1} \right)^{\left(\frac{\gamma-1}{\gamma}\right)} - 1 \right]
$$
Therefore, the isentropic and polytropic efficiencies are related according to the following relation:
$$
\eta_s = \frac{h_{2s} - h_1}{h_2 - h_1} = \frac{  \left( \frac{p_2}{p_1} \right)^{\left(\frac{\gamma-1}{\gamma}\right)} - 1 }{\left( \frac{p_2}{p_1} \right)^{\frac{1}{\eta_\text{p}} \left(\frac{\gamma-1}{\gamma}\right)} - 1}
$$
where for the case of a compression, it holds that $\eta_s < \eta_{\text{p}}$.



## Polytropic efficiency in a two-phase expansion process


The polytropic efficiency is be defied as the limit of the isentropic efficiency when the enthalpy difference tends to zero (i.e., very small pressure ratio):


$$\eta_\text{p} = \lim_{\Delta h_s \to 0} \frac{\Delta h}{\Delta h_s} = \frac{\text{d} h}{\text{d} h_s} =  \rho  \frac{\text{d}h}{\text{d}p}$$
where $h$ is the specific enthalpy of the mixture and $h_s$ is the enthalpy along an ideal isentropic process. Assuming thermodynamic equilibrium between phases, the specific enthalpy of the two-phase mixture is given by:
$$
h = y_\text{g} \,h_\text{g}(p, \, T) + y_\text{l}\, h_\text{l}(p, \, T)
$$
Taking the total derivative of enthalpy with respect to pressure and applying the polytropic relation yields:
$$
\begin{gather}
    \frac{\text{d}T}{\text{d}p} = \frac{1}{\rho c_p} (\rho \mu_T + \eta_p)
\end{gather}
$$
where mixture density $\rho$, specific heat capacity $c_p$, and Joule–Thomson coefficient $\mu_T$ are evaluated from the mass fractions and phase properties:
$$
\begin{gather}
\frac{1}{\rho} = \frac{y_\text{g}}{\rho_\text{g}} +  \frac{y_\text{l}}{\rho_\text{l}} \\
c_p = \left(\frac{\partial h}{\partial T}\right)_p =  y_\text{g} \, c_{p,\text{g}} + y_\text{l} \, c_{p,\text{l}} \\
\mu_T = \left(\frac{\partial h}{\partial p}\right)_T =  y_\text{g} \, \mu_{T,\text{g}} + y_\text{l} \, \mu_{T,\text{l}}
\end{gather}
$$
The evolution of temperature as a function of pressure is obtained by integrating Eq. XX from the inlet pressure and temperature to the outlet pressure. Once the pressure–temperature trajectory along the expansion is known, other thermodynamic properties of the two-phase mixture are computed as weighted averages of the phase properties. For example, the dynamic viscosity $\mu$ and the speed of sound $c$ are given by:
$$
\begin{gather}
    \mu = \phi_\text{g} \,\mu_\text{g}(p,T) + \phi_\text{l} \,\mu_\text{l}(p,T) \\
    \frac{1}{\rho c^2} = \frac{\phi_\text{g}}{\rho_\text{g} c_\text{g}^2(p,T)} + \frac{\phi_\text{l}}{\rho_\text{l} c_\text{l}^2(p,T)}
\end{gather}
$$
where the volumetric fractions of the gas and liquid phases, $\phi_\text{g}$ and $\phi_\text{l}$, are defined as:

#### Introducing the compressibility factors

The isothermal Joule-Thomson coefficient of the mixture can expressed as:

$$  
\mu_T = \frac{1}{\rho}(1- \beta T)
$$
where $\beta = -\frac{1}{\rho} \left(\frac{\partial \rho}{\partial T}\right)_p$ is the coefficient of isobaric expansion. Introducing this parameter into the polytropic differential equation for temperature we get:

$$
\begin{gather}
    \frac{\text{d}T}{\text{d}p} = \frac{1}{\rho c_p} (1 + \eta_p - \beta T)
\end{gather}
$$
which for a perfect gas ($\beta=1/T$) yields:
$$
\begin{gather}
    \frac{\text{d}T}{\text{d}p} = \frac{\eta_p}{\rho c_p} 
\end{gather}
$$

#### Explaining the change in density

I have not rigorously proven this, but my intuition tells me that:

$$
\begin{gather}
    \text{d}\rho = \rho \kappa \text{d}p + \rho\beta \text{d}T    \\
    \frac{1}{\rho}\frac{\text{d}\rho}{\text{d}p} = \kappa + \frac{\beta}{\rho c_p}(1 + \eta_p - \beta T)
    
\end{gather}
$$

Maybe it is possible to get the approximate compressibility factors when the liquid is incompressible and learn some things analytically. However, it seems that it is not critical for the time being and just a numeric sensitivity analysis is sufficient to show that there is no sensitivity.