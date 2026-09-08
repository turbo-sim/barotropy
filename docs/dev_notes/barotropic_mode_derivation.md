## Derivation of the Barotropic model equations for a single component system

The isentropic efficiency of an expansion process is defined as:
$$
\eta_{is} =  \dfrac{\Delta h}{\Delta h_{is}}
$$
Considering now the polytropic efficiency of an expansion process, this is defined as the efficiency of the efficiency of an expansion process composed of infinitesimal small stages. In other words, it is:
$$
 \eta_{pol} = \lim_{\Delta p \rightarrow 0} {\frac{\Delta h}{\Delta h_{is}}} =  \frac{d h}{d h_{is}}
$$
Considering the equation that relates the variation of entropy as a function of enthalpy and pressure:
$$
Tds = dh - \frac{1}{\rho}dp
$$
For an isentropic process we have that
$$
0 = dh_{is} - \dfrac{1}{\rho}dp \rightarrow dh_{is} =\dfrac{1}{\rho}dp
$$
Thus, substituting this result in the definition of polytropic efficiency becomes:
$$ 
\eta_{pol} = \rho \dfrac{dh}{dp} 
$$
With the final ODE of the barotropic model being:
$$
\dfrac{dh}{dp} = \dfrac{\eta_{pol}}{\rho}
$$
This latter equation allows to relate differences in pressure with differences in enthalpy, once the polytropic efficiency is known. With this strategy, it is possible to decouple the energy equation, simplifying the solution of the flow field.

As an example it is reported the calculation of the enthalpy at the step $n+1$ knowing all the variables at the step $n$
$$
\frac{h^{n+1} - h^n}{\Delta p} = \frac {\eta_{pol}}{\rho^n}
$$

$$
h^{n+1} = h^n + \Delta p \eta_{pol} \rho (h^n,p^n)
$$
where $p^{n+1} = p^n + \Delta p$ is the free variable.

