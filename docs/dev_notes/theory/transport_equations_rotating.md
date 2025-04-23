## Mass and momentum conservation equations in a rotating frame (CFX formulation)

### Mass conservation
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0
$$

### Momentum conservation
$$
\frac{\partial (\rho \vec{u})}{\partial t} + \nabla \cdot (\rho \vec{u} \otimes \vec{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} - 2\rho \vec{\Omega} \times \vec{u} - \rho \vec{\Omega} \times (\vec{\Omega} \times \vec{r})
$$

---

### Explanation of terms

- $\rho$: fluid density  
- $\vec{u}$: **relative velocity**, i.e., velocity observed in the rotating frame  
- $p$: static pressure  
- $\boldsymbol{\tau}$: viscous stress tensor for a Newtonian fluid:
  $$
  \boldsymbol{\tau} = \mu \left( \nabla \vec{u} + (\nabla \vec{u})^T \right) - \frac{2}{3} \mu (\nabla \cdot \vec{u}) \mathbf{I}
  $$
- $\vec{\Omega}$: angular velocity vector of the rotating frame  
- $\vec{r}$: position vector (from axis of rotation to control volume center)  
- $-2\rho \vec{\Omega} \times \vec{u}$: **Coriolis force**  
- $-\rho \vec{\Omega} \times (\vec{\Omega} \times \vec{r})$: **centrifugal force**

---

In this study, both steady-state and transient simulations are conducted using Ansys CFX pressure-based solver. The flow is described by the three-dimensional, compressible Reynolds-Averaged Navier-Stokes (RANS) equations for the two-phase mixture. The stator domain uses a stationary frame of reference, while the rotor domain is formulated a rotating frame of reference using the relative velocity formulation:
$$
\begin{gather}
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0 \\[1ex]

\frac{\partial (\rho \vec{u})}{\partial t} + \nabla \cdot (\rho \vec{u} \otimes \vec{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \vec{S}_{\Omega}
\end{gather}
$$
Here, $\rho$ and $p$ denote the density and pressure of the two-phase mixture, $\vec{u}$ is the velocity vector relative to the rotating frame, and $\boldsymbol{\tau}$ is the stress tensor, which includes both viscous and turbulent stresses. Turbulence is modeled using the two-equation SST $k$–$\omega$ model with automatic wall functions. The rotational source term of the momentum equation is given by:
$$
\vec{S}_{\Omega} = -2\rho \vec{\Omega} \times \vec{u} - \rho \vec{\Omega} \times (\vec{\Omega} \times \vec{x})
$$
where $\vec{\Omega}$ is the angular velocity of the rotating frame and $\vec{x}$ is the position vector. The first component of $\vec{S}_{\Omega}$ corresponds to the Coriolis force, and the second to the centrifugal force. Both terms vanish in the stationary nozzle domain where $\vec{\Omega} = 0$.

