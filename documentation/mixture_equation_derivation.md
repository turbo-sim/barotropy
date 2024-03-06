

## Abstract
- Using a full multiphase mode not practical if the mixture has a wide distribution
	- Particle size
	- Component densities
- Approximations are possible
	- Multiphase represented by homogeneous single-phase system
	- The influence of the phases is taken into account in the values of physical properties
	- Interaction between phases (such as drag due to difference sin velocity) can be modeled in the flow equations of the equivalent single-phase system
- The multi-phase natura cannot be avoided when the concentration gradients are large and alter the fluid dynamic behavior of the phases
- In many practical/engineering applications the mixture model is a sufficiently accurate approximation
- The report derives the mixture model equations and closure relations
	- Mass and momentum equations for each phase
	- Combination to derive the mixture equation
	- Similar to the single-phase system, but using the density, viscosity and velocity of the mixture
	- Additional terms arise from the slop of the dispersed phases (source terms)
	- In the mixture model the volume fraction of each phase is solved from a phase continuity transport equation
- Closure for the mixture model equations
	- Algebraic relation for the velocity of a dispersed phase with respect to the continuous phase
	- Underlying assumption is the the dispersed phase particles (droplets or bubbles) reach the terminal velocity in a short time period compared to the characteristic time scale of the flow of the mixture
	- **I have to compare the timescales for a supersonic nozzle and the terminal velocity time-scale**
	- The viscous and turbulent stress terms are usually combined and modeled together as a generalized stress
- In multiphase mixtures, acceleration due to gravity/centrifugal/advection tends to cause differences in fluid velocity.
	- There are different formulations and naming conventions
		- Drift flux model
		- Mixture model
		- Algebraic-slip model
		- Others
	- The common factor is that the model equations are
		- One continuity equation for each phase
		- One momentum equation for the mixture
		- The momentum equation contains an additional term representing the effect of velocity differences between phases on the hydrodynamic behavior of the mixture
		- A closure equation modeling the force balance for the dispersed phases is required
		- The constitutive equations for the relative velocities varies in the different mixture models
		- The underlying assumption is local equilibrium is stablished over a short spatial length scale (the particles of the suspended phase reach their terminal velocity rapidly)
- Even if a full multiphase model is theoretically more advanced, the uncertainties in the closure relations can make them less reliable than the simple mixture model
- **I guess the barotropic model goes one step further in reliability, which is essential for complex turbomachinery flows where the main interest is the computation of integral performance parameters**
- Modeling approaches
	- Boltzmann distribution averaging
	- Lagrangian approach
		- Medium is a continuum
		- Dispersed particles are discrete entities
	- Eulerian approach
		- Both medium and phases are a continuun
		- Two-fluid model
		- Mixture model
		- Barotropic model

## Single-phase flow equations
The mass and momentum equations for a single-phase system are given by:
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0
$$
$$
\frac{\partial}{\partial t} (\rho \mathbf{v}) + \nabla \cdot (\rho \mathbf{v} \otimes \mathbf{v}) = -\nabla p + \nabla \cdot (\mathbf{T}_v + \mathbf{T}_t) + \rho \, \mathbf{g}
$$
Where:
- $\rho$ is the density of the fluid.
- $\mathbf{v}$ is the velocity vector of the fluid.
- $p$ is the pressure in the fluid.
- $\mathbf{T}_v$ is the viscous stress tensor.
- $\mathbf{T}_t$ is the turbulence stress tensor in the context of a RANS model.
- $\mathbf{g}$ is a body force field acting on the fluid (e.g., gravity acceleration).

## Multi-phase component equations

In a multi-phase system, the mass and momentum equations for each phase $k$ are given by
$$
\begin{gather}
\frac{\partial}{\partial t}\left(\alpha_k \rho_k\right) + \nabla \cdot (\alpha_k \rho_k \mathbf{v}_k) = \Gamma_k = \sum_{j=1}^n (\Gamma_{k,j}^+ - \Gamma_{j,k}^+)  \\
\frac{\partial}{\partial t}  \left( \alpha_k  \rho_k \mathbf{v}_k  \right) + \nabla \cdot  (\alpha_k \rho_k \mathbf{v}_k \otimes \mathbf{v}_k) = -  \alpha_k \nabla p_k  + \nabla \cdot  [\alpha_k (\mathbf{T}_{v,k} + \mathbf{T}_{t,k})] + \alpha_k \rho_k \, \mathbf{g} + \sum_{j=1}^n (\Gamma_{k,j}^+ \mathbf{u}_{j} - \Gamma_{j,k}^+ \mathbf{u}_{k}) + \mathbf{M}_k
\end{gather}
$$
Where:
- $\alpha_k$ is the volume fraction of phase $k$.
- $\rho_k$ is the density of phase $k$.
- $\mathbf{v}_k$ is the velocity vector of phase $k$.
- $\Gamma_k = \sum_{j=1}^n (\Gamma_{k,j}^+ - \Gamma_{j,k}^+)$ is a represents the rate net mass transfer to phase $k$  from all other phases. The term $\Gamma_{k,j}^+\geq 0$ represents the rate of mass transfer from phase $j$ to phase $k$.
- $\sum_{j=1}^n (\Gamma_{k,j}^+ \mathbf{u}_{j} - \Gamma_{j,k}^+ \mathbf{u}_{i})$ represents the net rate of momentum transfer induced by interphase mass transfer when the phases have different velocities.
- $p_k$ is the pressure of phase $k$.
- $\mathbf{T}_{v,k}$ is the viscous stress tensor of phase $k$.
- $\mathbf{T}_{t,k}$ is the turbulence stress tensor of phase $k$.
- $\mathbf{g}$ is a body force field acting on all phases.
- $\mathbf{M}_k$ is a source source term representing the rate of momentum generation at the interface. This term includes drag forces between the phases and surface tension due to the shape of the interphase. 

The momentum equation can be casted into non-conservative form using the product rule on the left-hand side terms and subtracting the mass transport equation:
$$
 \alpha_k \rho_k \left(\frac{\partial \mathbf{v}_k}{\partial t}  + ( \mathbf{v}_k \cdot \nabla) \mathbf{v}_k \right)= -  \alpha_k \nabla p_k  + \nabla \cdot  [\alpha_k (\mathbf{T}_{v,k} + \mathbf{T}_{t,k})] + \alpha_k \rho_k \, \mathbf{g} + \sum_{j=1}^n \Gamma_{k,j}^+ (\mathbf{u}_{j} - \mathbf{u}_{k}) + \mathbf{M}_k
$$

In order to solve this system of equations we need constitutive relations for:
- $\Gamma_k$  - mass transfer between phases
- $\mathbf{M}_k$ - momentum transfer (i.e., forces) between phases
- $\mathbf{T}_{v,k}$ - viscous stress tensor
- $\mathbf{T}_{t,k}$ - turbulence stress tensor

Even if a full multiphase model is the most sophisticated modeling approach, the uncertainties in the closure relations can make them less reliable than the simple mixture model. In many practical/engineering applications the mixture model is a sufficiently accurate approximation
 

## Multi-phase mixture equations
#### Assumptions:
- Mixture with $n$ phases
- One phase is a continuous fluid
- The other phases are dispersed phases consisting of particles, bubbles, or droplets

#### Mixture mass transport equation
The continuity equation for the mixture model is obtained summing the continuity equation of each phase:
$$
\frac{\partial}{\partial t} \sum_{k=1}^n \left(\alpha_k \rho_k\right) + \nabla \cdot \sum_{k=1}^n (\alpha_k \rho_k \mathbf{v}_k) = \sum_{k=1}^n \Gamma_k
$$
Since the total mass of the system is conserved (e.g., during evaporation the mass of fluid leaving the liquid phase is equal to the mass of fluid entering the vapor phase) the mass source term over all phases is zero:
$$\sum_{k=1}^n \Gamma_k = 0 $$
Defining the mixture density as:
$$ \rho_m = \sum_{k=1}^n \alpha_k \rho_k$$
and the mixture mass flux $\rho_m\mathbf{v}_m$ as:
$$ \mathbf{v}_m = \frac{1}{ \rho_m} \sum_{k=1}^n \alpha_k \rho_k \mathbf{v}_k =   \sum_{k=1}^n y_k \mathbf{v}_k$$
we can write the mass transport equation for the mixture as:
$$
\frac{\partial \rho_m}{\partial t} + \nabla \cdot (\rho_m \mathbf{v}_m) = 0
$$

#### Mixture momentum transport equation
The momentum transport equation for the mixture is also obtained summing the momentum equation of each phase:
$$
\begin{align}
\frac{\partial}{\partial t}  \sum_{k=1}^n\left( \alpha_k  \rho_k \mathbf{v}_k  \right) + \nabla \cdot  \sum_{k=1}^n(\alpha_k \rho_k \mathbf{v}_k \otimes \mathbf{v}_k) =& - \sum_{k=1}^n (\alpha_k \nabla p_k)  + \nabla \cdot  \sum_{k=1}^n[\alpha_k (\mathbf{T}_{v,k} + \mathbf{T}_{t,k})] + \\ 
& +\sum_{k=1}^n (\alpha_k \rho_k \, \mathbf{g}) + \sum_{k=1}^n \sum_{j=1}^n (\Gamma_{k,j}^+ \mathbf{u}_{j} - \Gamma_{j,k}^+ \mathbf{u}_{i}) + \sum_{k=1}^n\mathbf{M}_k
\end{align}
$$

Assuming that the effects of surface tension negligible and that all phases have the same pressure (i.e., mechanical equilibrium) we have that:
$$\nabla p_m = \sum_{k=1}^n \alpha_k \nabla p_k = \nabla p$$
In the absence of surface tension the last term of Eq XX vanishes to zero because any force exerted by one phase is balanced by a reaction force of equal magnitude and opposite sign (Newton's third law):
$$\sum_{k=1}^n\mathbf{M}_k = 0$$
Similarly, the net momentum transfer induced by mass transfer is also zero:
$$ \sum_{k=1}^n \sum_{j=1}^n (\Gamma_{k,j}^+ \mathbf{u}_{j} - \Gamma_{j,k}^+ \mathbf{u}_{i}) = 0$$
In addition, the second term of equation XX can be rewritten in terms of the mixture variables using the concept of diffusion/drift velocity:
$$\mathbf{v}_{d,k} = \mathbf{v}_k - \mathbf{v}_m$$
Introducing the drift velocity definition into the advection term we get:
$$\nabla \cdot  \sum_{k=1}^n(\alpha_k \rho_k \mathbf{v}_k \otimes \mathbf{v}_k) = \nabla \cdot  \sum_{k=1}^n\alpha_k \rho_k (\mathbf{v}_m +  \mathbf{v}_{d,k} ) \otimes (\mathbf{v}_m +  \mathbf{v}_{d,k} )$$
Expanding the outer product and using the fact that $\sum_{k=1}^n\alpha_k \rho_k \mathbf{v}_{d,k}=0$ we can simplify this equation to:
$$\nabla \cdot  \sum_{k=1}^n(\alpha_k \rho_k \mathbf{v}_k \otimes \mathbf{v}_k) = \nabla \cdot (\rho_m \mathbf{v}_m \otimes \mathbf{v}_m) +  \nabla \cdot \sum_{k=1}^n (\alpha_k \rho_k \mathbf{v}_{d,k} \otimes \mathbf{v}_{d,k} )$$
Using these results, we can rewrite the mixture momentum equation as:
$$
\frac{\partial}{\partial t} (\rho_m \mathbf{v}_m) + \nabla \cdot (\rho_m \mathbf{v}_m \otimes \mathbf{v}_m) = -\nabla p + \nabla \cdot (\mathbf{T}_{v,m} + \mathbf{T}_{t,m} + \mathbf{T}_{d,m}) + \rho_m \, \mathbf{g}
$$
where the three stress tensors are given by:
$$\mathbf{T}_{v,m} = \sum_{k=1}^n \alpha_k \mathbf{T}_{v,k}$$
$$\mathbf{T}_{t,m} = \sum_{k=1}^n \alpha_k \mathbf{T}_{t,k}$$
$$\mathbf{T}_{d,m} =  \sum_{k=1}^n \alpha_k \rho_k \mathbf{v}_{d,k} \otimes \mathbf{v}_{d,k} $$
The term $\mathbf{T}_{d,m}$ represents the diffusion of momentum of due to the relative motion (i.e., slip) between phase velocity and mixture velocity.

#### Volume fraction transport equation
We can obtain the continuity equation for a phase $k$ introducing the definition of diffusion/drift velocity to eliminate the phase velocity:
$$
\frac{\partial}{\partial t}\left(\alpha_k \rho_k\right) + \nabla \cdot (\alpha_k \rho_k \mathbf{v}_m) = \Gamma_k - \nabla \cdot (\alpha_k \rho_k \mathbf{v}_{d,k})
$$
Or in terms of mass fractions:
$$
\frac{\partial}{\partial t}\left(y_k \rho_m \right) + \nabla \cdot (y_k \rho_m \mathbf{v}_m) = \Gamma_k - \nabla \cdot (y_k \rho_m \mathbf{v}_{d,k})
$$

#### Summary of the mixture model equations
The mass, momentum and volume fraction transport equations for the mixture model are summarized as:
$$
\begin{gather}
\frac{\partial \rho_m}{\partial t} + \nabla \cdot (\rho_m \mathbf{v}_m) = 0  \\
\frac{\partial}{\partial t} (\rho_m \mathbf{v}_m) + \nabla \cdot (\rho_m \mathbf{v}_m \otimes \mathbf{v}_m) = -\nabla p + \nabla \cdot (\mathbf{T}_{v,m} + \mathbf{T}_{t,m} + \mathbf{T}_{d,m}) + \rho_m \, \mathbf{g}  \\
\frac{\partial}{\partial t}\left(\alpha_k \rho_k\right) + \nabla \cdot (\alpha_k \rho_k \mathbf{v}_m) = \Gamma_k - \nabla \cdot (\alpha_k \rho_k \mathbf{v}_{d,k})  \\
\frac{\partial}{\partial t}\left(y_k \rho_m \right) + \nabla \cdot (y_k \rho_m \mathbf{v}_m) = \Gamma_k - \nabla \cdot (y_k \rho_m \mathbf{v}_{d,k})
\end{gather}
$$
The mass, momentum and volume fraction transport equations for the mixture model can also be written in non-conservative form after some algebraic manipulations:
$$
\begin{gather}
\frac{\partial \rho_m}{\partial t} + \nabla \cdot (\rho_m \mathbf{v}_m) = 0  \\

\rho_m \left(\frac{\partial \mathbf{v}_m}{\partial t}  + ( \mathbf{v}_m \cdot \nabla) \mathbf{v}_m \right) = -\nabla p + \nabla \cdot (\mathbf{T}_{v,m} + \mathbf{T}_{t,m} + \mathbf{T}_{d,m}) + \rho_m \, \mathbf{g}  \\

\rho_m \left(\frac{\partial y_k}{\partial t} + ( \mathbf{v}_m \cdot \nabla) y_k \right) = \Gamma_k - \nabla \cdot (y_k \rho_m \mathbf{v}_{d,k})
\end{gather}
$$
**Observations:**
- The mixture mass transport equation has the same form as the single-phase equation.
- The mixture momentum transport equation has the same form as the single-phase equation with an additional source term due to the difference in velocity between the phases of the system.
- In order to close the system of equations of the mixture we need to provide constitutive relations for:
	- $\Gamma_k$  - mass transfer between phases
	- $\mathbf{v}_{d,k}$ - drift velocity between phases
	- $\mathbf{T}_{v,k}$ - viscous stress tensor
	- $\mathbf{T}_{t,k}$ - turbulence stress tensor

## Constitutive relation for the diffusion velocity
The drift in slip in velocity is usually caused by density differences, which result in forces on the dispersed phases different from those in the continuous phase. This causes differences in velocity that are balanced by drag force between the phases.

The constitutive relations to compute the drift velocity are based on:
- Rigorous analytical manipulation of the flow equations
- Approximation of some terms of the flow equations
- Neglection some terms of the flow equations
- Assuming that the dispersed phase consists of spherical particles of a single, average size and that the drag between particles and continuous phase can be modeled with semi-empirical correlations.
#### Relation between diffusion/drift velocity and slip velocity
The constitutive relations for the relative velocity between the phases are usually expressed in terms of a slip velocity, which is defined as the velocity of one of the dispersed phases relative to the velocity of the continuous phase:
$$\mathbf{v}_{s,k} = \mathbf{v}_{k} - \mathbf{v}_{c}$$
The drift velocity can be computed as a function of the slip velocity according to:
$$
\begin{align*}
\mathbf{v}_{d,k} &= \mathbf{v}_k - \mathbf{v}_m \\
&= (\mathbf{v}_{s,k} + \mathbf{v}_c) - \mathbf{v}_m \\
&= \mathbf{v}_{s,k} + \sum_{i=1}^n y_i( \mathbf{v}_c-\mathbf{v}_i)\\
&= \mathbf{v}_{s,k} - \sum_{i=1}^n y_i\mathbf{v}_{s,i} \quad
\end{align*}
$$
In the above derivation, we've employed the following identities:
$$\sum_{i=1}^n y_i\mathbf{v}_c=\mathbf{v}_c$$
$$ \sum_{i=1}^n y_i \mathbf{v}_i = \mathbf{v}_m$$
#### Drag force experienced by a sphere
The physical reasoning to derive a constitutive equation for the drift velocity is based on a simple model problem, namely the motion of a spherical particle in a fluid at rest under a gravitational field $g$.

The force balance equation for the sphere states that the  product acceleration experienced by the particle and its mass is equal to the sum of forces acting on it:
- Weight force due to the gravity field
- Buoyancy force due to the volume of fluid displaced by the particle
- Drag force due to the friction in the boundary layers surrounding the particle

$$m_p \frac{\mathrm{d}v_{p}}{\mathrm{d}t} = g\,(\rho_p - \rho_c) V_p  - \frac{1}{2}\rho_c A_p C_D \,v_{p}^2$$

where:
- $m_p$ is the mass of the particle
- $v_{p}$ is the velocity of the particle
- $V_p = \frac{4}{3}\pi r_p^3$ is the volume of the particle
- $A_p = \pi r_p^2$ is the cross-section area of the particle
- $\rho_p$ and $\rho_c$ are the densities of the particle and the surrounding fluid, respectively
- $C_D$ is the drag coefficient

The terminal velocity is the velocity achieved by the falling particle when the weight and buoyancy forces are balanced by the drag force such that the acceleration is zero and the particle falls with constant velocity given by:
$$v_{\mathrm{eq}}^2 = \frac{4}{3} \left( \frac{g \, d_p}{C_D} \right) \left(\frac{\rho_p - \rho_c}{\rho_c} \right)$$
In the limit of low-speed, laminar flow the drag coefficient has an analytical solution given by:
$$C_D = \frac{24}{\mathrm{Re}_p}$$
where:
$$\mathrm{Re}_p = \frac{\rho_c \, v_p \, d_p}{\mu_c}$$
Under these conditions, the equation describing the motion of the fluid is a linear, first-order ODE with an analytical solution given by an exponential. The time constant of the exponential solution is given by:
$$ \tau_p = \frac{\rho_p d_p^2}{18 \mu_c}$$
The drag coefficient given by eq ZZ is only valid under the assumption of Stokes flow (very low Reynolds numbers). The equation is not accurate when the Reynolds number increases and a wake is formed downstream the sphere.

Schiller and Nauman proposed a semi-empirical equation to correlate the dependency of the drag coefficient with the Reynolds number over a wide range of Reynolds numbers:
$$ C_D = \begin{cases} 
\frac{24}{\mathrm{Re}_p}\cdot(1 + 0.15 \mathrm{Re}_p^{0.687}) & \text{if } \mathrm{Re}_p < 1000 \\ 
0.44 & \text{if } \mathrm{Re}_p \geq 1000 \end{cases} $$
When the drag coefficient is given by this nonlinear relation the equation does not longer have an analytic solution and has to be solved by means of numerical integration.

This equation is valid for a single particle falling in a still fluid.
In the case of the multi-phase mixture model the velocity of the particle is replaced by the relative velocity between the phase $k$ and the continuous phase $k=c$ (i.e., the slip velocity)

In a multi-phase system the motion of the particle is influenced by the flow field distortion caused by the other particles. There are several empirical correlations to account for this effect that are based on modifications of the drag coefficient formula or the definition of the viscosity employed in the Reynolds number.

#### Rate of momentum generated at the interface
The momentum equation for phase $k=p$ has a source term $\mathbf{M}_k$ representing the rate of momentum generation at the interface. This term accounts for the drag force caused by the relative velocity between the dispersed phase and the continuous phase and it is sometimes referred to as the drag-induced momentum transfer.

In the previous section we considered a model problem and discussed how to calculate the drag force and terminal velocity for a sphere moving in a fluid at rest under the acceleration of gravity. We can use the same physical principles to estimate the drag-induced momentum transfer. This source represents the sum of drag forces experienced by each particle within in a differential volume of fluid:
$$\mathbf{M}_k = \frac{\sum \mathbf{F}_{D,p}}{\delta V} = \frac{N_p \cdot \mathbf{F}_{D,p}}{\delta V}$$
where $N_p$ is the number of particles contained in the volume  $\delta V$. The drag force experienced by each particle has a magnitude proportional to the square of relative velocity between particle and fluid and direction opposite to the relative velocity:
$$\mathbf{F}_{D,p} = -\frac{1}{2} \rho_c A_p C_D \, |\mathbf{v}_p - \mathbf{v}_c|(\mathbf{v}_p - \mathbf{v}_c)$$
The volume fraction of phase $p$ is defined as the ratio of the volume occupied by the particles to total volume:
$$\alpha_p = \frac{\delta V_p}{\delta V} = \frac{N_p \cdot V_p}{\delta V}$$
Using this identity we can rewrite the source term using familiar flow variables:
$$\mathbf{M}_k = \frac{\alpha_p}{V_p} \cdot F_{D,p} = \frac{\alpha_p}{V_p} \cdot \mathbf{F}_{D,p}$$
$$\mathbf{M}_k = -\frac{3}{4} \left(\frac{\alpha_p \, \rho_c}{d_p} \right) C_D \, |\mathbf{v}_p - \mathbf{v}_c|(\mathbf{v}_p - \mathbf{v}_c)$$
Or in terms of the slip velocity:
$$\mathbf{M}_k = -\frac{3}{4} \left(\frac{\alpha_p \, \rho_c}{d_p} \right) C_D \, |\mathbf{v}_{s,k} | \, \mathbf{v}_{s,k}$$

The same relation can also be derived by considering the drag force per unit of volume by dividing the drag force felt by a particle by the volume of the particle and scaling according to the volume fraction:
$$\mathbf{M}_k  =\frac{\mathbf{F}_{D,p}}{V_p} \frac{\delta V_p}{\delta V} = - \left( \frac{ \alpha_p \rho_c  A_p}{2 V_P}\right)\, C_D\, |\mathbf{v}_p - \mathbf{v}_c|(\mathbf{v}_p - \mathbf{v}_c)$$
See Manninen et al. (1996) and Ishii & Mishima (1984) for more details.

#### Constitutive relation for the drift velocity
The closure equation to calculate the relative (diffusion/drift) velocity can be rigorously derived combining the momentum equations for the dispersed phase and the mixture. The momentum equation in non-conservative form for phase $k = p$ is given by
$$
\alpha_p  \rho_p \left(\frac{\partial \mathbf{v}_p}{\partial t} + (\mathbf{v}_p \cdot \nabla ) \,\mathbf{v}_p \right)= - \alpha_p\nabla p + \nabla \cdot [\alpha_p (\mathbf{T}_{v,p} + \mathbf{T}_{t,p})] + \alpha_p \rho_p \, \mathbf{g} +\mathbf{M}_p +\sum_{j=1}^n \Gamma_{k,j}^+ (\mathbf{u}_{j} - \mathbf{u}_{k})  
$$
The mixture momentum equation multiplied by $\alpha_p$ is given by:
$$
\alpha_p \rho_m \left(\frac{\partial  \mathbf{v}_m}{\partial t} + ( \mathbf{v}_m \cdot \nabla )  \mathbf{v}_m \right) = - \alpha_p\nabla p + \alpha_p \nabla \cdot (\mathbf{T}_{v,m} + \mathbf{T}_{t,m} + \mathbf{T}_{d,m}) + \alpha_p\rho_m \, \mathbf{g}
$$
Subtracting the mixture equation from the phase equation and solving for the drag-induced momentum source term we obtain:
$$
\begin{align}
\mathbf{M}_p =\; & + \alpha_p \left[  \rho_p\frac{\partial \mathbf{v}_p}{\partial t} -  \rho_m\frac{\partial \mathbf{v}_m}{\partial t}  \right]  \\
& +\alpha_p \left[  \rho_p (\mathbf{v}_p \cdot \nabla ) \,\mathbf{v}_p) - \rho_m (\mathbf{v}_m \cdot \nabla ) \,\mathbf{v}_m)\right] \\
& - \nabla \cdot [\alpha_p (\mathbf{T}_{v,p} + \mathbf{T}_{t,p})] + \alpha_p \nabla \cdot (\mathbf{T}_{v,m} + \mathbf{T}_{t,m} + \mathbf{T}_{d,m}) \\
& - \alpha_p (\rho_p -\rho_m) \, \mathbf{g} \\
& +\sum_{j=1}^n \Gamma_{k,j}^+ (\mathbf{u}_{j} - \mathbf{u}_{k})

\end{align}
$$
Next, we make several approximations to simplify this equation:
- We neglect the viscous, diffusion, and turbulence stresses
- We neglect the momentum transfer induced by interphase mass transfer
- We assume that the advective velocity of the particle is equal to the advective velocity of the mixture
$$ (\mathbf{v}_p \cdot \nabla ) \,\mathbf{v}_p \approx (\mathbf{v}_m \cdot \nabla ) \,\mathbf{v}_m$$
Under these assumptions we have that:
$$
\begin{align}
\mathbf{M}_p =\; & + \alpha_p \left[  \rho_p\frac{\partial \mathbf{v}_p}{\partial t} -  \rho_m\frac{\partial \mathbf{v}_m}{\partial t}  \right]  \\
& +\alpha_p (\rho_p  - \rho_m) (\mathbf{v}_m \cdot \nabla ) \,\mathbf{v}_m) \\
& - \alpha_p (\rho_p -\rho_m) \, \mathbf{g}
\end{align}
$$

Introducing the drift velocity to replace the time derivative of the particle velocity we get

$$
\begin{align}
\mathbf{M}_p =\; & + \alpha_p \left[  \rho_p\frac{\partial \mathbf{v}_{d,p}}{\partial t} +  (\rho_p - \rho_m)\frac{\partial \mathbf{v}_m}{\partial t}  \right]  \\
& +\alpha_p (\rho_p  - \rho_m) (\mathbf{v}_m \cdot \nabla ) \,\mathbf{v}_m) \\
& - \alpha_p (\rho_p -\rho_m) \, \mathbf{g} \\\
\end{align}
$$
If furthermore we assume that the particle phase is in local equilibrium (i.e., the particles reach the terminal velocity) we have that $\frac{\partial \mathbf{v}_{d,p}}{\partial t}=0$, leading to:
$$

\mathbf{M}_p = \alpha_p \, (\rho_p - \rho_m) \left[\frac{\partial \mathbf{v}_m}{\partial t}  +  (\mathbf{v}_m \cdot \nabla ) \,\mathbf{v}_m) - \mathbf{g} \right]

$$
This can be interpreted as the rate of momentum transfer to phase $k$ being proportional to:
- the void fraction of phase $p$
- the difference of density between phase $p$ and the mixture
- the local acceleration of the mixture

If we compare the rate of momentum transfer obtained from this derivation and the one obtained from the drag force equation we obtain the constitutive relation for the diffusion/drift velocity:
$$ |\mathbf{v}_{s,k} | \, \mathbf{v}_{s,k} = -\frac{4}{3} \left(\frac{d_p}{C_D}\right) \left(\frac{\rho_p - \rho_m}{\rho_c}\right) \left[\frac{\partial \mathbf{v}_m}{\partial t}  +  (\mathbf{v}_m \cdot \nabla ) \,\mathbf{v}_m) - \mathbf{g} \right]$$
This is the final equation that lets us compute the drift velocity

Note that this is an implicit equation because the drag coefficient $C_D$ is itself a function of the drift velocity through its dependence on the Reynolds number.

## Constitutive relation for viscous shear stress

Not trivial how to define these stresses
It is common to assume one-phase stress tensor with Stokes hypothesis using the molecular viscosity obtained as the volumetric average of the phase viscosities.

> Alternatively, the closure law for the viscous stress can be formulated by considering the mixture as a single-phase fluid and determining the viscous stress tensor analogously to single-fluid flow in terms of the mixture parameters.

## Constitutive relation for turbulence shear stress

> The turbulent (Reynolds) stresses are caused by the fluctuations of the velocity relative to the mean velocity. The turbulent term is in general more important in multiphase flow than in single-phase flow. Even if the freestream turbulence is negligible, the flow around individual particles can generate velocity fluctuations

Boussinesq’s assumption for the Reynolds stresses is usually extended to multiphase systems.

Application of the generalised stress is most straightforward, especially in turbulent flows. Uncertainties in the influence of turbulence are naturally significant because of our limited capability of modelling turbulent flows in general.