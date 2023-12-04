


## Polynomial representation of the barotropic model

In the context of CFD simulations, the direct evaluation of fluid properties throught the equations of state implemented in software libraries like REFPROP or CoolProp often leads to prohibitive computational costs. To address this limitation, various techniques have been proposed to approximate the fluid properties. Among these, interpolating values from pre-computed look-up tables is a commonly used method [citation?]. Another approach is the adoption of surrogate models tailored to specific regions of the equation of state [citation?]. For the barotropic model, where all fluid properties depend only on pressure, single-variable polynomials emerge as a suitable surrogate to approximate the fluid properties. The benefits of this approach include:

- **Accuracy:** A high-order polynomial fit can closely match the equation of state with minimal deviations.
- **Smoothness**: Polynomials can produce a smooth curve over the entire domain. In contrast, some interpolation methods might produce curves that are not smooth at the data points.
- **Computational Efficiency**: The computational cost of evaluation a polynomial is usually negligible compared to the evaluation of multi-parameter equations of state.
- **Ease of Implementation**: Polynomials expressions can be easily integrated into CFD codes, for example, through user-defined functions.


Based on these advantages, we opted to use a polynomial representation of the fluid properties. 
The procedure to generate the polynomial expressions for the barotropic model is described below.


### Fluid Property Calculation
The first step is integrating Equation XX as described in section X. Through this integration, we obtain the sequence of fluid property values across the pressure interval from $p_1$ to $p_N$:

$$
\phi_i = \{\phi_1,\, \phi_2,\, \phi_3, \dots,\, \phi_N\}
$$

Here, $\phi$ represents any fluid property, such as density or viscosity. During this step, fluid properties are evaluated directly from the equation of state, for instance, using REFPROP.

### Polynomial Fitting

Once the fluid properties over the specified pressure interval are computed, we proceed to approximate them with polynomials.

For any property, $\phi$, a polynomial segment in _power form_ is given by:

$$
P_{\phi}(x) = \sum_{k=0}^{n} a_k x^k =  a_n x^n + a_{n-1} x^{n-1} + \ldots + a_2 x^2 + a_1 x + a_0
$$

where $a_k$ are the fitting coefficients, $x = \frac{p}{p_{\mathrm{ref}}}$ is a normalized pressure, and $p_{\mathrm{ref}}$ is a reference pressure (like the inlet pressure or the critical pressure). Utilizing a normalized pressure, rather than the absolute pressure, helps to improve conditioning of the polynomial fit. The coefficients are determined by minimizing the squared differences between the calculated fluid properties and the polynomial expression:

$$
E =  \sum_{i=1}^{N} \left( P_{\phi}(x_i) - \phi_i \right)^2
$$

This least-squares fitting is achieved using the `numpy.polyfit` function.


### Polynomial evaluation
In practice, the polynomial expressions are not implemented in _power form_, but in _Horner's form_:

$$ P(x) = a_n x^n + x \left( a_{n-1} x^{n-1} + x \left( \ldots + x \left( a_2 x^2 + x (a_1 x + a_0) \right) \ldots \right) \right) $$

Evaluating a polynomial in Horner's form has several computational advantages:
- **Efficient Evaluation**: The primary advantage of Horner's form is that it allows for a more efficient evaluation of a polynomial at a given point. In the standard polynomial representation, evaluation requires n multiplications and n-1 additions. In Horner's form, it requires n multiplications and n additions, which can be faster in many cases, especially when n is large.
- **Numerical Stability**: When evaluating polynomials using floating-point arithmetic, Horner's method tends to be more numerically stable. This is especially important for polynomials of higher degrees, where the standard method might introduce significant errors.
- **Numerical Accuracy**: For certain polynomials, especially when coefficients have very different magnitudes, using Horner's form can reduce rounding and truncation errors, leading to more accurate results.


### Polynomial extrapolation

CFD solvers may occasionally yield negative pressures in some cells of the flow domain. Although these negative pressures are unphysical and should not appear in the final solution, they can manifest during the convergence process. To ensure a stable solution of the flow equations, the barotropic model must provide meaningful property values, even when evaluated at negative pressures. Consequently, the polynomial expressions were extended to pressures below $p_1$ using an exponential extrapolation:

$$
\phi(x) = 
\begin{cases} 
\alpha e^{\beta\,(x-x_1)} & \text{if } x < x_1 \\
P_{\phi}(x) & \text{if } x \geq x_1
\end{cases}
$$

In the equation above $x_1$ represents the reduced pressure corresponding to $p_1$ and the parameters $\alpha$ and $\beta$ are determined such that the function $\phi$ has first-order smoothness at the breakpoint $x_1$:

$$
\begin{gather*}
\alpha = P_{\phi}(x_1) \\
\beta = \frac{\dot{P_{\phi}}(x_1)}{P_{\phi}(x_1)}
\end{gather*}
$$

The rationale behind choosing exponential extrapolation is clear when considering the density of the fluid. This method ensures the density remains positive across all pressure values, approaching zero as the pressure tends towards negative infinity. Furthermore, it guarantees a monotonous decrease in density with a decline in pressure, leading to a well-defined speed of sound.



### Phase change and metastable states

In the case of single-component systems with phase change (e.g., flashing of liquid into vapor) fluid properties are fitted using a piecewise polynomial with two segments:

$$
P_{\phi}(x) = 
\begin{cases} 
P_{\phi}^{2p}(x) & \text{if } x < x^* \\
P_{\phi}^{1p}(x) & \text{if } x \geq x^*
\end{cases}
$$

where:
 - $P_{\phi}^{1p}(x)$ is the polynomial for single-phase region (including metastable states).
 - $P_{\phi}^{2p}(x)$ is the polynomial for the two-phase region.
 - $x^*$ is the reduced pressure at which the phase change is triggered (somewhere between the saturation pressure and the spinodal pressure)

Using two-polynomial segments in this case is necessary to achieve good fitting accuracy because the variation of fluid properties is significantly different inside and outside the two-phase region.

To do:
- Add links to metastable property calculations
- Explain the blending of the metastable spates