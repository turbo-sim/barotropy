import os
import barotropy as bpy
import matplotlib.pyplot as plt
 
import jaxprop as jxp

# Initialize script
bpy.print_package_info()
bpy.set_plot_options(grid=False)
save_figures = True
DIR_OUT = "./output"
os.makedirs(DIR_OUT, exist_ok=True)

# Define inlet state and outlet pressure
fluid_name = "cyclopentane"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")
p_in = 0.15 * fluid.critical_point.p
state_in = fluid.get_state(jxp.PQ_INPUTS, p_in, 0.1)
# state_in = fluid.get_state(jxp.PT_INPUTS, p_in, state_in.T-5)
# p_out = 1.2*fluid.critical_point.p
p_out = state_in.p / 20
 
# Create barotropic model object
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    p_in=state_in["p"],
    rho_in=state_in["rho"],
    p_out=p_out,
    efficiency=1.00,
    calculation_type="equilibrium",
    HEOS_solver="hybr",
    ODE_solver="LSODA",
    ODE_tolerance=1e-9,
    polynomial_degree=[8],
    polynomial_format="horner",
    output_dir=DIR_OUT,
    # polynomial_variables=["density"],
    polynomial_variables=["density", "viscosity", "speed_of_sound", "void_fraction", "vapor_quality"],
    
)
 
# Evaluate barotropic model and export polynomial expressions
model.solve()
model.fit_polynomials()
model.export_expressions_fluent(output_dir=DIR_OUT)
model.export_expressions_cfx(output_dir=DIR_OUT)
model.poly_fitter.plot_phase_diagram(
    fluid=fluid,
    var_x="s",
    var_y="T",
    savefig=True,
    showfig=True,
    plot_spinodal_line=False,
    dT_crit=3.0
)
for var in model.poly_fitter.variables:
    model.poly_fitter.plot_polynomial_and_error(
        var=var,
        savefig=True,
        showfig=True,
    )


# Show figures
if not os.environ.get("DISABLE_PLOTS"):
    plt.show()

 