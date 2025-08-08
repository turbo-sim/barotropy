
import os
import barotropy as bpy
import matplotlib.pyplot as plt
 
bpy.print_package_info()
bpy.set_plot_options(grid=False)
 
# Create folder to save results
DIR_OUT = "./barotropic_model"
save_figures = True
if not os.path.isdir(DIR_OUT):
    os.makedirs(DIR_OUT)
 
# Define inlet state and outlet pressure
fluid_name = "MDM"
fluid = bpy.Fluid(name=fluid_name, backend="REFPROP")
T_in = 269+273.15
p_in = 9.02e5
state_in= bpy.compute_properties_coolprop(fluid._AS, bpy.PT_INPUTS, p_in, T_in)

# p_out = fluid.triple_point_liquid.p*25
# p_out = 1e+5/60
p_out = 0.1e4
 
# Create barotropic model object
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    backend="REFPROP",
    p_in=p_in,
    rho_in=state_in["rho"],
    p_out=p_out,
    efficiency=1,
    calculation_type="equilibrium",
    HEOS_solver="hybr",
    ODE_solver="LSODA",
    ODE_tolerance=1e-12,
    polynomial_degree=12,
    polynomial_format="horner",
    output_dir=DIR_OUT,
    polynomial_variables=["density", "viscosity", "speed_sound", "void_fraction", "vapor_quality"],
 
   
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
    dT_crit=2,
)
 
for var in model.poly_fitter.variables:
    model.poly_fitter.plot_polynomial_and_error(
        var=var,
        savefig=True,
        showfig=True,
    )
 
# Show figure
plt.show()