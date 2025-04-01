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
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name, backend="HEOS")
T_in = 37.00 + 273.15
p_in = 91e5
p_out = 1.1*fluid.triple_point_liquid.p

# Create barotropic model object
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    T_in=T_in,
    p_in=p_in,
    p_out=p_out,
    efficiency=0.8,
    calculation_type="equilibrium",
    HEOS_solver="hybr",
    ODE_solver="LSODA",
    ODE_tolerance=1e-9,
    polynomial_degree=5,
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
    plot_spinodal_line=True,
)
for var in model.poly_fitter.variables:
    model.poly_fitter.plot_polynomial_and_error(
        var=var,
        savefig=True,
        showfig=True,
    )

# Show figure
plt.show()
