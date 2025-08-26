# %%
import os
import matplotlib.pyplot as plt


import barotropy as bpy

bpy.print_package_info()

bpy.set_plot_options()


# ------------------------------ #
# ------ Barotropic Model ------ #
# ------------------------------ #

dir_barotropic = "barotropic_kaist"
if not os.path.isdir(dir_barotropic):
    os.makedirs(dir_barotropic)

# Create fluid object
fluid_name = "CO2"
fluid = bpy.Fluid(name=fluid_name)
polytropic_efficiency = 1.0
calculation_type = "blending"
T_in = 35.06 +273.15
p_in = 76.00e5
p_out = 89.76e5
q_onset = 0.98
q_transition = 0.03

# Create barotropic model object
model = bpy.BarotropicModel(
    fluid_name=fluid_name,
    T_in=T_in,
    p_in=p_in,
    p_min=1.2*fluid.triple_point_liquid.p,
    p_max=1.5*p_in,
    process_type="compression",
    efficiency=polytropic_efficiency,
    calculation_type=calculation_type,
    blending_onset=q_onset,
    blending_width=q_transition,
    HEOS_solver="hybr",
    ODE_solver="LSODA",
    ODE_tolerance=1.0e-9,
    polynomial_degree=[4, 4, 4],
    polynomial_format="horner",
    polynomial_variables=["density", "viscosity", "speed_sound", "void_fraction", "vapor_quality"]
)

# Evaluate barotropic model and export polynomial expressions
model.solve()
model.fit_polynomials()
model.export_expressions_fluent(output_dir=dir_barotropic)
model.export_expressions_cfx(output_dir=dir_barotropic)
model.poly_fitter.plot_phase_diagram(
    fluid=fluid,
    var_x="s",
    var_y="T",
    savefig=True,
    showfig=True,
)
for var in model.poly_fitter.variables:
    model.poly_fitter.plot_polynomial_and_error(
        var=var, savefig=True, showfig=True
    )


# Show figures
plt.show()
