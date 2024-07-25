import barotropy as bpy


# State calculation
fluid = bpy.Fluid(name="Water", backend="HEOS")
state = fluid.set_state(bpy.PT_INPUTS, 101325, 300)
print(f"Water density is {state.rho:0.2f} kg/m3 at p={state.p:0.2f} Pa and T={state.T:0.2f} K")

# State calculation
fluid = bpy.Fluid(name="Air", backend="HEOS")
state = fluid.set_state(bpy.PT_INPUTS, 101325, 300)
print(f"Air heat capacity ratio is {state.gamma:0.2f} at p={state.p:0.2f} Pa and T={state.T:0.2f} K")

