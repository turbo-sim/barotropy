import CoolProp as CP
import barotropy as bpy


# # Superheating examples
# abstract_state = CP.AbstractState("HEOS", "water")
# abstract_state.update(CP.PT_INPUTS, 101325, 120+273.15)
# superheating = bpy.calculate_superheating(abstract_state)
# print(f"Superheating is {superheating:+0.3f} K")

# abstract_state = CP.AbstractState("HEOS", "water")
# abstract_state.update(CP.PQ_INPUTS, 101325, 0.95)
# superheating = bpy.calculate_superheating(abstract_state)
# print(f"Superheating is {superheating:+0.3f} K")


# Subcooling examples
abstract_state = CP.AbstractState("HEOS", "water")
abstract_state.update(CP.PT_INPUTS, 101325, 25+273.15)
subcooling = bpy.calculate_subcooling(abstract_state)
print(f"Subcooling is {subcooling:+0.3f} K")

abstract_state = CP.AbstractState("HEOS", "water")
abstract_state.update(CP.PQ_INPUTS, 101325, 0.05)
subcooling = bpy.calculate_subcooling(abstract_state)
print(f"Subcooling is {subcooling:+0.3f} K")

