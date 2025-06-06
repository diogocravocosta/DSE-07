from pyfluids import Fluid, FluidsList, Input

hydrogen = Fluid(FluidsList.Hydrogen).with_state(
    Input.temperature(300), Input.pressure(101325)
)

print(f"Hydrogen Density: {hydrogen.density} kg/m^3")
print(f"Hydrogen Specific Heat: {hydrogen.specific_heat} J/kg-K")
