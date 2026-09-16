from pyfluids import Fluid, FluidsList, Input

coolant = Fluid(FluidsList.Ethanol).with_state(
    Input.pressure(6000000),     
    Input.temperature(25.0),   
)

cp_coolant = coolant.specific_heat
