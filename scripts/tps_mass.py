def estimate_thickness(heat_load, surface_area, entry_speed) -> float:
    """
    Estimate the thickness of a TPS (Thermal Protection System) using the heat load and entry speed.

    Parameters:
    heat_load (float): The heat load in kW/m^2.
    entry_speed (float): The entry speed in m/s.

    Returns:
    float: The estimated thickness of the PICA TPS in m.
    """
    # Calculate the mass using Mooij's formula
    thickness = 1.8696*(heat_load * surface_area / entry_speed**2)**0.1879
    
    # Convert thickness to meters
    thickness = thickness * 1e-2  # Convert from cm to m

    return thickness

if __name__ == "__main__":
    # Example values
    heat_load = 32e6  # J/m^2
    entry_speed = 7200  # m/s

    PICA_DENSITY = 270  # kg/m^3

    # calculate surface area, assume a circle
    diameter = 7.0  # m
    surface_area = 3.14 * (diameter**2) / 4  # m^2

    # convert heat load to J/cm^2 and entry speed to km/s as needed by the model
    heat_load = heat_load / 1e4  # J/cm^2
    entry_speed = entry_speed / 1000  # km/s

    # calculate thickness
    thickness = estimate_thickness(heat_load, surface_area, entry_speed)
    surface_area = 2 * 3.14 * (1.5**2)  # m^2

    # Ensure the thickness is not less than 3.27 cm
    if thickness < 3.27e-2:
        thickness = 3.27e-2

    print(thickness)
    print("Estimated mass of PICA TPS: ", thickness * surface_area * PICA_DENSITY, "kg")
