import numpy as np

#material density [kg/m3]
#material strength in cryogenic conditions [MPa]
materials_properties = {
    "Al-Li": {"density": 3000, "strength": 560},
    "Annealed 304L Stainless Steel": {"density": 7800, "strength": 190},
    "Annealed Ti-6Al-4V": {"density": 4400, "strength": 910},
    "Inconel 718": {"density": 8300, "strength": 1190}
}

propellant_mass = 0
propellant_density = 0
propellant_pressure = 0
tank_diameter = 7

def tank_volume():
    ##
# calculating tank thickness
def calculate_tank_thickness(tank_diameter, tank_length, propellant_pressure, allowable_stress):
    # Assuming a thin-walled cylinder for simplicity
    # Using the formula: t = (P * D) / (2 * σ)
    # where P is the internal pressure, D is the diameter, and σ is the allowable stress
    thickness = (propellant_pressure * tank_diameter) / (2 * allowable_stress)
    return thickness

def calculate_tank_mass(tank_diameter, tank_length, thickness, material_density):
    # Calculate the volume of the tank material (assuming a cylindrical shape)
    outer_radius = tank_diameter / 2 + thickness
    inner_radius = tank_diameter / 2
    volume = np.pi * (outer_radius**2 - inner_radius**2) * tank_length

    mass = volume * material_density
    return mass

