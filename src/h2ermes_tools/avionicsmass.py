from h2ermes_tools.variables import total_dry_mass

def avionicsmass(drymass):
    """
    Sizing from Barry
    All in kg
    """
    return 10 * (drymass ** 0.361)

def wiringmass(drymass, length):
    """
    Sizing from Barry
    Mass in kg, length in m
    """
    return 1.058 * (drymass ** 0.5) * (length ** 0.25)

if __name__ == "__main__":
    print(f"Avionics mass: {avionicsmass(total_dry_mass.value)} kg")
    print(f"Wiring mass: {wiringmass(total_dry_mass.value, 20)} kg")