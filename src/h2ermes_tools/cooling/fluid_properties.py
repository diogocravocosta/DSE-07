import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

# --- Single Temperature Calculation ---

temperature = 20  # K

vapor_pressure = CP.PropsSI('P', 'T', temperature, 'Q', 0, 'H2')

print(f"Temperature: {temperature} K")
print(f"Vapor Pressure: {vapor_pressure:,.2f} Pa")
print(f"Vapor Pressure: {vapor_pressure / 100000:,.4f} bar")
print(f"Vapor Pressure: {vapor_pressure / 101325:,.4f} atm")

print("\n" + "="*30 + "\n")

temperatures = np.linspace(13.9, 33.145, 100)  # K
pressures = []

for T in temperatures:
    P = CP.PropsSI('P', 'T', T, 'Q', 0, 'H2')
    pressures.append(P)

plt.figure(figsize=(10, 6))
plt.plot(temperatures, np.array(pressures) / 100000) # Convert Pa to bar for the plot
plt.title('Vapor Pressure of Liquid Hydrogen')
plt.xlabel('Temperature (K)')
plt.ylabel('Vapor Pressure (bar)')
plt.grid(True)
plt.show()
