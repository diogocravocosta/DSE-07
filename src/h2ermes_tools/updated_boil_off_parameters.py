import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import pytest
from data import material as mat
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
#------------------------------------------------
#Functions
#------------------------------------------------

def mass_hydrogen_vapor_after_insertion(tank_volume,m_h2_propellant_left,rho_h2,rho_gaseous_hydrogen_launchpad):
    volume_vapor_before_insertion = 0.1*tank_volume
    volume_vapor_after_insertion = tank_volume - m_h2_propellant_left/rho_h2
    mass_hydrogen_vapor_before_insertion = 215
    mass_hydrogen_vapor_after_2 = volume_vapor_after_insertion*mass_hydrogen_vapor_before_insertion/volume_vapor_before_insertion
    mass_hydrogen_vapor_after_insertion = rho_gaseous_hydrogen_launchpad*volume_vapor_after_insertion

    return mass_hydrogen_vapor_after_insertion

mass_hydrogen_vapor_launch = mass_hydrogen_vapor_after_insertion(tank_volume=887,m_h2_propellant_left=24674,rho_h2=67.7,rho_gaseous_hydrogen_launchpad=2.49)

def pressure_change_during_refueling(m_gh2_before_refuelling,a_coefficient,b_coefficient,tank_volume,volume_vapor_before_refuelling,pressure_before_refuelling,m_h2_after_refuelling,rho_lh2_before_refuelling):
    #assume no boil off
    volume_vapor_after_refuelling = tank_volume - m_h2_after_refuelling / rho_lh2_before_refuelling
    n_gh2_before_refuelling = m_gh2_before_refuelling / 0.002
    
    pressure_term_before_refuelling = pressure_before_refuelling + a_coefficient * n_gh2_before_refuelling **2 / volume_vapor_before_refuelling **2
    volume_term_before_refuelling = volume_vapor_before_refuelling - n_gh2_before_refuelling * b_coefficient
    volume_term_after_refuelling = volume_vapor_after_refuelling - n_gh2_before_refuelling * b_coefficient
    
    pressure_term_after_refuelling = pressure_term_before_refuelling * volume_term_before_refuelling / volume_term_after_refuelling
    pressure_after_refuelling = pressure_term_after_refuelling - a_coefficient * n_gh2_before_refuelling **2 / volume_vapor_after_refuelling **2

    return pressure_after_refuelling

pressure_after_refuelling = pressure_change_during_refueling(m_gh2_before_refuelling = 377, 
                                                             a_coefficient = 0.0246,
                                                             b_coefficient = 2.661e-5,
                                                             tank_volume = 887,
                                                             volume_vapor_before_refuelling = 434.88, 
                                                             pressure_before_refuelling = 697000, 
                                                             m_h2_after_refuelling = 5798, 
                                                             rho_lh2_before_refuelling = 63.3)


print("pressure inside tank when no boil off occurs",pressure_after_refuelling)

#4437 mass of vapor right after refuelling 

# 1300 kg mass of vapor after launch

def mass_vapor_after_refuelling(tank_volume,rho_gh2_after_refuelling,m_hydrogen_after_refuelling,rho_lh2_after_refuelling):
    #assume no pressure change

    volume_vapor_after_refuelling = tank_volume - m_hydrogen_after_refuelling/rho_lh2_after_refuelling

    mass_hydrogen_vapor_after_refuelling = volume_vapor_after_refuelling*rho_gh2_after_refuelling#volume_vapor_after_refuelling*mass_hydrogen_vapor_before_refuelling/volume_vapor_before_refuelling
    return mass_hydrogen_vapor_after_refuelling

payload_mass_delievered = [16753.25,16979.6,17206.4,17480,17685.2]
time_to_refuel = [6,8,10,12,14]

def quadractic_regression(payload_mass_delievered,time_to_refuel):
    time_to_refuel = np.array(time_to_refuel).reshape(-1, 1)
    payload_mass = np.array(payload_mass_delievered)  # Non-linear target

# Generate quadratic features (degree 2)
    poly = PolynomialFeatures(degree=2)
    X_poly = poly.fit_transform(time_to_refuel)  # adds x^0, x^1, x^2

    # Fit the model
    model = LinearRegression()
    model.fit(X_poly, payload_mass)

    # Predict and compute MSE
    payload_mass_pred = model.predict(X_poly)
    mse = mean_squared_error(payload_mass, payload_mass_pred)

    # Print coefficients and error
    print("Coefficients:", model.coef_)        # [bias, x^1 coeff, x^2 coeff]
    print("Intercept:", model.intercept_)      # same as model.coef_[0] in some contexts
    print("Mean Squared Error:", mse)
    # Plot data and regression curve
    plt.figure(figsize=(8, 5))
    plt.scatter(time_to_refuel, payload_mass, color='blue', label='Data points')

    # Create a smooth curve for plotting
    x_range = np.linspace(min(time_to_refuel), max(time_to_refuel), 300).reshape(-1, 1)
    x_range_poly = poly.transform(x_range)
    y_range_pred = model.predict(x_range_poly)

    plt.plot(x_range, y_range_pred, color='red', label='Quadratic Fit')
    plt.xlabel("Time to Refuel")
    plt.ylabel("Payload Mass Delivered")
    plt.title("Quadratic Regression: Payload vs Refuel Time")
    plt.legend()
    plt.grid(True)
    plt.show()
    return model

def linear_regression(payload_mass_delievered,time_to_refuel):
    time_to_refuel = np.array(time_to_refuel).reshape(-1, 1)
    payload_mass = np.array(payload_mass_delievered)

    model = LinearRegression()
    model.fit(time_to_refuel, payload_mass)

    payload_mass_pred = model.predict(time_to_refuel)
    mse = mean_squared_error(payload_mass, payload_mass_pred)

    print("Coefficient:", model.coef_)
    print("Intercept:", model.intercept_)
    print("Mean Squared Error:", mse)
    print("R2 Score:", model.score(time_to_refuel, payload_mass))

    plt.figure(figsize=(8, 5))
    plt.scatter(time_to_refuel, payload_mass, color='blue', label='Data points')
    plt.plot(time_to_refuel, payload_mass_pred, color='green', label='Linear Fit')
    plt.xlabel("Time to refuel [hours]")
    plt.ylabel("Payload mass [kg] delivered to depot")
    plt.title("Linear Regression: Payload vs Refuel Time")
    plt.legend()
    plt.grid(True)
    plt.show()
    return model


payload_mass_delievered = [16753.25,16979.6,17206.4,17480,17685.2]
time_to_refuel = [6,8,10,12,14]

payload_mass = [item /0.95 for item in payload_mass_delievered]

boil_off = [24674-item-3000 +215 for item in payload_mass]
print(boil_off)
p_insertion = 2

V_gh2_insertion = (56700) 
#p_start = p_insertion*V_gh2_insertion/V_gh2_start


#linear_model = linear_regression(payload_mass_delievered,time_to_refuel)