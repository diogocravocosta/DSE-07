# /// script
# [tool.marimo.runtime]
# auto_instantiate = false
# ///

import marimo

__generated_with = "0.13.6"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    return (np,)


@app.function
# Zandbergen Methods

def zandbergen(option,S, mp, n_s,eng_bool, f_v, n_e): #input which option you want, faring surface, mass of propellant, number of stages, take into account engines True/False, optional if engine true: vaccuum thrust, number of engines
    #initial inputs of values
    
    options = ['no', 'Steel', 'Composite'] #choose no for no particular material, Steel or composite if you know a material
    choice = options[option]
    
    engine = eng_bool # if true, use more granular method
    
    S_fairing = S #m2 S faring in range 15-250 m2
    m_prop = mp #kg, the mass of propellant, may take into account residual propellant (not being used)
    n_stages = n_s #[-] the number of stages in your vehicle
    f_vac = f_v #N the vaccuum thrust for each engine
    n_engine = n_e #[-] number of engines, if you want to calculate using the granular method
    
    
    
    #--------------------------------------------------------------------------------------
    # Results Calculation Zone
    #--------------------------------------------------------------------------------------
    
    m_fairing = 10.3*S_fairing #kg/m2, the mass per meter squared of surface area of the faring
    
    if not engine:
        if choice == 'no':
            m_dry = (0.1011*m_prop/1000+1.201) *1000 #kg | valid for 8-985 tons, r2 = 0.9914, RSE = 26% | cryogenic stage dry mass N_data = 29
        elif choice == 'Steel':
            m_dry = 0.1726*m_prop*1000.+1173 #kg | n_data = 11, r2 = 0.9907, RSE = 20.1% | 
        elif choice == 'Composite':
            m_dry == 0.1173*m_prop #kg | n_data = 11, r2 = 0.9833, RSE = 24.1% | 
    elif engine:
        c_mass = (0.0872*m_prop/1000 + 1.013) * 1000 #kg | n_data = 11, r2 = 0.9833, RSE = 24.1% | construction mass, dry mass minus total engine mass
        m_eng = 0.0016*f_vac +36.1
        m_dry = c_mass+m_eng
    
    VEB = 0.345*(m_dry)**0.703 #kg | n_data = 11, r2 = 0.9381, RSE = 72% | vehicle equipment bay


@app.function
# Maryland 'murica Methods (yee haw)

#mass of lox, mass of lh2, surface area of lox and lh2 tanks, mass of gas (if gas tank used), surface area of fairing, engine thrust, area expansion ratio, estimation of total rocket mass


def maryland(m_lox, m_lh2, S_lox, S_lh2, m_gas, S_fairing, f_eng, exp_ratio, m_tot):
    m_lox_tank = 0.0152*m_lox+318 #kg
    m_lox_insulation = 1.123*S_lox #kg
    m_lh2_tank = 0.0694*m_lh2+363 #kg
    m_lh2_insulation = 2.88*S_lh2 #kg
    if m_gas:
        m_gas_tank = 2 #kg
    m_eng = 7.81*10**-4 * f_eng + 3.37 *10**-5 * f_eng * exp_ratio**0.5 + 59 #kg
    m_thrust_structure = 2.55*10**-4 * f_eng #kg
    m_fairing = 13.3 * S_fairing
    m_avionics = 10*(m_tot)**0.361
    m_wiring = 1.058 * m_tot **(1/8)


@app.cell
def _(m_recovery, np):
    # Apel Methods


    def model_d(delta_v, isp, m_prop, rho_prop, v_tank, f_vac, m_0, m_recovered, m_landing, m_prop_usable): # delta v, isp, propellant mass, propellant mass density, tank volume, vaccuum thrust, initial vehicle mass, vehicle mass recovered, vehicle mass at landing, usable propellant mass
        m_struct = 1000*(0.1*((m_prop/1000)/(6*rho_prop))**0.95 + 0.02*m_0) #kg
        m_tank = (2.903+46.25/(np.log(v_tank)))*v_tank #kg
        m_eng = 8.09*10**-3 * f_vac**0.9 #kg zandbergen says this relation is sus
        m_avionics = 1000*(5.938-0.00404*(m_struct+m_eng))*(0.01*(m_struct+m_eng))
        m_parachute = 0.1*m_recovery
        m_landing_gear = 0.03* m_landing
        m_prop_residual = 0.01*m_prop_usable
    return


if __name__ == "__main__":
    app.run()
