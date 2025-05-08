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
    return


@app.cell
def _():
    # Zandbergen Methods

    #assuming new glenn with barge landing

    #initial inputs of values

    options = ['no', 'Steel', 'Composite'] #choose no for no particular material, Steel or composite if you know a material

    engine = True # if true, use more granular method

    S_fairing = 1 #m2 S faring in range 15-250 m2
    m_prop = 930*1000 #kg, the mass of propellant, may take into account residual propellant (not being used)
    n_stages = 2 #[-] the number of stages in your vehicle
    f_vac = 1 #N the vaccuum thrust for each engine

    choice = options[0]

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
    

    VEB = 0.345*(m_dry)**0.703 #kg | n_data = 11, r2 = 0.9381, RSE = 72% | vehicle equipment bay





    return


@app.cell
def _():





    return


if __name__ == "__main__":
    app.run()
