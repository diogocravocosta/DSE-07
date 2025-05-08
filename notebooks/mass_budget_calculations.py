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

    m_prop = 930 #ton, the mass of propellant, may take into account residual propellant (not being used)
    #m_dry = 1 #ton, the vehicle dry mass
    n_stages = 2 #[-] the number of stages in your vehicle
    f_vac = 1 #N the vaccuum thrust for each engine
    m_fairing = 10.3 #kg/m2, the mass per meter squared of surface area of the faring, S faring in range 15-250 m2
    return (m_prop,)


@app.cell
def _(m_prop):
    #cryogenic stage dry mass N_data = 29

    m_dry = 0.1011*m_prop+1.201 #valid for 8-985 tons, r2 = 0.9914, RSE = 26%





    return


if __name__ == "__main__":
    app.run()
