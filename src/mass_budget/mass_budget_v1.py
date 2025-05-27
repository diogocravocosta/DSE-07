import numpy as np


#define design variables
propellant = "cryogenic"
rho_propellant = 1.0 #density of propellant in TONS/m^3
n_engines = 2 #number of engines
parachute = True #True or False depending on the re-entry method used
landing_gear = True #True or False depending on the landing method used
fairing = False #True or False depending on the use of a fairing
if fairing == True:
    M_fairing = 30
else:
    M_fairing = 0


#model to be used for mass budgeting (either TU Delft or University of Maryland)
model = "TU Delft"


#CRYOGENIC ROCKET MASS BUDGETING
if propellant == "cryogenic":
    #structural mass in TONS, propellant mass in TONS, M_o in TONS and density in TONS/m^3
    M_structure = 0.1*(M_propellant/(6*rho_propellant))**0.95 + 0.02*M_o

    #tankage mass (aluminum tanks) in kg, V_tank in m^3
    M_tank = (2.903+46.25/np.log(V_tank))*V_tank

    #mass of pump-feed rocket engine in kg, vacuum thrust in Newtons
    M_engine = 8.09*(10**(-3))*F_vac**(0.9)
    M_engines = M_engine*n_engines
    M_engines = M_engines/1000 #convert to TONS

    #avionics mass in TONS, M_structure in TONS, M_engines in TONS
    M_avionics = (5.938-0.00404*(M_structure+M_engines))*(0.01*(M_structure+M_engines))

    #Mass of parachute
    if parachute == true:
        M_parachute = 0.1*M_recovery
    else:
        M_parachute = 0.0
    
    #Mass of landing gear 
    if landing_gear == true:
        M_landing_gear = 0.025*M_landing
    else:
        M_landing_gear = 0.0
    
    #Mass of residual and reserve propellant
    M_propellant_residual = 0.01*M_propellant_usable

else:
    if model == "TU Delft":
        total_mass = M_stage_dry + M_propellant + M_fairing
        if propellant == "liquid oxygen and liquid hydrogen":
            M_stage_dry = 0.1011*M_propellant + 1.201 #in tons, valid in propellant range 8-985 tons

        elif propellant == "semi-cryogenic":
            M_stage_dry = 0.0668*M_propellant + 1.468 #in tons, valid in propellant range 3.5-2040 tons
        
        elif propellant == "liquid hydrazine and liquid nitrogen-tetroxide":
            M_stage_dry = 0.0701 * M_propellant + 0.768 #in tons, valid in propellant range 0.5-420 tons
        
        elif propellant == "solid propellant":
            M_stage_dry = 0.1554 * M_propellant #in kg, valid in propellant range 2-500 tons
        
        M_VEB = 0.345*(M_dry)**0.703
        
        






