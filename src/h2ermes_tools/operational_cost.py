import numpy as np
import pandas as pd

#inputs needed from the integration
M_0 = 200 #gross-take off mass in Mg
M_p = 100000 #propellant mass in kg

def calc_OPS(LpA, f_c, f_8, f_v, Q_N, f_11, L_0, W, N, M_p, M_0, M_press, r, I, P, c_payl, F, T_s, S, c_f, c_ox, c_press):
    #Direct Operations Cost (DOC)
    C_ground = (W * 8 * (M_0**0.67) * (LpA**(-0.9)) * (N ** 0.7) * f_c * f_v * L_0 *f_8 * f_11)/1000 #ground operations in k€
    print('C_ground'+str(C_ground/1000))
    C_prop = ((M_p/(r+1))*c_f + (M_p - (M_p/(r+1)))* c_ox + M_press * c_press)/1000 #propellant cost in k€
    print('C_prop'+str(C_prop/1000))
    C_FM = (W * 20 * Q_N * (LpA ** 0.65) * L_0 * f_8)/1000 #flight and mission cost in k€
    print('C_FM'+str(C_FM/1000))
    C_trans = T_s * M_0 #transportation cost in k€
    print('C_trans', C_trans/1000)
    C_FI = I + F + (c_payl * P)/1000 #fees & insurance costs in k€
    print('C_FI', C_FI/1000)
    C_DOC = C_ground + C_prop + C_FM + C_trans + C_FI
    print('C_DOC', C_DOC/1000) 
    #Indirect Operations Cost (IOC)
    C_IOC = (40 * S + 24)* (LpA ** (-0.379)) * W / 1000
    print('C_IOC', C_IOC/1000)
    C_ops = C_DOC + C_IOC + (11367850.00/1000)
    return C_ops

#https://www.ark-invest.com/newsletters/issue-335 refurbishment cost, on the code it has been converted to euros
#INPUTS FOR OPS COSTS
LpA = 25 #launches per year
W = 286.425 #work-year costs in k€
N = 1 #number of stages
f_c = 0.85 #assembly and integration factor
f_v = 1 #launch vehicle type factor, 0.8 for storable propellants, 1 for cryogenic propellants
L_0 = 0.64 #average learning factor operations
f_8 = 1 #country productivity factor 
f_11 = 0.55 #commercial factor
c_ox = 0.12 #according to https://www.sciencedirect.com/science/article/pii/S036031992400627X
P = 12500 #payload mass in kg
c_payl = 5.51 #payload charge site fee in eur per kg
F = 1220 #launch site fee in k€
T_s = 5.365 #specific transportation cost in eur per kg
S = 0.2 #percentage of work subcontracted out 
Q_N = 0.4 * N #vehicle complexity factor
c_press = 35.62 #cost per kg of helium
I = 100 #public damage insurance in M€

#M_0 is the gross-take off mass in Mg (tons)
#M_p is the propellant mass in kg

if __name__ == "__main__":
    ops_cost_per_flight = calc_OPS(LpA, f_c, f_8, f_v, Q_N, f_11, L_0, W, N, M_p, M_0, M_press=0, r=5, I=I, P=P, c_payl=c_payl, F=F, T_s=T_s, S=S, c_f=7.08, c_ox=c_ox, c_press=c_press)
    print(f"Operational cost per flight: {ops_cost_per_flight/1000} M€")