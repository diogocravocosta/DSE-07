import numpy as np
import pandas as pd

#Load excel file with mass inputs
df = pd.read_excel('mass_input_cost_model.xlsx')



#Nomenclature
#DEV = Development cost (engineering cost)
#MAIT = Manufacturing, Assembly, Integration and Test
#ENG = Engineering 
#FM1 = First flight model
#M/PA = Management and Production Assurance cost
#c_p = profit retention cost reduction factor
#STH = total system test hardware
#L_d = development learning factor
#HW = number of times a subsytem is re-used on the launch vehicle

#FIXED VALUES FOR DEV COSTS 
L_d = 0.9 #as a default in Drenthe’s SOLSTICE cost estimating model
s_BAU = 0.6 #scope of subcontracted work under Business As Usual
s_COM = 0.2 #scope of subcontracted work under commercial development
q = 0.08 #average subcontractor profit
c_p = (s_COM*q+1)/(s_BAU*q+1)
M_PA_perc = 5.25 
concept_number = input("Enter concept number (e.g., 1, 2, 3): ")

if concept_number == "1" or concept_number == "2":
    delta_TRL = 4
elif concept_number == "3" or concept_number == "4":
    delta_TRL = 2
elif concept_number == "5" or concept_number == "6":
    delta_TRL = 0
else: 
    print ("Invalid concept number. Please enter a number between 1 and 6.")

#print(f"Development learning factor = {L_d}")
#print(f"Profit retention cost reduction factor = {c_p}")


#LISTS
equipment_names = ["Pressurizant Tank",
    "Fuel Tank",
    "Oxidizer Tank",
    "Thrust Cone",
    "Skirt",
    "Thermal Control",
    "Engine(s)",
    "Thrust Vector Control",
    "Pressurizant System",
    "Pipes",
    "Valves",
    "Stage Harness",
    "Payload Adapter",
    "Payload Fairing",
    "Comms",
    "Power",
    "Data Handling",
    "GNC",
    "Avionics Harness",
    "Attitude Control Module",
    "Interstage Structure"]
a = [19.99465, 19.99465, 19.99465, 2.799300, 2.799300, 2.799300, 31.48271, 33.90978, 11.50618, 8.958770, 8.958770, 27.45211, 26.01794, 23.59239, 51.11253, 42.01174, 141.6820, 69.05491, 27.45211, 257.8420, 6.70369]
b = [0.71253, 0.71253, 0.71253, 0.91199, 0.91199, 0.91199, 0.78811, 0.60977, 1.06948, 0.68815, 0.68815, 0.44623, 0.44623, 0.70000, 0.80000, 0.80000, 0.80000, 0.82458, 0.44623, 0.75000, 0.68041]

if concept_number == "1":
    M = df['CONCEPT 1'].iloc[1:22].tolist()
elif concept_number == "2":
    M = df['CONCEPT 2'].iloc[1:22].tolist()
elif concept_number == "3":
    M = df['CONCEPT 3'].iloc[1:22].tolist()
elif concept_number == "4":
    M = df['CONCEPT 4'].iloc[1:22].tolist()
else: 
    print ("Invalid concept number. Please enter a number between 1 and 4.")
HW = [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25] 
#print (M)
def calc_T1 (a, b, M):
    T1 = []
    for i in range(21):
        C = a[i] * M[i]**b[i]
        T1.append(C)
    T1_total = np.sum(T1)
    return T1

T1 = calc_T1(a, b, M)

def calc_DEV(T1, M_PA_perc, delta_TRL, HW, L_d, c_p):
    DEV_list = []
    for i in range(len(a)):
        print('T1:'+str(T1[i]))
        DD = 3*T1[i]
        print("dd:"+str(DD))
        DM = 0.3 * T1[i]
        print("dm:"+str(DM))
        EM = 1.3 * T1[i]
        print("em:"+str(EM))
        PFM = 1.5 * T1[i]
        print("pfm:"+str(PFM))
        DD_prime = DD + delta_TRL
        print("dd_prime:"+str(DD_prime))
        FM1 = 1 - M_PA_perc/100
        print("fm1:"+str(FM1))
        ENG = DD_prime * FM1
        print("eng:"+str(ENG))
        STH = DM + EM + PFM
        print("sth:"+str(STH))
        MAIT = FM1*STH*L_d*HW[i]
        print("mait:"+str(MAIT))
        DEV = c_p*((ENG+(MAIT+ENG)*M_PA_perc)+(FM1*STH*L_d*HW[i]))
        print("dev:"+str(DEV))
        DEV_list.append(DEV)
        print("dev list:"+str(DEV_list))
        DEV_total = np.sum(DEV_list)
        print("dev total:"+str(DEV_total))
    return DEV_total




def calc_OPS(LpA, f_c, f_8, f_v, Q_N, f_11, L_0, W, N, M_p, M_0, M_press, r, I, P, c_payl, F, T_s, S, c_f, c_ox, c_press):
    #Direct Operations Cost (DOC)
    C_ground = (W * 8 * (M_0**0.67) * (LpA**(-0.9)) * (N ** 0.7) * f_c * f_v * L_0 *f_8 * f_11)/1000 #ground operations in k€
    C_prop = ((M_p/(r +1))*c_f + (M_p - (M_p/(r+1)))* c_ox + M_press * c_press)/1000 #propellant cost in k€
    C_FM = (W * 20 * Q_N * (LpA ** 0.65) * L_0 * f_8)/1000 #flight and mission cost in k€
    C_trans = T_s * M_0 #transportation cost in k€
    C_FI = I + F + (c_payl * P)/1000 #fees & insurance costs in k€
    C_DOC = C_ground + C_prop + C_FM + C_trans + C_FI 
    #Indirect Operations Cost (IOC)
    C_IOC = (40 * S + 24)* (LpA ** (-0.379)) * W / 1000

    C_ops = C_DOC + C_IOC
    return C_ops


#INPUTS FOR OPS COSTS
LpA = 50 #launches per year
W = 301200 #work-year costs in k€
N = 2 #number of stages
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

if concept_number == "1":
    M_p = df.at[25, 'CONCEPT 1'] #mass of propellant and oxidiser in kg
    M_0 = (df.at[24, 'CONCEPT 1'])/1000 #gross take-off-mass Mg
    r = df.at[28, 'CONCEPT 1'] #mass mixture ratio
    c_f = 7.08 #cost of liquefied hydrogen per kg, from https://www.sciencedirect.com/science/article/pii/S2949908923002789
    M_press = df.at[26, 'CONCEPT 1'] #in kg
elif concept_number == "2":
    M_p = df.at[25, 'CONCEPT 2']
    M_0 = (df.at[24, 'CONCEPT 2'])/1000
    r = df.at[28, 'CONCEPT 2']
    c_f = 1.56 #cost of liquefied methane per kg, from https://www.sciencedirect.com/science/article/pii/S1875510021002845
    M_press = df.at[26, 'CONCEPT 2']
elif concept_number == "3":
    M_p = df.at[25, 'CONCEPT 3']
    M_0 = (df.at[24, 'CONCEPT 3'])/1000
    r = df.at[28, 'CONCEPT 3']
    c_f = 7.08
    M_press = df.at[26, 'CONCEPT 3']
elif concept_number == "4":
    M_p = df.at[25, 'CONCEPT 4']
    M_0 = (df.at[24, 'CONCEPT 4'])/1000
    r = df.at[28, 'CONCEPT 4']
    c_f = 1.56
    M_press = df.at[26, 'CONCEPT 4']
else:
    print ("Invalid concept number. Please enter a number between 1 and 6.")

print (T1)
print ("Development cost (engineering cost) in M€= ", (calc_DEV(T1, M_PA_perc, delta_TRL, HW, L_d, c_p))/1000)
print ("Operational cost per flight in M€ = ", (calc_OPS(LpA, f_c, f_8, f_v, Q_N, f_11, L_0, W, N, M_p, M_0, M_press, r, I, P, c_payl, F, T_s, S, c_f, c_ox, c_press))/1000)

