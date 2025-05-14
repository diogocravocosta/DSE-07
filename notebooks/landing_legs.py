import numpy as np


#constants
g = 9.80665


#model inputs 
h_v = 5 #height of vehicle body including margin for clearance, in m
r_v= 0.2286 #raius of vehicle body, in m
alpha = 0.15 #maximum centre of mass deflection, in deg
m_p = 362.874 #propellant mass left, in kg
m_dry = 181.437 #dry mass of the vehicle, in kg
m = m_p + m_dry #total mass of the vehicle, in kg
h = 0.1524 #drop height, in m
V_o = 0 #initial velocity, in m/s
t = 1 
inner_radius = 0.0254
outer_radius = 0.0508
rho = 2700 #density of aluminum alloy, in kg/m^3
E = 68.9 * 10**9 #modulus of elasticity of aluminum alloy, in N/m^2
K = 1 

I = (np.pi/4) * (outer_radius**4 - inner_radius**4) #moment of inertia for a HOLLOW CYLINDRICAL CROSS-SECTION, in m^4
h_COM = h_v/2 #height of center of mass, in m

h1 = 1.381 #height of primary strut connection, in m
h2 = 1.024 #height of secondary strut connection, in m

V_f = np.sqrt(2*g*h)
F = ((m * (V_o - V_f))/t) #force of impact on legs during landing, in N
P_cr = (F/4)*4 #includes safety margin of 2.11 
f = 2 * h_COM * alpha * 1.67 #footfrint diameter including SF=1.67, in m

s1 = np.sqrt (h1**2+(f/2-r_v)**2) #length of primary strut, in m
s2 = np.sqrt (h2**2+(f/2)**2) #length of secondary strut, in m

I = (K*s1**2)/(E*np.pi**2)

r = ((4*I)/(15*np.pi))**0.25 #inner radius of the strut, in m
R = 2 * r 


mass = 4 *((np.pi*(R**2 - r**2) * s1 * rho)+(np.pi*(R**2 - r**2)*s2*rho)) #mass of the strut, in kg
print (s1, s2, r, R, mass)



