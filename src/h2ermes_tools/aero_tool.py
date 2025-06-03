import numpy as np
import matplotlib.pyplot as plt

class bluntBody:
    def __init__(self, B, c, vel, rho_at, Rc, a, velocities, rho_atmos, aoa, m):
        self.a = a           # Cone length
        self.B = B           # Cone maximum radius
        self.c = c           # Cone to sphere base (arc height of hemisphere)
        self.vel = vel       # velocity for aerodynamic calculations
        self.rho_at = rho_at # atmospheric density
        self.Rc = Rc         # Cone minimum radius
        self.vels = velocities
        self.rho_atmos = rho_atmos
        self.aoa = aoa #angle of attack during glide
        self.m = m #mass

        #derived geometery
        self.R = (2*B)**2 / (8*c) + c/2          #Radius of the hemisphere
        assert self.R >= B, "Hemisphere radius is too small, must be >= B"
        self.mu_b = np.arccos((self.R - c)/self.R)
        self.alpha_max = np.pi/2 - self.mu_b
        self.w = np.arctan(B / ((B*a) / (B-Rc)))


    def compute_aerodynamics(self, mode):
        if mode == 1: #computing aerodynamics for the structure with respect to all angles of attacks
            self.alpha = np.linspace(-1 * self.alpha_max, self.alpha_max, 100)
            self.alpha = np.linspace(0, self.alpha_max, 100)
            
            self.C_x = 0.5 * (np.sin(self.alpha)**2) * (np.sin(self.mu_b)**2) + \
                (1 + np.cos(self.mu_b)**2) * (np.cos(self.alpha)**2)
            self.C_y = np.sin(self.alpha) * np.cos(self.alpha) * (np.sin(self.mu_b)**2)

            self.C_x_tilde = []
            self.C_y_tilde = []
            for a in self.alpha:
                rho_str = np.arccos(-1 * np.tan(self.w) / np.tan(a))
                print(a, self.w)
                if np.abs(a) >= self.w:
                    C_x_t = -1*((np.sin(a)**2)*(np.cos(self.w)**2)/np.pi)*(((1 + (2*(np.cos(rho_str)**2)))*(np.pi - rho_str)) + (3*np.sin(rho_str)*np.cos(rho_str)))
                    self.C_x_tilde.append(C_x_t)
                    C_y_t = (2 * (np.sin(a))**2 * (np.cos(self.w))**2 / np.pi)*(((np.pi - rho_str)*(np.cos(rho_str))) + ((np.sin(rho_str)/3)*(2 + (np.cos(rho_str))**2)))*(1/np.tan(self.w))
                    self.C_y_tilde.append(C_y_t)
                else:
                    self.C_x_tilde.append(0)
                    self.C_y_tilde.append(0)
            
            self.D_x = 0.5 * self.rho_at * self.vel**2 * np.pi*self.B**2 * (self.C_x + self.C_x_tilde)
            self.L_y = 0.5 * self.rho_at * self.vel**2 * np.pi*self.B**2 * (self.C_y + self.C_y_tilde)

            self.D_x_tot = self.C_x_tilde + self.C_x
            self.C_y_tot = self.C_y + self.C_y_tilde

            #calculating lift curve slope
            self.Cl_alpha = (self.C_y_tot[1] - self.C_y_tot[-1]) / (self.alpha[1] - self.alpha[-1])
            print("Cl_alpha = ", self.Cl_alpha)
        elif mode == 2: #computing aerodynamics for a specific angle of attack for a range of velocities and densities
            self.C_x_aoa = 0.5 * (np.sin(self.aoa)**2) * (np.sin(self.mu_b)**2) + \
                (1 + np.cos(self.mu_b)**2) * (np.cos(self.aoa)**2)
            self.C_y_aoa = np.sin(self.aoa) * np.cos(self.aoa) * (np.sin(self.mu_b)**2)
            rho_str = np.arccos(-1 * np.tan(self.w) / np.tan(self.aoa))
            if self.aoa >= self.w:
                self.C_x_cone_aoa = -1*((np.sin(self.aoa)**2)*(np.cos(self.w)**2)/np.pi)*(((1 + (2*(np.cos(rho_str)**2)))*(np.pi - rho_str)) + (3*np.sin(rho_str)*np.cos(rho_str)))
                self.C_y_cone_aoa = (2 * (np.sin(self.aoa))**2 * (np.cos(self.w))**2 / np.pi)*(((np.pi - rho_str)*(np.cos(rho_str))) + ((np.sin(rho_str)/3)*(2 + (np.cos(rho_str))**2)))*(1/np.tan(self.w))

            self.drags = []
            for i in range(len(self.vels)):
                self.drags.append(0.5 * self.rho_atmos[i] * self.vels[i]**2 * np.pi*self.B**2 * (self.C_x_aoa + self.C_x_cone_aoa))
            

    
    def plots(self, mode):
        if mode == 1:
            plt.plot(self.alpha * 180/np.pi, self.C_x + self.C_x_tilde, label="Cx coefficient", color="blue")
            plt.plot(self.alpha * 180/np.pi, self.C_y + self.C_y_tilde, label="Cy coefficient", color="red")
            #plt.plot(self.alpha * 180/np.pi, self.D_x, label = "Drag force", color = "black")
            #plt.plot(self.alpha * 180/np.pi, self.L_y, label = "Lift force", color = "black")

            plt.legend()
            plt.show()
        elif mode == 2:
            plt.plot(self.drags, self.rho_atmos, label="Drag through flight regime", color="red")
            plt.legend()
            plt.show()

    def draw_geom(self):
        #hemisphere meshgrid
        mu = np.linspace(0, self.mu_b, 100)
        rho = np.linspace(0, 2*np.pi)
        mu,rho = np.meshgrid(mu, rho)

        # Hemisphere parametric equations (front)
        x_sphere = self.R * (1 - np.cos(mu))  # Modified to start at x=R and go backward
        y_sphere = self.R * np.sin(mu) * np.cos(rho)
        z_sphere = self.R * np.sin(mu) * np.sin(rho)

        # Cone meshgrid
        x_cone = np.linspace(c, c+a, 100)
        rho_cone = np.linspace(0, 2*np.pi, 100)
        x_cone, rho_cone = np.meshgrid(x_cone, rho_cone)

        #cone parameteric equations
        r_cone = B - (self.Rc/a)*(x_cone-c)
        y_cone = r_cone * np.cos(rho_cone)
        z_cone = r_cone * np.sin(rho_cone)

        # Plotting
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')

        # Hemisphere (front)
        ax.plot_surface(x_sphere, y_sphere, z_sphere, color='skyblue', alpha=0.9)

        # Cone (extending backward)
        ax.plot_surface(x_cone, y_cone, z_cone, color='lightcoral', alpha=0.9)

        # Labels and aspect
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title('Blunt Body (Hemisphere + Cone)')
        ax.set_box_aspect([2, 1, 1])  # Wider x to show elongation

        plt.show()

    def moments(self):
        Y = np.array(self.C_y)
        Y_tilde = np.array(self.C_y_tilde)
        self.xb = self.R*(1 - np.cos(self.mu_b))
        x_tilde = self.xb + (self.a - (2*self.a / 3*(np.cos(self.w))**2))
        self.Y_cap = 0.5 * self.rho_at * self.vel**2 * np.pi*self.B**2 * Y
        self.Y_cone = 0.5 * self.rho_at * self.vel**2 * np.pi*self.B**2 * Y_tilde
        self.x_moment_arm = ((self.R*self.Y_cap) + (x_tilde*self.Y_cone)) / (self.Y_cap + self.Y_cone)
        
        plt.xlabel("Moment Arm")
        plt.xlabel("Moments Nm")
        plt.ylabel("Angle of Attack")
        #plt.plot(self.x_moment_arm, self.alpha * 180/np.pi, label = "moment arm", color = "black")
        plt.plot((self.R*self.Y_cap), self.alpha * 180/np.pi, label = "moment on cap", color = "red")
        plt.plot((x_tilde*self.Y_cone), self.alpha * 180/np.pi, label = "moment on cone", color = "blue")
        plt.legend()
        plt.show()
    
    def dynamic_stability(self):
        self.I = self.B**2 * self.m * ((3*self.a + 8*self.B)/(10*(self.a + 2*self.B)))
        print(self.I)
        self.rg = np.sqrt(self.I / self.m)
        print(self.rg)
        self.Cmq_Cmadot = ((4*self.I) / (self.m*self.B**2)) * self.Cl_alpha
        print(self.Cmq_Cmadot)
        self.damping = self.D_x_tot - self.Cl_alpha + (2*self.B / self.rg)*self.Cmq_Cmadot
        plt.plot(self.alpha * 180/np.pi, self.damping, label='damping coefficient')
        plt.show()
        #print(self.damping)


if __name__ == "__main__":
    B = 4.94
    c = 2.5
    vel = 4000
    rho_at = 3.31E-04
    Rc = 2.47
    a = 10
    vels = [1000,2000,3000,4000,5000,6000,7000]
    rho_atmos = [0.1, 0.5, 0.7, 0.8, 0.9, 1.1, 1.5]
    aoa = 30 * np.pi/180
    m = 10000

    body = bluntBody(B, c, vel, rho_at, Rc, a, vels, rho_atmos, aoa, m
    )

    body.compute_aerodynamics(mode = 1)

    body.plots(mode = 1)

    body.dynamic_stability()

    #body.draw_geom()