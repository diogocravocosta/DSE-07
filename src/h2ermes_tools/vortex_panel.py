#imports
import numpy as np
from aerodynamic_coeffs import BluntBody
import matplotlib.pyplot as plt

class Vortex:
    '''
    Class: Implements the vortex panel method for 2d subsonic flows
    Author: Soham Katewale

    Initialised Variables:
        freestream_v = freestream velocity
        alpha_attack = angle of attack
    '''

    def __init__(self,
                 alpha_attack,
                 freestream_v
                 ):
        
        self.alpha_attack = alpha_attack
        self.freestream_v = freestream_v

    def geometry(self):
        panel_coordinates_upper = [(0,0),
                             (1,1),
                             (2,2),
                             (3,3),
                             (4,4),
                             (5,0),
                             ]
        panel_coordinates_lower = [(0,0),
                                   (1,-1),
                                   (2,-2),
                                   (3,-3),
                                   (4,-4),
                                   (5,0)
                                ]
        panel_coordinates = np.asarray(panel_coordinates_upper[::-1] + panel_coordinates_lower[1:])
        midpoints = []
        for i in range(len(panel_coordinates) - 1):
            x_mid = (panel_coordinates[i][0]+ panel_coordinates[i+1][0]) /2
            y_mid = (panel_coordinates[i][1] + panel_coordinates[i+1][1]) /2
            midpoints.append((x_mid, y_mid))

        #velocity vector
        self.v_inf = np.array([self.freestream_v*np.cos(self.alpha_attack), self.freestream_v*np.sin(self.alpha_attack)])
        
        #generate vectors for each panel
        self.panels = []
        for i in range(len(panel_coordinates) - 1):
            self.panels.append(np.array(panel_coordinates[i+1]) - np.array(panel_coordinates[i]))
        self.panels = np.array(self.panels)

        #generate angles for each panel
        x_axis = [1, 0]
        self.theta = []

        for i in self.panels:
            theta_i = np.arctan2(i[1], i[0])  # angle between vector i and x-axis
            self.theta.append(theta_i)

        self.theta = np.array(self.theta)

        #generate normal vector for each panel
        self.norms = []
        for i in self.panels:
            norm = (-i[1], i[0])
            self.norms.append(norm)
        
        self.norms = np.array(self.norms)
        print(self.panels)

            
        





        #plotting
        nodes_x = [coord[0] for coord in panel_coordinates]
        nodes_y = [coord[1] for coord in panel_coordinates]
        mid_x = [coord[0] for coord in midpoints]
        mid_y = [coord[1] for coord in midpoints]

        plt.plot(nodes_x, nodes_y, label='panels', marker='o', color='red')
        plt.scatter(mid_x, mid_y, label='mids', marker='x', color='blue')
        plt.legend()
        plt.show()



if __name__ == '__main__':

    cone_length = 13.95
    cone_max_radius = 4.95
    cone_min_radius = 2.46
    base_arc_height = 2
    mass = 28000

    aero = Vortex(alpha_attack = 10, freestream_v = 100)
    aero.geometry()
