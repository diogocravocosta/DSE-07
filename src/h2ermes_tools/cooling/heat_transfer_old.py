import numpy as np
import thermo
import scipy.optimize 

# TODO
# output setperate heat fluxes for inner and outer channel
# write a coolant wrapper class that includes Re, v etc



class HeatTransfer():
    def __init__(self, cea, gas, geometry, material, coolant_inner, coolant_outer, cooling_geometry, m_dot, m_dot_coolant, output, settings2D, model='standard-bartz', cool_model='gnielinski', eta_c_star=0.92):
        """[summary]
        Class to calculate heat tranfer coefficients and radiative heat transfer, pass those onto the 2D simulation and 
        advance to the next chamber section. Several options are available for coolant and heat transfer models. 
        Coolant flow properties stored as total conditions, and use a 'thermo' Mixture object to update coolant properties.        
        Hot gas properties from CEA and taken at total chamber conditions. Temperature and pressure of the gas are determined 
        using isentropic expansion.        
        Material thermal condutivity can be temperature dependent.
        USE SI UNITS (sorry freedom lovers)
        """
        self.cea = cea                       			# NASA CEA object 
        self.gas = gas                      			# isentropic gas object 
        self.geometry = geometry             			# engine contour
        self.material = material             			# material object
        self.coolant_inner = coolant_inner              # themro Mixture object
        self.coolant_outer = coolant_outer            	# themro Mixture object
        self.m_dot = m_dot                   			# total mass flow
        self.m_dot_coolant = m_dot_coolant   			# coolant mass flow [kg/s]
        self.model = model                   			# hot gas side heat transfer model
        self.cool_model = cool_model		 			# coolant side heat transfer model
        self.eta_c_star = eta_c_star       			    # c* efficiency
        self.cooling_geometry = cooling_geometry	 	# cooling channel geometry class

        self.A_c_inner = 2 * cooling_geometry.Ai_arr * cooling_geometry.N       # area of inner cooling channels [m^2]
        self.A_c_outer = 2 * cooling_geometry.Ao_arr * cooling_geometry.N       # area of outer cooling channels [m^2]
        self.D_h_inner = cooling_geometry.dhi_arr			                    # cooling channel hydraulic diameter [m]
        self.D_h_outer = cooling_geometry.dho_arr			                    # cooling channel hydraulic diameter [m]
        self.m_dot_c_inner = self.m_dot_coolant * self.A_c_inner[-1] / (self.A_c_inner[-1] + self.A_c_outer[-1])
        self.m_dot_c_outer = self.m_dot_coolant * self.A_c_outer[-1] / (self.A_c_inner[-1] + self.A_c_outer[-1])

        self.T_amb = 288.0								# ambient temperature [K], must be float!!!
        self.boltzmann = 5.67e-8						# stefan boltzmann constant

        self.out = output								# output class
        self.settings2D = settings2D					# settings for 2D thermal sim


    def heat_trans_coeff_gas(self, T_wall, idx):
        # get hot gas properties for this chamber location. CEA properties assumed constant, with isentropic gas properties being taken at local point
        gamma = self.cea.gamma
        mu    = self.cea.mu
        cp    = self.cea.Cp
        Pr    = self.cea.Pr
        T_s   = self.gas.T_s[idx]
        p_s   = self.gas.p_s[idx]
        M     = self.gas.M[idx]
        T_aw  = self.gas.T_aw[idx]

        # local geometrical properties 
        local_area  = np.pi * self.geometry[idx,1]**2
        throat_area = np.pi * np.amin(self.geometry[:,1])**2
        D_t         = 2 * np.amin(self.geometry[:,1])
        
        cstar = self.cea.chamber_pressure * throat_area / self.m_dot
        
        # hot gas temperature estimate used for the cinjarev correlation
        T_hg       = T_s + 0.9*(T_aw * self.eta_c_star**2 - T_s)

        # Nusselt number correlation for the combustion gases
        if self.model == 'standard-bartz':
            # standard bartz correlation, mostly applicable for larger engines. Tends to oveerpredict heat flux near the troat, especially in smaller engines
            sigma  = ((0.5 * T_wall/self.cea.Tc * (1+(gamma-1)/2 * M**2) + 0.5)**(0.68) * (1+(gamma-1)/2 * M**2)**(0.12))**(-1)
            halpha = 0.026/(D_t**0.2) * (mu**0.2)*cp/(Pr**0.6) * (self.gas.p_t/cstar)**0.8 * (throat_area/local_area)**0.9 * sigma

        elif self.model == 'modified-bartz':
            # modified bartz correlation, mostly applicable for larger engines. Tends to oveerpredict heat flux near the troat, especially in smaller engines
            T_f    = 0.5 * T_wall + 0.28 * T_s + 0.22 * T_aw
            G      = self.m_dot/local_area
            halpha = 0.026 * G**0.8/(D_t)**0.2 * mu**0.2*cp/Pr**0.6 * (self.cea.Tc/T_f)**0.68

        elif self.model == 'cinjarev':
            # cinjarev correlation, more useful for medium to small scale motors
            halpha = 0.01975 * self.cea.k**0.18 * (self.m_dot*cp)**0.82 / (2*self.geometry[idx,1])**1.82 * (T_hg/T_wall)**0.35

        else:
            raise ValueError('Invalid heat transfer method. Select: "standard-bartz", "modified-bartz" or "cinjarev"')

        return halpha, T_hg
        
        
    def pressure_drop(self, idx, coolant, D_h, Re, v_coolant):
        # using Darcy Friction Factor to determine pressure drop. Use surface roughness given in Material Class
        surface_roughness = self.material.roughness
        fd = scipy.optimize.fsolve(lambda f: -2 * np.log10(surface_roughness / (3.7 * D_h[idx]) + 2.51 / (Re * np.sqrt(f))) -1 / np.sqrt(f), 0.001)
        dp = fd * self.section_length[idx] / D_h[idx] * 0.5 * coolant.rho * v_coolant**2 
        
        return dp, fd

        
    def radiation(self, idx):
        # simple radiation model based on h2o and co2 emissions
        R_c = self.geometry[idx,1]					# chamber radius

        p_h2o = self.cea.mole_fractions[1]['H2O'][0] * self.gas.p_s[idx]
        # radiative heat flux of water molecules 
        q_rad = 5.74 * (p_h2o/1e5*R_c)**0.3 * (self.gas.T_s[idx]/100)**3.5

        # check if CO2 is present in the exhaust
        if '*CO2' in self.cea.mole_fractions[1]:

            p_co2 = self.cea.mole_fractions[1]['*CO2'][0] * self.gas.p_s[idx]
            # radiative heat flux of co2 molecules 
            q_rad += 4 * (p_co2/1e5*R_c)**0.3 * (self.gas.T_s[idx]/100)**3.5

        return q_rad


    def heat_trans_coeff_coolant_inner(self, T_wall_coolant, idx):
        # coolant flow velocity based on bulk properteis
        self.v_coolant_inner = self.m_dot_c_inner / (self.coolant_inner.rho * self.A_c_inner[idx])
        self.Re_inner        = self.coolant_inner.rho * self.v_coolant_inner * self.D_h_inner[idx] / self.coolant_inner.mu
			
		# thermodynamic properties of the near wall fluid, Pr implementation of thermo does not work reliably for mixtures near their critical point. Uses gaseouse Pr if liquid Pr returns 'None'
        Pr = self.coolant_inner.Pr
        Cp = self.coolant_inner.Cp
        mu = self.coolant_inner.mu

        # coolant thermal conductivity
        k  = Cp * mu / Pr

        # Nusselt number correlations for the coolant side 
        if self.cool_model == 'dittus-boelter':
            # Dittus-Boelter equation (valid for Re > 1e4 and 0.7 < Pr < 16700)
            T_avg = (T_wall_coolant - self.coolant_inner.T) / np.log(T_wall_coolant / self.coolant_inner.T) 
            near_wall_coolant = thermo.Mixture(IDs=self.coolant_inner.IDs, ws=self.coolant_inner.ws, T=T_avg, P=self.coolant_inner.P)
            Nu = 0.027 * self.Re_inner**0.8 * Pr**0.33 * (near_wall_coolant.mu / mu)**(0.14)		

        elif self.cool_model == 'dittus-boelter-simple':
            # Dittus-Boelter equation without viscosity correction (valid for Re > 1e4 and 0.7 < Pr < 160)
            Nu = 0.023 * self.Re_inner**0.8 * Pr**0.4

        elif self.cool_model == 'gnielinski':
            # Gnielinski correlation (valid for 3000 < Re < 5e6 and 0.5 < Pr < 2000)
            _, fd = self.pressure_drop(idx, self.coolant_inner, self.D_h_inner, self.Re_inner, self.v_coolant_inner)
            Nu = (fd/8) * (self.Re_inner - 1000)*Pr / (1 + 12.7*(fd/8)**0.5 * (Pr**(2/3) - 1))

        else:
            raise ValueError('Invalid heat transfer method. Select: "dittus-boelter", "dittus-boelter-simple" or "gnielinski"')
        
        # convert Nusselt number to heat transfer coefficient
        halpha = Nu * k / self.D_h_inner[idx]

        if halpha < 0:
            raise ValueError('Negative heat transfer coefficient, check applicability of cooling model to Reynolds number range')

        # check if halpha is an array or a float (artifact of the pressure drop dependence or fluid model)
        if isinstance(halpha, float):
            return halpha
        else:
            return halpha[0]

    
    def heat_trans_coeff_coolant_outer(self, T_wall_coolant, idx):
        # coolant flow velocity based on bulk properteis
        self.v_coolant_outer = self.m_dot_c_outer / (self.coolant_outer.rho * self.A_c_outer[idx])
        self.Re_outer        = self.coolant_outer.rho * self.v_coolant_outer * self.D_h_outer[idx] / self.coolant_outer.mu
			
		# thermodynamic properties of the near wall fluid, Pr implementation of thermo does not work reliably for mixtures near their critical point. Uses gaseouse Pr if liquid Pr returns 'None'
        Pr = self.coolant_outer.Pr
        Cp = self.coolant_outer.Cp
        mu = self.coolant_outer.mu

        # coolant thermal conductivity
        k  = Cp * mu / Pr

        # Nusselt number correlations for the coolant side 
        if self.cool_model == 'dittus-boelter':
            # Dittus-Boelter equation (valid for Re > 1e4 and 0.7 < Pr < 16700)
            T_avg = (T_wall_coolant - self.coolant_outer.T) / np.log(T_wall_coolant / self.coolant_outer.T) 
            near_wall_coolant = thermo.Mixture(IDs=self.coolant_outer.IDs, ws=self.coolant_outer.ws, T=T_avg, P=self.coolant_outer.P)
            Nu = 0.027 * self.Re_outer**0.8 * Pr**0.33 * (near_wall_coolant.mu / mu)**(0.14)		

        elif self.cool_model == 'dittus-boelter-simple':
            # Dittus-Boelter equation without viscosity correction (valid for Re > 1e4 and 0.7 < Pr < 160)
            Nu = 0.023 * self.Re_outer**0.8 * Pr**0.4

        elif self.cool_model == 'gnielinski':
            # Gnielinski correlation (valid for 3000 < Re < 5e6 and 0.5 < Pr < 2000)
            _, fd = self.pressure_drop(idx, self.coolant_outer, self.D_h_outer, self.Re_outer, self.v_coolant_outer)
            Nu = (fd/8) * (self.Re_outer - 1000)*Pr / (1 + 12.7*(fd/8)**0.5 * (Pr**(2/3) - 1))

        else:
            raise ValueError('Invalid heat transfer method. Select: "dittus-boelter", "dittus-boelter-simple" or "gnielinski"')
        
        # convert Nusselt number to heat transfer coefficient
        halpha = Nu * k / self.D_h_outer[idx]

        if halpha < 0:
            raise ValueError('Negative heat transfer coefficient, check applicability of cooling model to Reynolds number range')

        # check if halpha is an array or a float (artifact of the pressure drop dependence or fluid model)
        if isinstance(halpha, float):
            return halpha
        else:
            return halpha[0]


    def section_heat_flux(self, idx):
        self.q_rad = self.radiation(idx)
        
        # solve 2D section temperature profile, passes functions for heat transfer coefficients for temperature dependence of boundary conditions
                
        print('SOLVING Section Number	', len(self.geometry[:,1]) - idx, ' / ', len(self.geometry[:,1]))
        
        solver = HeatEquationSolver(idx, self.gas, self.material, self.cooling_geometry, self.heat_trans_coeff_gas, self.heat_trans_coeff_coolant_inner, self.heat_trans_coeff_coolant_outer, self.q_rad, self.coolant_inner.T, self.coolant_outer.T, self.T_amb, self.out.folder_path, self.settings2D)
        solver.run_sim()

        self.T_wall_i  = solver.inner_wall_temp	    # inner wall temperature taken as maximum temperature of the 2D profile
        self.halpha    = solver.halpha			    # heat transfer coefficients taken as the of the boundary conditions of the 2D solver at the last timestep
        self.halpha_c_bottom  = solver.halpha_c_bottom		# average heat transfer coefficient at bottom wall
        self.halpha_c_top  = solver.halpha_c_top
        self.halpha_c_side1  = solver.halpha_c_side1
        self.halpha_c_side2  = solver.halpha_c_side2
        self.T_hg      = solver.T_hg

        # total heat flux into the inner chamber wall
        self.q 		   = np.mean(solver.q_chamber)

        # radiation leaving to ambient 
        self.q_rad_out = np.mean(solver.q_outer)
	
        # ensure that FiPy output object is actually a float, multiplied by the number of channels for total area
        Q_inner = float(solver.Q_c_inner) * self.section_length[idx] * self.cooling_geometry.N 
        Q_outer = float(solver.Q_c_outer) * self.section_length[idx] * self.cooling_geometry.N 
 
        # update cooling fluid and pressure. Expected temperature rise and pressure drop 
        dT_inner             = Q_inner / (self.m_dot_c_inner * self.coolant_inner.Cp) 
        dp_inner, _          = self.pressure_drop(idx, self.coolant_inner, self.D_h_inner, self.Re_inner, self.v_coolant_inner)
        dT_outer             = Q_outer / (self.m_dot_c_outer * self.coolant_outer.Cp) 
        dp_outer, _          = self.pressure_drop(idx, self.coolant_outer, self.D_h_outer, self.Re_outer, self.v_coolant_outer)
        
        # recalcualte cooling fluid properties, such as density, Pr, Cp etc. 
        self.coolant_inner   = thermo.Mixture(IDs=self.coolant_inner.IDs, ws=self.coolant_inner.ws, T=(self.coolant_inner.T+dT_inner), P=(self.coolant_inner.P-dp_inner))
        self.coolant_outer   = thermo.Mixture(IDs=self.coolant_outer.IDs, ws=self.coolant_outer.ws, T=(self.coolant_outer.T+dT_outer), P=(self.coolant_outer.P-dp_outer))
    
        
    def run(self):
        # arrays for section length and area in between the 2D sections solved in the thermal sim 
        self.section_length     = np.ndarray(len(self.geometry[:,1]))

        # determine the section length and inner chamber surface area in each section
        for i in range(len(self.geometry[:,1])):
        # section number 0 at the injector!!!
			
            if self.settings2D.start_idx == -1:
				# for coolant entering at the bottom of the cooling channels, use idx as inverse of i 
                idx = len(self.geometry[:,1]) - i - 1
			
            elif self.settings2D.start_idx == 0:
				# in case coolant starts at injector side
                idx = i
			
            else:
                raise ValueError('Invalid starting point for cooling fluid, select -1 or 0 for nozzle or injector side, respectively')

            if i == 0:
                self.section_length[idx] = 0
                
            else:
                x_len = (self.geometry[idx+1,0] - self.geometry[idx,0])**2
                y_len = (self.geometry[idx+1,1] - self.geometry[idx,1])**2
                
                self.section_length[idx] = np.sqrt(x_len + y_len)

            # solve section heat transfer and temperature field and update coolant properties
            self.section_heat_flux(idx)
            
            # write to output class
            self.out.halpha[idx]          = self.halpha
            self.out.q_rad[idx]           = self.q_rad
            self.out.q[idx]               = self.q
            self.out.T_wall_i[idx]        = self.T_wall_i
            self.out.halpha_c_bottom[idx] = self.halpha_c_bottom
            self.out.halpha_c_top[idx]    = self.halpha_c_top
            self.out.halpha_c_side1[idx]  = self.halpha_c_side1
            self.out.halpha_c_side2[idx]  = self.halpha_c_side2
            self.out.T_c_inner[idx]       = self.coolant_inner.T
            self.out.P_c_inner[idx]       = self.coolant_inner.P
            self.out.T_c_outer[idx]       = self.coolant_outer.T
            self.out.P_c_outer[idx]       = self.coolant_outer.P
            self.out.Re_inner[idx]	      = self.Re_inner
            self.out.Re_outer[idx]	      = self.Re_outer
            self.out.T_hg[idx]	          = self.T_hg 
            self.out.v_coolant_inner[idx] = self.v_coolant_inner
            self.out.v_coolant_outer[idx] = self.v_coolant_outer
            self.out.T_c_avg[idx]         = (self.m_dot_c_outer * self.coolant_outer.T + self.m_dot_c_inner * self.coolant_inner.T) / (self.m_dot_c_outer + self.m_dot_c_inner)

        # write to file in output folder
        self.out.write_csv()

