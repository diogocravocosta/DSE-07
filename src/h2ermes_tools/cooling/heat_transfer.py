import numpy as np
import scipy.optimize

from data.constants import boltzmann
from h2ermes_tools.cooling.coolant import Coolant
from pyfluids import Fluid, FluidsList, Input


class HeatShield:
    def __init__(
        self,
        incident_heat_flux,
        contour,
        material,
        coolant,
        cooling_model="taylor",
    ):
        """
        Class to calculate heat tranfer coefficients and radiative heat transfer, pass those onto the 2D simulation and
        advance to the next chamber section. Several options are available for coolant and heat transfer models.
        Coolant flow properties stored as total conditions, and use a 'thermo' Mixture object to update coolant properties.
        Hot gas properties from CEA and taken at total chamber conditions. Temperature and pressure of the gas are determined
        using isentropic expansion.
        Material thermal condutivity can be temperature dependent.
        USE SI UNITS (sorry freedom lovers)

        This class combines heat shield cooling calculations. It calculates coolant heat transfer coefficients, pressure drops and
        radiative heat transfer.
        """
        self.ihf = incident_heat_flux  # incident heat flux from the combustion gases
        self.contour = contour  # heat shield contour
        self.material = material  # material object
        self.coolant = coolant  # Coolant object

    def pressure_drop(self, idx, coolant, hydraulic_diameter, reynolds_number, coolant_speed):
        # using Darcy Friction Factor to determine pressure drop. Use surface roughness given in Material Class
        surface_roughness = self.material.Ra  # effective roughness height in m
        f_D = scipy.optimize.fsolve(
            lambda f: -2
            * np.log10(surface_roughness / (3.7 * hydraulic_diameter[idx]) + 2.51 / (reynolds_number * np.sqrt(f)))
            - 1 / np.sqrt(f),
            0.001,
        )
        dp = f_D * self.section_length[idx] / hydraulic_diameter[idx] * 0.5 * coolant.rho * coolant_speed**2

        return dp, f_D  # pressure drop and Darcy friction factor

    def heat_transfer_coefficient_coolant(self, idx, coolant_wall_temperature):
        # coolant flow velocity based on bulk properteis
        self.v_coolant_inner = self.m_dot_c_inner / (
            self.coolant_inner.rho * self.A_c_inner[idx]
        )
        self.Re_inner = (
            self.coolant_inner.rho
            * self.v_coolant_inner
            * self.D_h_inner[idx]
            / self.coolant_inner.mu
        )

        # thermodynamic properties of the near wall fluid, Pr implementation of thermo does not work reliably for mixtures near their critical point. Uses gaseouse Pr if liquid Pr returns 'None'
        Pr = self.coolant.prandtl
        Cp = self.coolant.specific_heat
        mu = self.coolant.dynamic_viscosity

        # coolant thermal conductivity
        k = Cp * mu / Pr

        # Nusselt number correlations for the coolant side
        if self.cooling_model == "dittus-boelter":
            # Dittus-Boelter equation (valid for Re > 1e4 and 0.7 < Pr < 16700)
            T_avg = (coolant_wall_temperature - self.coolant_inner.T) / np.log(
                coolant_wall_temperature / self.coolant_inner.T
            )
            near_wall_coolant = self.coolant.with_state(
                Input.temperature(T_avg), Input.pressure(self.coolant_inner.P)
            )
            near_wall_coolant = thermo.Mixture(
                IDs=self.coolant_inner.IDs,
                ws=self.coolant_inner.ws,
                T=T_avg,
                P=self.coolant_inner.P,
            )
            Nu = (
                0.027
                * self.Re_inner**0.8
                * Pr**0.33
                * (near_wall_coolant.mu / mu) ** (0.14)
            )

        elif self.cooling_model == "dittus-boelter-simple":
            # Dittus-Boelter equation without viscosity correction (valid for Re > 1e4 and 0.7 < Pr < 160)
            Nu = 0.023 * self.Re_inner**0.8 * Pr**0.4

        elif self.cooling_model == "gnielinski":
            # Gnielinski correlation (valid for 3000 < Re < 5e6 and 0.5 < Pr < 2000)
            _, fd = self.pressure_drop(
                idx,
                self.coolant_inner,
                self.D_h_inner,
                self.Re_inner,
                self.v_coolant_inner,
            )
            Nu = (
                (fd / 8)
                * (self.Re_inner - 1000)
                * Pr
                / (1 + 12.7 * (fd / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
            )

        elif self.cooling_model == "taylor":
            Nu = self.coolant.get_nusselt_number_taylor(
                self.Re_inner
            )

        else:
            raise ValueError(
                'Invalid heat transfer method. Select: "dittus-boelter", "dittus-boelter-simple" or "gnielinski"'
            )

        # convert Nusselt number to heat transfer coefficient
        h_c = Nu * k / self.D_h_inner[idx]

        if h_c < 0:
            raise ValueError(
                "Negative heat transfer coefficient, check applicability of cooling model to Reynolds number range"
            )

        # check if h_c is an array or a float (artifact of the pressure drop dependence or fluid model)
        if isinstance(h_c, float):
            return h_c
        else:
            return h_c[0]

    def heat_trans_coeff_coolant_outer(self, T_wall_coolant, idx):
        # coolant flow velocity based on bulk properteis
        self.v_coolant_outer = self.m_dot_c_outer / (
            self.coolant_outer.rho * self.A_c_outer[idx]
        )
        self.Re_outer = (
            self.coolant_outer.rho
            * self.v_coolant_outer
            * self.D_h_outer[idx]
            / self.coolant_outer.mu
        )

        # thermodynamic properties of the near wall fluid, Pr implementation of thermo does not work reliably for mixtures near their critical point. Uses gaseouse Pr if liquid Pr returns 'None'
        Pr = self.coolant_outer.Pr
        Cp = self.coolant_outer.Cp
        mu = self.coolant_outer.mu

        # coolant thermal conductivity
        k = Cp * mu / Pr

        # Nusselt number correlations for the coolant side
        if self.cool_model == "dittus-boelter":
            # Dittus-Boelter equation (valid for Re > 1e4 and 0.7 < Pr < 16700)
            T_avg = (T_wall_coolant - self.coolant_outer.T) / np.log(
                T_wall_coolant / self.coolant_outer.T
            )
            near_wall_coolant = thermo.Mixture(
                IDs=self.coolant_outer.IDs,
                ws=self.coolant_outer.ws,
                T=T_avg,
                P=self.coolant_outer.P,
            )
            Nu = (
                0.027
                * self.Re_outer**0.8
                * Pr**0.33
                * (near_wall_coolant.mu / mu) ** (0.14)
            )

        elif self.cool_model == "dittus-boelter-simple":
            # Dittus-Boelter equation without viscosity correction (valid for Re > 1e4 and 0.7 < Pr < 160)
            Nu = 0.023 * self.Re_outer**0.8 * Pr**0.4

        elif self.cool_model == "gnielinski":
            # Gnielinski correlation (valid for 3000 < Re < 5e6 and 0.5 < Pr < 2000)
            _, fd = self.pressure_drop(
                idx,
                self.coolant_outer,
                self.D_h_outer,
                self.Re_outer,
                self.v_coolant_outer,
            )
            Nu = (
                (fd / 8)
                * (self.Re_outer - 1000)
                * Pr
                / (1 + 12.7 * (fd / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
            )

        else:
            raise ValueError(
                'Invalid heat transfer method. Select: "dittus-boelter", "dittus-boelter-simple" or "gnielinski"'
            )

        # convert Nusselt number to heat transfer coefficient
        halpha = Nu * k / self.D_h_outer[idx]

        if halpha < 0:
            raise ValueError(
                "Negative heat transfer coefficient, check applicability of cooling model to Reynolds number range"
            )

        # check if halpha is an array or a float (artifact of the pressure drop dependence or fluid model)
        if isinstance(halpha, float):
            return halpha
        else:
            return halpha[0]

    def section_heat_flux(self, idx):
        self.q_rad = self.radiation(idx)

        # solve 2D section temperature profile, passes functions for heat transfer coefficients for temperature dependence of boundary conditions

        print(
            "SOLVING Section Number	",
            len(self.geometry[:, 1]) - idx,
            " / ",
            len(self.geometry[:, 1]),
        )

        solver = HeatEquationSolver(
            idx,
            self.gas,
            self.material,
            self.cooling_geometry,
            self.heat_transfer_coefficient_gas,
            self.heat_trans_coeff_coolant_inner,
            self.heat_trans_coeff_coolant_outer,
            self.q_rad,
            self.coolant_inner.T,
            self.coolant_outer.T,
            self.T_amb,
            self.out.folder_path,
            self.settings2D,
        )
        solver.run_sim()

        self.T_wall_i = (
            solver.inner_wall_temp
        )  # inner wall temperature taken as maximum temperature of the 2D profile
        self.halpha = solver.halpha  # heat transfer coefficients taken as the of the boundary conditions of the 2D solver at the last timestep
        self.halpha_c_bottom = (
            solver.halpha_c_bottom
        )  # average heat transfer coefficient at bottom wall
        self.halpha_c_top = solver.halpha_c_top
        self.halpha_c_side1 = solver.halpha_c_side1
        self.halpha_c_side2 = solver.halpha_c_side2
        self.T_hg = solver.T_hg

        # total heat flux into the inner chamber wall
        self.q = np.mean(solver.q_chamber)

        # radiation leaving to ambient
        self.q_rad_out = np.mean(solver.q_outer)

        # ensure that FiPy output object is actually a float, multiplied by the number of channels for total area
        Q_inner = (
            float(solver.Q_c_inner) * self.section_length[idx] * self.cooling_geometry.N
        )
        Q_outer = (
            float(solver.Q_c_outer) * self.section_length[idx] * self.cooling_geometry.N
        )

        # update cooling fluid and pressure. Expected temperature rise and pressure drop
        dT_inner = Q_inner / (self.m_dot_c_inner * self.coolant_inner.Cp)
        dp_inner, _ = self.pressure_drop(
            idx, self.coolant_inner, self.D_h_inner, self.Re_inner, self.v_coolant_inner
        )
        dT_outer = Q_outer / (self.m_dot_c_outer * self.coolant_outer.Cp)
        dp_outer, _ = self.pressure_drop(
            idx, self.coolant_outer, self.D_h_outer, self.Re_outer, self.v_coolant_outer
        )

        # recalcualte cooling fluid properties, such as density, Pr, Cp etc.
        self.coolant_inner = thermo.Mixture(
            IDs=self.coolant_inner.IDs,
            ws=self.coolant_inner.ws,
            T=(self.coolant_inner.T + dT_inner),
            P=(self.coolant_inner.P - dp_inner),
        )
        self.coolant_outer = thermo.Mixture(
            IDs=self.coolant_outer.IDs,
            ws=self.coolant_outer.ws,
            T=(self.coolant_outer.T + dT_outer),
            P=(self.coolant_outer.P - dp_outer),
        )

    def run(self):
        # arrays for section length and area in between the 2D sections solved in the thermal sim
        self.section_length = np.ndarray(len(self.geometry[:, 1]))

        # determine the section length and inner chamber surface area in each section
        for i in range(len(self.geometry[:, 1])):
            # section number 0 at the injector!!!

            if self.settings2D.start_idx == -1:
                # for coolant entering at the bottom of the cooling channels, use idx as inverse of i
                idx = len(self.geometry[:, 1]) - i - 1

            elif self.settings2D.start_idx == 0:
                # in case coolant starts at injector side
                idx = i

            else:
                raise ValueError(
                    "Invalid starting point for cooling fluid, select -1 or 0 for nozzle or injector side, respectively"
                )

            if i == 0:
                self.section_length[idx] = 0

            else:
                x_len = (self.geometry[idx + 1, 0] - self.geometry[idx, 0]) ** 2
                y_len = (self.geometry[idx + 1, 1] - self.geometry[idx, 1]) ** 2

                self.section_length[idx] = np.sqrt(x_len + y_len)

            # solve section heat transfer and temperature field and update coolant properties
            self.section_heat_flux(idx)

            # write to output class
            self.out.halpha[idx] = self.halpha
            self.out.q_rad[idx] = self.q_rad
            self.out.q[idx] = self.q
            self.out.T_wall_i[idx] = self.T_wall_i
            self.out.halpha_c_bottom[idx] = self.halpha_c_bottom
            self.out.halpha_c_top[idx] = self.halpha_c_top
            self.out.halpha_c_side1[idx] = self.halpha_c_side1
            self.out.halpha_c_side2[idx] = self.halpha_c_side2
            self.out.T_c_inner[idx] = self.coolant_inner.T
            self.out.P_c_inner[idx] = self.coolant_inner.P
            self.out.T_c_outer[idx] = self.coolant_outer.T
            self.out.P_c_outer[idx] = self.coolant_outer.P
            self.out.Re_inner[idx] = self.Re_inner
            self.out.Re_outer[idx] = self.Re_outer
            self.out.T_hg[idx] = self.T_hg
            self.out.v_coolant_inner[idx] = self.v_coolant_inner
            self.out.v_coolant_outer[idx] = self.v_coolant_outer
            self.out.T_c_avg[idx] = (
                self.m_dot_c_outer * self.coolant_outer.T
                + self.m_dot_c_inner * self.coolant_inner.T
            ) / (self.m_dot_c_outer + self.m_dot_c_inner)

        # write to file in output folder
        self.out.write_csv()
