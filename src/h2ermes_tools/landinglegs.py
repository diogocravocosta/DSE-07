from data import (
    constants,
    material
)

import numpy as np


class LandingLegs:
    class Vehicle:
        """
        Internal class to handle the vehicle properties. Stores mass without legs, center of gravity, and geometry.
        """

        def __init__(self, mass_land: float, phi: float, r_bottom: float) -> None:
            """
            Internal vehicle class constructor. For X-Y plane, coordinate system defined as X along the length of
            the vehicle.
            Args:
                mass_land (float): Landing mass of the vehicle (without landing legs) [kg]
                x_cog (float): X location of the center of the gravity [m]
                y_cog (float): Y location of the center of the gravity [m]
                phi (float): Phi angle of the vehicle [deg]
                r_bottom (float): Radius of the bottom of the vehicle [m]
            """
            self.mass_land = mass_land
            self.phi = phi
            self.r_bottom = r_bottom

    class SafetyFactors:
        """
        Internal class to handle safety factors. Contains safety factors for:
            - Metal (general)
            - Joints
            - Buckling
        For all three, both yield and ultimate safety factors are present.
        """

        class FailureMode:
            """
            Internal class for the safety factors to handle different failure modes.
            """

            def __init__(self, y: float, u: float) -> None:
                """
                Args:
                    y: Yield safety factor [-]
                    u: Ultimate safety factor [-]
                """
                self.y = y
                self.u = u

        def __init__(self) -> None:
            """
            Constructor for safety factors. Single safety factors can be called as SafetyFactors.metal.y as an example.
            """
            self.metal = self.FailureMode(1.5, 2.0)
            self.joint = self.FailureMode(2.0, 2.0)
            self.buckling = self.FailureMode(1.5, 2.0)

    class ShockAbsorber:
        """
        Internal class for the shock absorber. Oleo-pneumatic design.
        """

        def __init__(self, ll) -> None:
            self.ll = ll

            self._stroke_length()

            self.outer_radius = self.ll.radius - 0.001
            self.bore_thickness = 0.002
            self._piston_area()

            self._gas_volume()

            self.rod_radius = self.outer_radius - self.bore_thickness - 0.002
            self.rod_length = self.stroke_length * 1.1  # overlap of 10%
            self.rod_thickness = 0.002

            self.barrel_length = self.stroke_length + 0.5 + self.gas_volume / self.piston_area

            self._mass()

        def _gas_volume(self) -> float:
            """
            Calculates the gas volume
            Returns: Gas volume [m^3]
            """
            p1 = self.ll.mass_land / self.ll.n_legs / self.piston_area
            p2 = 2 * p1

            volume_fraction = (p1 / p2) ** 1.4  # pv^gamma = constant

            delta_volume = self.stroke_length * self.piston_area
            self.gas_volume = delta_volume / (1 - volume_fraction)
            return self.gas_volume

        def _stroke_length(self) -> float:
            """
            Calculates the stroke length required for the shock absorber.
            Returns: Stroke length [m]
            """
            mass = self.ll.mass_land / self.ll.n_legs
            a_max = 2  # g
            v_land = 5  # m/s
            kinetic_energy = 0.5 * mass * v_land ** 2 * 2  # safety factor of 2
            force_average = a_max * mass

            self.stroke_length = kinetic_energy / force_average
            return self.stroke_length

        def _piston_area(self) -> float:
            """
            Calculates the piston area required for the shock absorber.
            Returns: Piston area [m^2]
            """
            self.piston_area = np.pi * self.outer_radius ** 2
            return self.piston_area

        def _mass(self) -> float:
            """
            Calculates the mass of the shock absorber.
            Returns: Mass [kg]
            """
            gas_mass = self.gas_volume * 1.25
            piston_mass = np.pi * self.rod_radius * self.rod_thickness * self.rod_length * self.ll.material.rho
            barrel_mass = self.barrel_length * np.pi * self.outer_radius * self.bore_thickness * self.ll.material.rho
            self.mass = gas_mass + piston_mass + barrel_mass + 2.5
            return self.mass

    class Deployment:
        """
        Internal class for the pneumatic deployment mechanism.
        """

        def __init__(self, ll) -> None:
            self.ll = ll

            self._actuator()

            self._mass()

        def _actuator(self) -> None:
            """
            Sizes the actuator.
            """
            self.radius = self.ll.radius - 0.001
            self.length = self.ll.length + self.ll.shockabsorber.stroke_length
            self.thickness = 0.003

            self.rod_radius = self.radius - self.thickness - 0.002
            self.rod_thickness = self.thickness

        def _mass(self) -> float:
            """
            Calculates the mass of the deployment mechanism.
            Returns: Mass [kg]
            """

            tank_mass = 2.5
            rod_mass = np.pi * self.rod_radius * self.rod_thickness * self.length * self.ll.material.rho
            barrel_mass = np.pi * self.radius * self.thickness * self.length * self.ll.material.rho

            self.mass = tank_mass + rod_mass + barrel_mass + 3  # 3 for locking, mounts, etc.
            return self.mass

    class AeroCover:
        """
        Internal class for the aerocover of the landing legs. Mainly used for sizing as a dependency of the legs.
        """

        def __init__(self, ll) -> None:
            self.ll = ll

            self.length = 1.1 * self.ll.length
            self.radius_a = 2.5 * self.ll.radius
            self.radius_b = 1.5 * self.ll.radius
            self.thickness = 0.001
            self.material = self.ll.material
            self._mass()

        def _mass(self) -> float:
            area = np.pi / 2 * (self.radius_a * self.radius_b - (self.radius_b - self.thickness) * (
                    self.radius_a - self.thickness))
            volume = area * self.length
            self.mass = volume * self.material.rho
            return self.mass

    def __init__(self, n_legs: int, mass_land: float, phi: float, r_bottom: float,
                 material: material.Material, clearance_height: float) -> None:
        """
        Landing legs constructor.
        Args:
            n_legs (int): Number of legs [-]
            mass_land (float): Landing mass of the vehicle (without landing legs) [kg]
            x_cog (float): X location of the center of the gravity [m]
            y_cog (float): Y location of the center of the gravity [m]
            phi (float): Phi angle of the vehicle [deg]
            r_bottom (float): Radius of the bottom of the vehicle [m]
            material (Material): Material selected for the landing legs
            clearance_height (float): Clearance height required by the landing legs, defined as the distance between the
                                      bottom point and the most bottom point of the vehicle (on the TPS most likely) [m]
        """
        # Define number of legs, vehicle class, safety factors, material, and clearance height
        self.n_legs = n_legs
        self.v = self.Vehicle(mass_land, phi, r_bottom)
        self.sf = self.SafetyFactors()
        self.material = material
        self.h_clearance = clearance_height

        # Computes the minimum length of the legs based on the required clearance and the phi angle
        self.min_length = self.h_clearance / np.cos(np.radians(self.v.phi))

        # Initialize geometry parameters
        self.radius = 0
        self.length = 0
        self.thickness = 0
        self.area = self._area()
        self.ixx = self._area_moment_inertia()
        self.mass_land = self.v.mass_land + self.n_legs * self._mass()

        # Initialize failure modes
        self.buckling_failure = False
        self.yield_failure = False
        self.interaction_failure = False
        self.bending_failure = False
        self.compressive_failure = False

        # Initialize possible configurations
        self.configs = []

    def __reset_failures(self) -> None:
        """
        Private method to reset the failure modes for a new configuration.
        """
        self.buckling_failure = False
        self.yield_failure = False
        self.interaction_failure = False
        self.bending_failure = False
        self.compressive_failure = False

    def _set_geometry(self, length: float, radius: float, thickness: float) -> None:
        """
        Internal method to set the geometry parameters.
        Args:
            length (float): Length of the legs [m]
            radius (float): Radius of the legs [m]
            thickness (float): Thickness of the legs [m]
        """
        self.length = length
        self.radius = radius
        self.thickness = thickness

    def _area(self) -> float:
        """
        Internal method to calculate the cross-sectional area.
        Returns: Cross-sectional area [m^2]
        """
        self.area = np.pi * (self.radius ** 2 - (self.radius - self.thickness) ** 2)
        return self.area

    def _area_moment_inertia(self) -> float:
        """
        Internal method to calculate the second moment of inertia.
        Returns: Second moment of inertia [m^4]
        """
        self.ixx = np.pi / 4 * (self.radius ** 4 - (self.radius - self.thickness) ** 4)
        return self.ixx

    def _mass(self) -> float:
        """
        Internal method to calculate the mass.
        Returns: Mass of a single leg [kg]
        """
        self.mass = self._area() * self.length * self.material.rho
        return self.mass

    def _determine_loads(self):
        """
        Internal method to determine the different loads based on the geometry and mass
        """
        self.force = self.mass_land * constants.g_0 / self.n_legs
        self.force_transverse = self.force * np.sin(np.radians(self.v.phi))
        self.force_axial = self.force * np.cos(np.radians(self.v.phi))

    def _euler_buckling(self) -> float:
        """
        Internal method to calculate the critical Euler buckling stress.
        Returns: Critical Euler buckling stress, sigma_cr [Pa]
        """
        # k = 1  # Pinned-pinned connection; conservative approach
        # self.stress_critical_euler_buckling = np.pi ** 2 * self.material.E * self._area_moment_inertia() / (
        # (k * self.length) ** 2 * self._area())
        self.stress_critical_euler_buckling = np.pi ** 2 * self.material.E / self.slenderness_ratio ** 2
        return self.stress_critical_euler_buckling

    def _johnson_buckling(self) -> float:
        """
        Internal method to calculate the critical Johnson buckling stress.
        Returns: Critical Johnson buckling stress, sigma_cr [Pa]
        """
        self.stress_critical_johnson_buckling = self.material.ys - 1 / self.material.E * (
                self.material.ys / (2 * np.pi)) ** 2 * self._slenderness_ratio() ** 2
        return self.stress_critical_johnson_buckling

    def _slenderness_ratio(self) -> float:
        """
        Internal method to calculate the slenderness ratio for buckling analysis. Also determines the critical
        slenderness ratio.
        Returns: Slenderness ratio [-]
        """
        k = 1
        self.slenderness_ratio = k * self.length * np.sqrt(self._area() / self._area_moment_inertia())
        self.slenderness_ratio_critical = np.sqrt(2 * np.pi ** 2 * self.material.E / self.material.ys)
        return self.slenderness_ratio

    def _buckling_analysis(self) -> None:
        """
        Internal method to perform the buckling analysis.
        Determines which buckling to use based on the critical slenderness ratio.
        """
        self._slenderness_ratio()

        # Johnson buckling better for lower slenderness ratio, else Euler
        if self.slenderness_ratio < self.slenderness_ratio_critical:
            self.stress_critical_buckling = self._johnson_buckling()
        else:
            self.stress_critical_buckling = self._euler_buckling()

        # Set the failure to true in case the leg buckles
        if self.sf.buckling.y * self.stress_axial >= self.stress_critical_buckling:
            self.buckling_failure = True

    def _compressive_analysis(self) -> None:
        """
        Internal method to perform the compressive analysis.
        """
        self.stress_axial = self.force_axial / self._area()

        # Set the failure to true in case the leg yields under pure compression
        if self.sf.metal.y * self.stress_axial > self.material.ys:
            self.compressive_failure = True

    def _bending_analysis(self) -> None:
        """
        Internal method to perform the bending analysis. Uses the formula:
        sigma = (F*l*r) / I,
        where F is the transverse force (m*g*sin(phi)), l the length, r the radius, and I the area moment of inertia.
        """
        self.stress_bending = self.force_transverse * self.length * self.radius / self._area_moment_inertia()

        # Set the failure to true in case the leg yields under pure bending
        if self.sf.metal.y * self.stress_bending > self.material.ys:
            self.bending_failure = True

    def _interaction_analysis(self) -> None:
        """
        Internal method to perform the interaction analysis. Ideally, the interaction equation is:
        sigma_compressive / sigma_cr + sigma_bending / sigma_yield <= 1,

        with the first ratio for buckling and the second for bending. In case of a non-slender rod however, buckling
        can be assumed to be not limiting, and the interaction equation can better be:
        sigma_compressive / sigma_yield + sigma_bending / sigma_yield <= 1
        Returns:

        """
        self.interaction_value_yield = self.sf.metal.y * (self.stress_axial + self.stress_bending) / self.material.ys

        self.interaction_value_buckling = self.sf.buckling.y * self.stress_axial / self.stress_critical_buckling + self.sf.metal.y * self.stress_bending / self.material.ys

        # Set the failure to true in case the leg fails in interaction
        if (self.interaction_value_buckling >= 1) or (self.interaction_value_buckling >= 1):
            self.interaction_failure = True

    def _determine_failed(self) -> bool:
        """
        Internal method to evaluate every failure mode and see whether the leg failed any.
        Returns: True if leg failed, False otherwise
        """
        if sum([
            self.buckling_failure,
            self.yield_failure,
            self.interaction_failure,
            self.bending_failure,
            self.compressive_failure
        ]) > 0:
            self.failed = True
        else:
            self.failed = False

        return self.failed

    def _single_leg_landing(self) -> None:
        stress = self.mass_land * constants.g_0 / self._area() * self.sf.metal.y

        if stress > self.material.ys:
            self.compressive_failure = True

        self._slenderness_ratio()

        # Johnson buckling better for lower slenderness ratio, else Euler
        if self.slenderness_ratio < self.slenderness_ratio_critical:
            stress_critical_buckling = self._johnson_buckling()
        else:
            stress_critical_buckling = self._euler_buckling()

        # Set the failure to true in case the leg buckles
        if self.sf.buckling.y * stress >= stress_critical_buckling:
            self.buckling_failure = True

    def _failure_analysis(self):
        """
        Wrapper for all failure mode analysis methods.
        """
        self._area()
        self._area_moment_inertia()
        self._mass()

        self._determine_loads()

        self._compressive_analysis()
        self._bending_analysis()

        self._buckling_analysis()

        self._interaction_analysis()

        self._single_leg_landing()

        self._determine_failed()

    def run_sizing(self) -> None:
        """
        Method to run the sizing of the legs. Procedure:
            1) Determine a range for the three geometry parameters
                a) Length: From minimum calculated length until 5 m, steps of 0.1 m
                b) Radius: From 0.01 m until 1 m in steps of 0.01 m
                c) Thickness: From 1 mm until 5 cm in steps of 1 mm
            2) For every possible combination, analyze the leg with the loads
                a) For combinations where the thickness > radius, the analysis is not performed
            3) In case the leg is able to work, it is added to a list of possible configurations
            4) From the possible configs, the best one is selected as the lowest mass one
            5) The properties of this configuration are applied in case the stresses need to be printed
        """
        # Ranges for the geometric properties
        length_range = np.arange(round(self.min_length, 1) + 0.1, 5.1, 0.1)
        thickness_range = np.arange(0.001, 0.051, 0.001)
        radius_range = np.arange(0.01, 1.01, 0.01)

        # Iterate over each combination
        for length in length_range:
            for radius in radius_range:
                for thickness in thickness_range:
                    if thickness > radius:
                        # Impossible configuration, skipping
                        continue

                    # Perform the analysis for the configuration
                    self._set_geometry(length, radius, thickness)
                    self._failure_analysis()

                    # In case of no failure, save the configuration (geometry and mass)
                    if not self.failed:
                        self.configs.append({
                            'length': self.length,
                            'radius': self.radius,
                            'thickness': self.thickness,
                            'mass': self.mass
                        })

                    # Reset the failure attributes for the next configuration
                    self.__reset_failures()

        # Once all configs are analyzed and only the working ones are filtered out, we can get the lowest mass one
        self.best_config = min(self.configs, key=lambda config: config['mass'])

        # Set the properties of this config
        self._set_geometry(self.best_config['length'], self.best_config['radius'], self.best_config['thickness'])

        # Re-run the failure analysis to get correct loads, stresses, and interaction values for the config
        self._failure_analysis()

        # Create the shock absorber
        self.shockabsorber = self.ShockAbsorber(self)

        # Create the deployment mechanism
        self.deployment = self.Deployment(self)

        # Create the aerocover
        self.aerocover = self.AeroCover(self)

        # Update the total mass
        self.total_mass = self.n_legs * (
                self.mass + self.aerocover.mass + self.shockabsorber.mass + self.deployment.mass)


def size_landing_legs(n_legs: int, mass_land: float, phi: float, r_bottom: float, material: material.Material,
                      clearance_height: float) -> float:

    ll = LandingLegs(n_legs, mass_land, phi, r_bottom, material, clearance_height)
    ll.run_sizing()

    return ll.total_mass


if __name__ == '__main__':
    testmaterial = material.Material(youngs_modulus=70e9, yield_strength=345e6, density=2800)

    ll = LandingLegs(
        n_legs=4,
        mass_land=25000,
        phi=10,
        r_bottom=5,
        material=testmaterial,
        clearance_height=2.5
    )
    ll.run_sizing()

    print(f'BEST CONFIG ({4} legs)\n\n'
          f'L: {ll.length} m\n'
          f'R: {ll.radius} m\n'
          f't: {ll.thickness * 1000} mm\n'
          f'MASS: {ll.mass} kg')

    print(f'Total mass landing system: {ll.total_mass} ({ll.total_mass / ll.v.mass_land * 100:.2f}%)')
    print(f'Expected mass landing system (F9): {ll.v.mass_land * 0.085} (8.5%)')
