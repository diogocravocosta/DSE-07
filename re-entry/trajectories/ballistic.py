import numpy as np

class BallisticEntry:
    def __init__(self,
                 V_e,        # Entry velocity [m/s]
                 gamma_e,    # Entry flight‐path angle [rad] (negative for descent)
                 K,          # Ballistic parameter [N/m^2]
                 rho0=1.225, # Sea‐level density [kg/m^3]
                 Hs=7200.0,  # Scale height [m]
                 g=9.81      # Gravity [m/s^2]
                 ):
        """
        Parameters
        ----------
        V_e     : initial entry speed (m/s)
        gamma_e : entry angle below horizontal (rad, negative)
        K       : ballistic parameter = W / (C_D S) [N/m^2]
        rho0    : sea-level density [kg/m^3]
        Hs      : atmospheric scale height [m]
        g       : gravitational acceleration [m/s^2]
        """
        self.Ve    = V_e
        self.gamma = gamma_e
        self.K     = K
        self.rho0  = rho0
        self.Hs    = Hs
        self.g     = g
        # beta such that rho = rho0 * exp(-beta*h)
        self.beta  = 1.0 / Hs

    def density(self, h):
        """Exponential atmosphere: rho(h) = rho0 * exp(-beta * h)"""
        return self.rho0 * np.exp(-self.beta * h)

    def pressure(self, h):
        """
        Exponential pressure for isothermal atmosphere:
          p(h) = (g/beta) * rho(h) = p0 * exp(-beta h)
        """
        return (self.g / self.beta) * self.density(h)

    def V_ratio(self, h):
        """
        Normalized speed V/V_e as function of altitude h using Eq. (5.30):
          ln(V/Ve) = [1/(2 K sin(gamma))] * p(h)
        """
        p = self.pressure(h)
        exponent = p / (2.0 * self.K * np.sin(self.gamma))
        return np.exp(exponent)

    def velocity(self, h):
        """True entry speed V(h) = Ve * V_ratio(h)"""
        return self.Ve * self.V_ratio(h)

    def deceleration(self, h):
        """
        Analytical deceleration a(h) = g/(2K) * rho(h) * V(h)^2  (Eq. 5.39)
        Returns a in m/s^2.
        """
        rho = self.density(h)
        V   = self.velocity(h)
        return (self.g / (2.0 * self.K)) * rho * V**2

    def max_deceleration(self):
        """
        Calculates:
          V'/V_e = 1/sqrt(e) ≈ 0.606  (Eq. 5.46)
          a_max   = beta * (-sin gamma) * V_e^2 / 2e   (Eq. 5.47, signed positive)
          h'      = (1/beta) ln[ -rho0 g / (K beta sin gamma) ]  (Eq. 5.49)
        Returns
        -------
        Vp      : speed at max decel [m/s]
        apos    : max decel [m/s^2]
        h_prime : altitude of max decel [m]
        """
        # normalised velocity
        V_ratio_p = 1.0 / np.sqrt(np.e)               # Eq. (5.46)
        Vp        = self.Ve * V_ratio_p

        # maximum deceleration (positive value)
        a_max = - (self.beta * np.sin(self.gamma) * self.Ve**2) / (2.0 * np.e)
        # altitude of maximum deceleration
        arg = - self.rho0 * self.g / (self.K * self.beta * np.sin(self.gamma))
        if arg <= 0:
            h_p = np.nan
        else:
            h_p = (1.0 / self.beta) * np.log(arg)

        return Vp, a_max, h_p

# --------------------------------------------------------------------
# Example usage: Apollo‐like entry
# --------------------------------------------------------------------
if __name__ == "__main__":
    # Apollo CM from LEO:
    V_e    = 7_800.0        # m/s
    gamma  = np.deg2rad(-10)  # entry angle -10°
    K      = 2806.0         # N/m^2  (example from text)
    rho0   = 1.225          # kg/m^3
    Hs     = 7200.0         # m
    g      = 9.81           # m/s^2

    entry = BallisticEntry(V_e, gamma, K, rho0=rho0, Hs=Hs, g=g)

    # sample profile
    hs = np.linspace(0, 80_000, 201)   # from sea‐level up to 80 km
    Vs = entry.velocity(hs)
    As = entry.deceleration(hs) / g    # in g‐units

    # compute max deceleration parameters
    Vp, a_max, h_p = entry.max_deceleration()

    print(f"Max deceleration a_max = {a_max/g: .1f} g at h' = {h_p/1e3: .1f} km")
    print(f"Speed at max decel: V' = {Vp/1e3: .2f} km/s (ratio {Vp/V_e:.3f})")

    # (Optional) quick plot
    try:
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
        ax1.plot(Vs/V_e, hs/1e3)
        ax1.set_xlabel("V/V_e")
        ax1.set_ylabel("Altitude h [km]")
        ax1.set_title("Analytical h-V profile")
        ax1.grid(True)

        ax2.plot(As, hs/1e3)
        ax2.set_xlabel("a/g")
        ax2.set_ylabel("Altitude h [km]")
        ax2.set_title("Analytical deceleration")
        ax2.axhline(h_p/1e3, color='k', ls='--')
        ax2.grid(True)

        plt.tight_layout()
        plt.show()
    except ImportError:
        pass
