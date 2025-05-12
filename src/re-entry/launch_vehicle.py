from constants import Constants

class Launcher:
    def __init__(self,
                 name,
                 mass,
                 thrust):
        self.ct = Constants()
        self.mass = mass
        self.thrust = thrust
        self.weight = mass * self.ct.g
        self.flight_path
