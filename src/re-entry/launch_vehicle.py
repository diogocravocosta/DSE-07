class Launcher:
    def __init__(self,
                 name,
                 mass,
                 thrust):
        self.mass = mass
        self.thrust = thrust
        self.weight = mass * g
