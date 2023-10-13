import pandas as pd
import tank as tank

class MatBalProject:
    def __init__(self, name, tank_properties=[]):
        self.name = name
        self.tank_properties = tank_properties
        self.tanks = []

    def initialize_tanks(self):
        for tank_prop in self.tank_properties:
            self.tanks.append(self.get_tank(tank_prop))

    def get_tank(self, tank_properties):
        return tank.tank(**tank_properties)
