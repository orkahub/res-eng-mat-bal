from pvtpy import CorrA
import pandas as pd

class PVT:
    def __init__(self):
        self.pvt_table = None

    def import_pvt(self, filename):
        self.pvt_table = pd.read_csv(filename)

    def generate_pvt(self, correlation_name, **kwargs):
        self.pvt_table = CorrA.formation_total_volume_factor(**kwargs)
