import pandas as pd
import krpy


class SCAL:
    def __init__(self):
        self.scal_table = None

    def import_scal(self, filename):
        self.scal_table = pd.read_csv(filename)

    def generate_scal(self, correlation_name, **kwargs):
        self.scal_table = krpy.corey(**kwargs)