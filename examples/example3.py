# insert directory where libraries exits
import sys
sys.path.insert(1,'./')
import plotly.graph_objs as go
import pandas as pd
import os
import datetime
import rematbal.plots as plots
from plotly.subplots import make_subplots
from plotly.offline import plot
from rematbal.main import tank as tank_instance
from scipy import stats
import matplotlib.pyplot as plt


regress = False
regress_config = None

path = os.path.dirname(__name__)
file = os.path.join(path, '../data/lsu_matbal/OilMBEx.xls')
#keep_cols = ['Days', 'Np\nSTBO', 'Gp\nMCF', 'Wp\nSTBW', 'Meas P\npsia']
df_prod = pd.read_excel(file, sheet_name='Calculations', skiprows=12, nrows=73)
df_prod.rename(columns={'Np\nSTBO': 'np', 'Gp\nMCF': 'gp', 'Wp\nSTBW': 'wp', 'Meas P\npsia': 'pressure'}, inplace=True)
startdate = datetime.date(2010, 1, 1)
df_prod['datestamp'] = [datetime.timedelta(day)+startdate for day in df_prod['Days']]
df_prod['wi'] = 0
df_prod['gi'] = 0
df_prod['gp'] = df_prod['gp']*1000.0
#keep_cols = ['Pressure\npsia', 'Bo\nrb/stbo', 'Rs\nscf/stbo', 'Bg rb/mcf', 'Bt rb/stbo']
df_pvt_o = pd.read_excel(file, sheet_name='Oil PVT', skiprows=3, nrows=16)
df_pvt_o.rename(columns={'Pressure\npsia': 'pressure', 'Bo\nrb/stbo': 'oil_fvf', 'Rs\nscf/stbo': 'solution_gas'}, inplace=True)
#keep_cols = ['Pressure psia', 'z', 'Bg\nrb/mcf']
df_pvt_g = pd.read_excel(file, sheet_name='Z Factors', skiprows=3, nrows=12)
df_pvt_g.rename(columns={'Pressure\npsia': 'pressure', 'Bg\nrb/mcf': 'gas_fvf'}, inplace=True)
tank = {
    'initial_inplace': 13.8E6,
    'initial_gascap': 0,
    'initial_pressure': 10180,
    'wei': 141.8e6,
    'J': 0.93,
    'swi': 0.2,
    'cw': 2.5e-6,
    'cf': 3e-5,
    'Boi': 1.735,
    'Bgi': 0.6508
}
pvt_master = {
    'gor': 1720,
    'sat_press': 8227,
    'temperature': 219,
}

tank = {'initial_inplace': 5393487.759211203, 'initial_gascap': 0, 'initial_pressure': 10180, 'wei': 10999999999.9988,
        'J': 0.4446624977502875, 'swi': 0.2, 'cw': 2.5e-06, 'cf': 3e-05, 'Boi': 1.735, 'Bgi': 0.6508,
        'initial_inplace_sd': 141439.29428151343, 'wei_sd': 0.008334858350091836, 'J_sd': 0.012664217123234959}

tank1 = tank_instance()
tank1.tank_data = tank
tank1.prod_table = df_prod
tank1.oil_pvt_table = df_pvt_o
tank1.gas_pvt_table = df_pvt_g
tank1.pvt_master = pvt_master
#tank1.regress = True
#tank1.regress_config = None
#ts_res, tank_results = tank1.matbal_run()
#print(tank_results)



sampling = 100
binning = 20

plt.figure(figsize=(15, 15))
STOIIPs = stats.norm.rvs(size=sampling, loc=tank['initial_inplace'], scale=tank['initial_inplace_sd'])
Weis = stats.norm.rvs(size=sampling, loc=tank['wei'], scale=tank['wei_sd'])
Weis[Weis < 0] = 0
Js = stats.norm.rvs(size=sampling, loc=tank['J'], scale=tank['J_sd'])

tank1.regress = False

for x in range(sampling):
    tank['initial_inplace'] = STOIIPs[x]
    tank['wei'] = Weis[x]
    tank['J'] = Js[x]
    tank1.tank_data = tank
    ts_res, tank_results = tank1.matbal_run()
    plt.plot(ts_res['Time'], ts_res['Calculated Pressure'], '-')
plt.show()



