# insert directory where libraries exits
import sys
sys.path.insert(1,'./')

import libraries.main as main
import pandas as pd
import os
import datetime
import libraries.plots as plots
from plotly.subplots import make_subplots
from plotly.offline import plot
import webbrowser


regress = False
regress_config = None

path = os.path.dirname(__name__)
#file = os.path.join(path, '../data/lsu_matbal/OilMBEx.xls')
file = os.path.join(path, './data/lsu_matbal/OilMBEx.xls') # Linux single dot(.)

keep_cols = ['Days', 'Np\nSTBO', 'Gp\nMCF', 'Wp\nSTBW', 'Meas P\npsia']
df_prod = pd.read_excel(file, sheet_name='Calculations', skiprows=12, nrows=73)
df_prod.rename(columns={'Np\nSTBO': 'np', 'Gp\nMCF': 'gp', 'Wp\nSTBW': 'wp', 'Meas P\npsia': 'pressure'}, inplace=True)
startdate = datetime.date(2010, 1, 1)
df_prod['datestamp'] = [datetime.timedelta(day)+startdate for day in df_prod['Days']]
df_prod['wi'] = 0
df_prod['gi'] = 0
df_prod['gp'] = df_prod['gp']*1000.0
keep_cols = ['Pressure\npsia', 'Bo\nrb/stbo', 'Rs\nscf/stbo', 'Bg rb/mcf', 'Bt rb/stbo']
df_pvt_o = pd.read_excel(file, sheet_name='Oil PVT', skiprows=3, nrows=16)
df_pvt_o.rename(columns={'Pressure\npsia': 'pressure', 'Bo\nrb/stbo': 'oil_fvf', 'Rs\nscf/stbo': 'solution_gas'}, inplace=True)
keep_cols = ['Pressure psia', 'z', 'Bg\nrb/mcf']
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

ts, Pres_calc, ts_obs, reservoir_pressure_obs, DDI, SDI, WDI, CDI, N, Wei, J = main.matbal_run(tank, df_prod,
                        pvt_master, df_pvt_o, df_pvt_g, regress, regress_config)


# Dashboard
html_graphs=open("dashboard_ex1.html", 'w')
html_graphs.write("<html><head></head><body>"+"\n")
plot1 = plots.plot_pressure_match(ts, Pres_calc, ts_obs, reservoir_pressure_obs)
plot(plot1, filename='plot1.html', auto_open=False)
html_graphs.write("  <object data=\""+'plot1.html'+"\" width=\"650\" height=\"500\"></object>"+"\n")
plot2 = plots.plot_drive_indices(ts, DDI, SDI, WDI, CDI)
plot(plot2, filename='plot2.html', auto_open=False)
html_graphs.write("  <object data=\""+'plot2.html'+"\" width=\"650\" height=\"500\"></object>"+"\n")
html_graphs.write("</body></html>")
print('Open exported html dashboard')



