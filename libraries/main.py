#import libraries.initialising as init
import libraries.iterations as itera
import libraries.matbal as mb
import plotly
import json
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit




def matbal_run(dict_tank, df_prod, dict_pvtmaster, df_pvt_oil, df_pvt_gas, regress, regress_config=None):
    #####Material Balance
    data_dict = {
        'df_prod': df_prod,
        'dict_pvtmaster': dict_pvtmaster,
        'df_pvt_oil': df_pvt_oil,
        'df_pvt_gas': df_pvt_gas,
        'dict_tank': dict_tank
    }

    df_prod = data_dict['df_prod']
    # Pres_calc, ts_obs, reservoir_pressure_obs, ts = eval_mbal_input2(data_dict)
    DDI = [None] * len(df_prod['np'])
    SDI = [None] * len(df_prod['np'])
    WDI = [None] * len(df_prod['np'])
    CDI = [None] * len(df_prod['np'])
    if regress is False:
        Pres_calc, ts_obs, reservoir_pressure_obs, ts = itera.eval_mbal_input(data_dict)
    else:
        popt, sd = mbal_fit(data_dict)
        dict_tank['initial_inplace'] = popt[0]
        dict_tank['wei'] = popt[1]
        dict_tank['J'] = popt[2]
        data_dict['dict_tank'] = dict_tank
        Pres_calc, ts_obs, reservoir_pressure_obs, ts = itera.eval_mbal_input(data_dict)

    data_dict['Pres_calc'] = Pres_calc
    DDI, SDI, WDI, CDI = drive_indices(data_dict)
    # plot = match_plot(ts, Pres_calc, ts_obs, reservoir_pressure_obs)
    return ts, Pres_calc, ts_obs, reservoir_pressure_obs, DDI, SDI, WDI, CDI, \
        dict_tank['initial_inplace'], dict_tank['wei'], dict_tank['J']


def match_plot(ts, Pres_calc, ts_obs, reservoir_pressure_obs):
    dataseries = []
    act1 = dict(
        name='Observed Data',
        # fill='tozeroy',
        # stackgaps='infer zero',
        # orientation = 'vertical',
        mode='line',
        type='scatter',
        x=Pres_calc,
        y=ts,
        # stackgroup='one',
    )
    act2 = dict(
        name='Calculated Data',
        # fill='tozeroy',
        # stackgaps='infer zero',
        # orientation = 'vertical',
        mode='line',
        type='scatter',
        x=reservoir_pressure_obs,
        y=ts_obs,
        # stackgroup='one',
    )
    dataseries.append(act1)
    dataseries.append(act2)
    graphs = [
        dict(
            data=dataseries
        )
    ]

    graphJSON = json.dumps(graphs, cls=plotly.utils.PlotlyJSONEncoder)

    # forecast_kpi = pd.DataFrame(rows_list)

    # df_json = forecast_kpi.to_json(orient='records')
    # data = {'graphJSON': graphJSON}
    return graphJSON


def mbal_fit(dict):
    def fit_mbal_input(ts_obs, N, Wei, J):
        Pres_calc2 = []
        # Pres_calc.clear()
        dict_tank = dict['dict_tank']
        dict_tank['initial_inplace'] = N
        dict_tank['wei'] = Wei
        dict_tank['J'] = J
        dict['dict_tank'] = dict_tank
        Pres_calc2, ts_obs, reservoir_pressure_obs, ts = mb.eval_mbal_input2(dict)
        Pres_calc_obs = []
        ts_obs_vals = ts_obs.values
        for x in range(len(ts_obs_vals)):
            Pres_calc_obs.append(np.interp(ts_obs_vals[x], ts, Pres_calc2))
        return Pres_calc_obs

    df_prod = dict['df_prod']
    dates = df_prod['datestamp']
    ts = pd.to_numeric(dates - df_prod['datestamp'].min()) / 864e11
    reservoir_pressure_obs = df_prod['pressure']
    ts_obs = ts[reservoir_pressure_obs.notnull()]
    reservoir_pressure_obs = reservoir_pressure_obs[reservoir_pressure_obs.notnull()]
    reservoir_pressure_obs = reservoir_pressure_obs * 1.0
    popt, pcov = curve_fit(fit_mbal_input, ts_obs, reservoir_pressure_obs, bounds=([1E6, 0.00001, 0.0001], [1E9, 10E9, 10.0]))
    sd = np.sqrt(np.diag(pcov))
    return popt, sd




# def mbal_run_simulation():
#     regress = False
#     regress_config = None
#     instance_tank = dict(request.POST.lists())
#     prod_fk = instance_tank['production_fk'][0]
#     pvt_fk = instance_tank['pvt_fk'][0]
#     # aqu_fk = instance_tank['aquifer_fk']
#
#     instance_prod = ProductionDataSet_MatBal.objects.filter(definition_fk=prod_fk)
#     instance_pvt_o = Pvt_table_oil.objects.filter(definition_table=pvt_fk)
#     instance_pvt_g = Pvt_Table_Gas.objects.filter(definition_table=pvt_fk)
#     instance_pvt_master = PvtMaster.objects.get(pk=pvt_fk)
#     # instance_aquifer = AquiferMaster.objects.get(pk=aqu_fk)
#
#     df_prod = pd.DataFrame(read_frame(instance_prod))
#     df_pvt_o = pd.DataFrame(read_frame(instance_pvt_o))
#     df_pvt_g = pd.DataFrame(read_frame(instance_pvt_g))
#     # df_tank_master = pd.DataFrame(read_frame(instance_tank))
#     # df_pvt_master = pd.DataFrame(read_frame(instance_pvt_master))
#     regress_config = {}
#     list_regress = ['initial_inplace_regress', 'aquifer_size_regress', 'aquifer_pi_regress']
#     for var in list_regress:
#         try:
#             regress_config[var] = instance_tank[var][0]
#             regress = True
#         except:
#             pass
#
#     #regress = True
#     ts, Pres_calc, ts_obs, reservoir_pressure_obs, DDI, SDI, WDI, CDI, N, Wei, J = mbal2.matbal_run2(instance_tank, df_prod,
#                             instance_pvt_master, df_pvt_o, df_pvt_g, regress, regress_config)
#
#     return ts, Pres_calc, ts_obs, reservoir_pressure_obs, DDI, SDI, WDI, CDI, N, Wei, J
#
# if __name__ == 'main':
#     mbal_run_simulation()
