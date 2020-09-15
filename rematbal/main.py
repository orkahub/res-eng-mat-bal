import rematbal.iterations as itera
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


class tank():

    def __init__(self):
        self.prod_table = pd.DataFrame()
        self.oil_pvt_table = pd.DataFrame()
        self.gas_pvt_table = pd.DataFrame()
        self.regress = False
        self.regress_config = None
        self.tank_data = {
            'initial_inplace': None,
            'initial_gascap': None,
            'initial_pressure': None,
            'wei': None,
            'J': None,
            'swi': None,
            'cw': None,
            'cf': None,
            'Boi': None,
            'Bgi': None
        }
        self.pvt_master = {
            'gor': None,
            'sat_press': None,
            'temperature': None,
        }
        self.ts_results = pd.DataFrame()

    def matbal_run(self):
        #####Material Balance

        def mbal_fit(dict):
            def fit_mbal_input(ts_obs, N, Wei, J):
                Pres_calc2 = []
                # Pres_calc.clear()
                dict_tank = dict['dict_tank']
                dict_tank['initial_inplace'] = N
                dict_tank['wei'] = Wei
                dict_tank['J'] = J
                dict['dict_tank'] = dict_tank
                Pres_calc2, ts_obs, reservoir_pressure_obs, ts = itera.eval_mbal_input(dict)
                Pres_calc_obs = []
                ts_obs_vals = ts_obs.values
                for x in range(len(ts_obs_vals)):
                    Pres_calc_obs.append(np.interp(ts_obs_vals[x], ts, Pres_calc2))
                return Pres_calc_obs

            BOUNDS = ([1E6, 0.00001, 0.0001], [1E9, 10E9, 10.0])

            df_prod = dict['df_prod']
            dates = df_prod['datestamp']
            ts = pd.to_numeric(dates - dates.min()) / 864e11
            reservoir_pressure_obs = df_prod['pressure']
            ts_obs = ts[reservoir_pressure_obs.notnull()]
            reservoir_pressure_obs = reservoir_pressure_obs[reservoir_pressure_obs.notnull()]
            reservoir_pressure_obs = reservoir_pressure_obs * 1.0
            popt, pcov = curve_fit(fit_mbal_input, ts_obs, reservoir_pressure_obs,
                                   bounds=BOUNDS)
            sd = np.sqrt(np.diag(pcov))
            return popt, sd

        data_dict = {
            'df_prod': self.prod_table,
            'dict_pvtmaster': self.pvt_master,
            'df_pvt_oil': self.oil_pvt_table,
            'df_pvt_gas': self.gas_pvt_table,
            'dict_tank': self.tank_data
        }

        initial_container = [None] * len(self.prod_table['np'])
        ddi = sdi = wdi = cdi = initial_container

        if self.regress is False:
            solution_set = itera.eval_mbal_input(data_dict)
        else:
            popt, sd = mbal_fit(data_dict)
            self.tank_data['initial_inplace'] = popt[0]
            self.tank_data['wei'] = popt[1]
            self.tank_data['J'] = popt[2]
            data_dict['dict_tank'] = self.tank_data
            solution_set = itera.eval_mbal_input(data_dict)

        #data_dict['Pres_calc'] = pres_calc
        #ddi, sdi, wdi, cdi = self.drive_indices()
        # plot = match_plot(ts, Pres_calc, ts_obs, reservoir_pressure_obs)
        self.ts_results['Time'] = solution_set['ts']
        self.ts_results['Calculated Pressure'] = solution_set['pres_calc']
        self.ts_results['Depletion Drive Index'] = solution_set['ddi']
        self.ts_results['Segregation Drive Index'] = solution_set['sdi']
        self.ts_results['Water Drive Index'] = solution_set['wdi']
        self.ts_results['Compaction Drive Index'] = solution_set['cdi']
        self.ts_results['Aquifer Water Influx'] = solution_set['we']
        self.ts_results['Aquifer Pressure'] = solution_set['aquifer_pres']
        self.ts_results['oip'] = solution_set['oip']
        return self.ts_results, self.tank_data


