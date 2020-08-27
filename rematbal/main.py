import rematbal.iterations as itera
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


class tank():
    prod_table = pd.DataFrame()
    oil_pvt_table = pd.DataFrame()
    gas_pvt_table = pd.DataFrame()
    regress = False
    regress_config = None
    tank_data = {
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
    pvt_master = {
        'gor': None,
        'sat_press': None,
        'temperature': None,
    }
    ts_results = pd.DataFrame()

    def matbal_run(self):
        #####Material Balance

        data_dict = {
            'df_prod': self.prod_table,
            'dict_pvtmaster': self.pvt_master,
            'df_pvt_oil': self.oil_pvt_table,
            'df_pvt_gas': self.gas_pvt_table,
            'dict_tank': self.tank_data
        }

        dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
        Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
        ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type = itera.mbal_setup(data_dict)

        # df_prod = dict['df_prod']
        # dates = df_prod['datestamp']
        # ts = pd.to_numeric(dates - dates.min()) / 864e11
        reservoir_pressure_obs = self.prod_table['pressure']
        # ts_obs = dates - dates.min()

        ts = dates - dates.min()
        ts = [time.days for time in ts]
        ts_obs = pd.Series(ts)
        ts_obs = ts_obs[reservoir_pressure_obs.notnull()]
        reservoir_pressure_obs = reservoir_pressure_obs[reservoir_pressure_obs.notnull()]
        reservoir_pressure_obs = reservoir_pressure_obs * 1.0

        def mbal_fit(dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb,
                                    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, ts, VEH_aq_type,
                                    td_array, VEH_dp_array, r, rr, aq_type):
            tsx = ts


            def fit_mbal_input(ts_obs, N, Wei, J):
                #Pres_calc2 = []
                # Pres_calc.clear()
                #dict_tank = dict['dict_tank']
                #dict_tank['initial_inplace'] = N
                #dict_tank['wei'] = Wei
                #dict_tank['J'] = J
                #dict['dict_tank'] = dict_tank
                #Pres_calc2, ts_obs, reservoir_pressure_obs, ts
                #solution_set = itera.eval_mbal_input(dict)
                solution_set = itera.eval_mbal_input(dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo,
                                                     pvt_oil_Bg,
                                                     pvt_oil_Rs, Rsb, Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj,
                                                     Bwinj, Ginj, J, tsx, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)

                Pres_calc_obs = []
                ts_obs_vals = ts_obs.values
                    #solution_set['ts'].values
                Pres_calc2 = solution_set['pres_calc'].values
                ts = solution_set['ts'].values
                #ts_obs.values
                for x in range(len(ts_obs_vals)):
                    Pres_calc_obs.append(np.interp(ts_obs_vals[x], ts, Pres_calc2))
                return Pres_calc_obs


            #reservoir_pressure_obs = reservoir_pressure_obs[reservoir_pressure_obs.notnull()]
            #reservoir_pressure_obs = reservoir_pressure_obs * 1.0
            popt, pcov = curve_fit(fit_mbal_input, ts_obs, reservoir_pressure_obs,
                                   bounds=([1E6, 0.01, 0.01], [1E9, 11E9, 10.0]))
            sd = np.sqrt(np.diag(pcov))
            return popt, sd




        ddi = [None] * len(self.prod_table['np'])
        sdi = [None] * len(self.prod_table['np'])
        wdi = [None] * len(self.prod_table['np'])
        cdi = [None] * len(self.prod_table['np'])
        if self.regress is False:
            solution_set = itera.eval_mbal_input(dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg,
                                    pvt_oil_Rs, Rsb, Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
                                    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
        else:
            popt, sd = mbal_fit(dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
                                    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
                                    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
            self.tank_data['initial_inplace'] = popt[0]
            N = popt[0]
            self.tank_data['wei']= popt[1]
            Wei = popt[1]
            self.tank_data['J'] = popt[2]
            J = popt[2]
            self.tank_data['initial_inplace_sd'] = sd[0]
            self.tank_data['wei_sd'] = sd[1]
            self.tank_data['J_sd'] = sd[2]
            data_dict['dict_tank'] = self.tank_data
            solution_set = itera.eval_mbal_input(dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, \
                                    pvt_oil_Bg, pvt_oil_Rs, Rsb, \
                                    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
                                    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)

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


