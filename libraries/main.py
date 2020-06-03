import libraries.iterations as itera
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

            df_prod = dict['df_prod']
            dates = df_prod['datestamp']
            ts = pd.to_numeric(dates - dates.min()) / 864e11
            reservoir_pressure_obs = df_prod['pressure']
            ts_obs = ts[reservoir_pressure_obs.notnull()]
            reservoir_pressure_obs = reservoir_pressure_obs[reservoir_pressure_obs.notnull()]
            reservoir_pressure_obs = reservoir_pressure_obs * 1.0
            popt, pcov = curve_fit(fit_mbal_input, ts_obs, reservoir_pressure_obs,
                                   bounds=([1E6, 0.00001, 0.0001], [1E9, 10E9, 10.0]))
            sd = np.sqrt(np.diag(pcov))
            return popt, sd

        data_dict = {
            'df_prod': self.prod_table,
            'dict_pvtmaster': self.pvt_master,
            'df_pvt_oil': self.oil_pvt_table,
            'df_pvt_gas': self.gas_pvt_table,
            'dict_tank': self.tank_data
        }

        ddi = [None] * len(self.prod_table['np'])
        sdi = [None] * len(self.prod_table['np'])
        wdi = [None] * len(self.prod_table['np'])
        cdi = [None] * len(self.prod_table['np'])
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

    #
    # def drive_indices(self):
    #     dates = self.prod_table['datestamp']
    #     ts = pd.to_numeric(dates - dates.min()) / 864e11
    #     np = self.prod_table['np']
    #     gp = self.prod_table['gp']
    #     wp = self.prod_table['wp']
    #     reservoir_pressure_obs = self.prod_table['pressure']
    #     ts_obs = ts[reservoir_pressure_obs.notnull()]
    #     pres_calc = dict['Pres_calc']
    #     reservoir_pressure_obs = reservoir_pressure_obs[reservoir_pressure_obs.notnull()]
    #     reservoir_pressure_obs = reservoir_pressure_obs * 1.0
    #     n = float(self.tank_data['initial_inplace'])
    #     wei = float(self.tank_data['wei'])
    #     swi = float(self.tank_data['swi'])
    #     cw = float(self.tank_data['cw'])
    #     cf = float(self.tank_data['cf'])
    #     J = float(self.tank_data['J'])
    #     m = float(self.tank_data['initial_gascap'])
    #     winj = self.prod_table['wi']
    #     winj = winj.fillna(0)
    #     ginj = self.prod_table['gi']
    #     ginj = ginj.fillna(0)
    #     we = [None] * len(np)
    #     we[0] = 0
    #     #####General PVT
    #     # Tsc = 60  # F
    #     # Psc = 15.025  # psia
    #     # Tres = 219  # F
    #     # Pbp = df_pvtmaster['sat_press']  # psia
    #     Rsi = self.pvt_master['gor']  # scf/stb
    #
    #     Pi = float(self.tank_data['initial_pressure'])
    #     Boi = np.interp(Pi, self.prod_table['pressure'], self.oil_pvt_table['oil_fvf'])
    #     Bgi = np.interp(Pi, self.gas_pvt_table['pressure'], self.gas_pvt_table['gas_fvf']) / 1000
    #     Rsb = self.pvt_master['gor']
    #     Bti = mb.formation_total_volume_factor(Boi, Bgi, Rsb, Rsi)
    #     #####Water PVT
    #     Bw = 1.0  # dict_tank['Bw']
    #     Bwinj = 1.0
    #     #####Oil PVT
    #     pvt_oil_pressure = self.oil_pvt_table['pressure']
    #     pvt_oil_Bo = self.oil_pvt_table['oil_fvf']
    #     pvt_oil_Rs = self.oil_pvt_table['solution_gas']
    #     #####Gas PVT
    #     pvt_gas_pressure = self.gas_pvt_table['pressure']
    #     pvt_gas_Bg = self.gas_pvt_table['gas_fvf']
    #     pvt_gas_Bg = pvt_gas_Bg / 1000
    #     arr = np.array(pvt_oil_pressure)
    #     interpol = lambda P: np.interp(P, pvt_gas_pressure, pvt_gas_Bg)
    #     pvt_oil_Bg = interpol(arr)
    #     # pvt_oil_Bg = np.interp(P, pvt_gas_pressure, pvt_gas_Bg)
    #
    #     aquifer_pres = [None] * len(np)
    #     aquifer_pres[0] = Pi
    #     Pres_calc_empty = [None] * len(np)
    #     DDI = [None] * len(np)
    #     SDI = [None] * len(np)
    #     WDI = [None] * len(np)
    #     CDI = [None] * len(np)
    #
    #     for x in range(len(np)):
    #         if x == 0:
    #             DDI[x] = 0
    #             SDI[x] = 0
    #             WDI[x] = 0
    #             CDI[x] = 0
    #         else:
    #             P = pres_calc[x]
    #             Bo = np.interp(P, pvt_oil_pressure, pvt_oil_Bo)
    #
    #             Bg = np.interp(P, pvt_oil_pressure, pvt_oil_Bg)
    #             Bginj = Bg
    #             Rs = np.interp(P, pvt_oil_pressure, pvt_oil_Rs)
    #             Bt = mb.formation_total_volume_factor(Bo, Bg, Rsb, Rs)
    #             Eo = mb.dissolved_oil_and_gas_expansion(Bt, Bti)
    #             Eg2 = mb.gas_cap_expansion2(Bti, Bg, Bgi)
    #             Eg = mb.gas_cap_expansion(Bti, Bg, Bgi)
    #             dP = Pi - P
    #             Efw = mb.pore_volume_reduction_connate_water_expansion(m, Boi, cw, swi, cf, dP)
    #             F, produced_oil_and_gas, produced_water, injected_gas, injected_water = \
    #                 mb.production_injection_balance(np[x],Bt, Rs, Rsi, Bg, wp[x], Bw, winj[x], Bwinj, ginj[x], Bginj, gp[x])
    #             wex, aq_pres = itera.aquifer_influx(x, P, wei, we, ts, pres_calc, Pi, J, aquifer_pres)
    #             we[x] = wex
    #             aquifer_pres[x] = aq_pres
    #
    #             Ncalc = mb.oil_in_place(F, Eo, m, Eg, Efw, we[x], Bw, Bti)
    #             DDI[x] = Ncalc * Eo / F
    #             SDI[x] = Ncalc * m * Eg2 * (Boi / Bgi) / F
    #             WDI[x] = (we[x] * Bw - wp[x] * Bw) / F
    #             CDI[x] = Ncalc * (1 + m) * Efw * Boi / F
    #
    #     return DDI, SDI, WDI, CDI
    #


