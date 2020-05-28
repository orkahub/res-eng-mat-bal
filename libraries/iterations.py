import libraries.matbal as mb
import numpy as np
from scipy.optimize import fsolve
import pandas as pd
from scipy import interpolate


def mbal_inner_calc(dict, P, Pres_calc, We, aquifer_pres, step):
    Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type = mbal_setup(dict)
    Bo = np.interp(P, pvt_oil_pressure, pvt_oil_Bo)
    # print(Bo, P)
    Bg = np.interp(P, pvt_oil_pressure, pvt_oil_Bg)
    Bginj = Bg
    Rs = np.interp(P, pvt_oil_pressure, pvt_oil_Rs)
    Bt = mb.formation_total_volume_factor(Bo, Bg, Rsb, Rs)
    # print("Bt = ", Bt)
    Eo = mb.dissolved_oil_and_gas_expansion(Bt, Bti)
    Eg = mb.gas_cap_expansion(Bti, Bg, Bgi)
    dP = Pi - P
    #Wex = 0.0
    Efw = mb.pore_volume_reduction_connate_water_expansion(m, Boi, cw, Swi, cf, dP)
    Npx = Np[step]
    Wpx = Wp[step]
    Gpx = Gp[step]
    Winjx = Winj[step]
    Ginjx = Ginj[step]
    #We[0] = 0
    F, produced_oil_and_gas, produced_water, injected_gas, injected_water = mb.production_injection_balance(Npx, Bt, Rs,
                                                                                                            Rsi, Bg,
                                                                                                            Wpx,
                                                                                                            Bw, Winjx,
                                                                                                            Bwinj, Ginjx,
                                                                                                            Bginj, Gpx)
    if aq_type == 'VEH':
        #VEH_avg_pressure = data[39]
        Wex, VEH_avg_pressure, VEH_dp_array = mb.VEH_aquifer_influx(VEH_aq_type, step, ts, td_array, VEH_dp_array,
                                                                    r, rr, P, Pi, VEH_avg_pressure)
    if aq_type == 'Fetkovich':
        Wex, aq_pres = aquifer_influx(step, P, Wei, We, ts, Pres_calc, Pi, J, aquifer_pres)
        #aquifer_pres.append(aq_pres)
    We.append(Wex)
    return F, Eo, m, Eg, Efw, We, aq_pres, Bw, Bti, N


def obj_funtion2(P, *data):
    dict = data[0]
    Pres_calc = data[1]
    We = data[2]
    aq_pres = data[3]
    step = data[4]
    # P = P[0]
    # Npx = data[0]
    # Wpx = data[1]
    # # Wex = data[2]
    # Gpx = data[2]
    # N = data[3]
    # step = data[4]
    # Wei = data[5]
    # pvt_oil_pressure = data[6]
    # pvt_oil_Bo = data[7]
    # pvt_oil_Bg = data[8]
    # pvt_oil_Rs = data[9]
    # Rsb = data[10]
    # Bti = data[11]
    # Bgi = data[12]
    # Pi = data[13]
    # m = data[14]
    # Boi = data[15]
    # cw = data[16]
    # Swi = data[17]
    # cf = data[18]
    # Rsi = data[19]
    # Bw = data[20]
    # Winj = data[21]
    # Bwinj = data[22]
    # Ginj = data[23]
    # We = data[24]
    # Pres_calc = data[25]
    # aquifer_pres = data[26]
    # J = data[27]
    # ts = data[28]
    # # DDI = data[29]
    # # SDI = data[30]
    # # WDI = data[31]
    # # CDI = data[32]
    # VEH_aq_type = data[29]
    # td_array = data[30]
    # VEH_dp_array = data[31]
    # r = data[32]
    # rr =  data[33]
    # aq_type = data[34]
    F, Eo, m, Eg, Efw, We, aq_pres, Bw, Bti, N = mbal_inner_calc(dict, P, Pres_calc, We, aq_pres, step)
    Wex = We[step]


    #mbal_res = Mbal_resutls(Pres_calc, We, aq_pres)
    Ncalc = mb.oil_in_place(F, Eo, m, Eg, Efw, Wex, Bw, Bti)

    of = (N - Ncalc)
    return of


def pressure_calculation(data, Pres_calc):
    step = data[4]
    x0 = Pres_calc[step - 1] - 10.0
    res = fsolve(obj_funtion2, x0, args=data)
    return res


def aquifer_influx(step, P, Wei, We, ts, Pres_calc, Pi, J, aquifer_pres):
    We_prev = We[step - 1]
    ts_prev = ts[step - 1]
    tsx = ts[step]
    avg_pres = (Pres_calc[step - 1] + P) / 2
    aq_pres = aquifer_pressure(step, Wei, We, aquifer_pres, Pi)
    # print(step,aq_pres)
    Wex = (Wei / Pi) * (aq_pres - avg_pres) * (1 - np.exp(-J * Pi * (tsx - ts_prev) / Wei))
    Wex = We_prev + Wex
    return Wex, aq_pres


def aquifer_pressure(step, Wei, We, aquifer_pres, Pi):
    We_prev = We[step - 1]
    if step == 1:
        aq_pres = Pi
    else:
        aq_pres = Pi * (1 - We_prev / (Wei))
    #aquifer_pres.append(aq_pres)
    return aq_pres


def mbal_setup(dict):

    #dict_temp = dict
    #dict_temp['df_prod'] = None
    df_prod = dict['df_prod']
    dict_tank = dict['dict_tank']
    dict_pvtmaster = dict['dict_pvtmaster']
    df_pvt_gas = dict['df_pvt_gas']
    df_pvt_oil = dict['df_pvt_oil']

    dates = df_prod['datestamp']
    ts = pd.to_numeric(dates - dates.min()) / 864e11
    Np = df_prod['np']
    Gp = df_prod['gp']
    # Gp = Gp * 1000.0
    Wp = df_prod['wp']
    reservoir_pressure_obs = df_prod['pressure']
    ts_obs = ts[reservoir_pressure_obs.notnull()]
    We = [None] * len(Np)
    Pres_calc = []
    reservoir_pressure_obs = reservoir_pressure_obs[reservoir_pressure_obs.notnull()]
    reservoir_pressure_obs = reservoir_pressure_obs * 1.0
    N = float(dict_tank['initial_inplace'])

    Swi = float(dict_tank['swi'])
    cw = float(dict_tank['cw'])
    cf = float(dict_tank['cf'])
    # N = 1.00E+07
    # Wei = 5.00E+07

    m = float(dict_tank['initial_gascap'])
    # We = 0
    Winj = df_prod['wi']
    Winj = Winj.fillna(0)
    Ginj = df_prod['gi']
    Ginj = Ginj.fillna(0)
    #####General PVT
    # Tsc = 60  # F
    # Psc = 15.025  # psia
    # Tres = 219  # F
    # Pbp = df_pvtmaster['sat_press']  # psia
    Rsi = dict_pvtmaster['gor']  # scf/stb
    try:
        aq_type = dict_tank['aq_type']
    except:
        aq_type = "Fetkovich"
    VEH_aq_type = ""
    r = ""
    rr = ""
    td_array = ""
    VEH_dp_array = ""

    if aq_type == 'Fetkovich':
        Wei = float(dict_tank['wei'])
        J = float(dict_tank['J'])
    if aq_type == 'VEH':
        VEH_dp_array = [None] * len(Np)
        VEH_aq_type = dict_tank['VEH_aq_type']
        r = dict_tank['r']
        rr = dict_tank['rr']
        td_array = mb.VEH_td(VEH_aq_type, dict_tank['k'], ts, dict_tank['poro'], dict_tank['visc'], dict_tank['ct'], rr,
                          dict_tank['La'])
    else:
        Wei = float(dict_tank['wei'])
        J = float(dict_tank['J'])



    Pi = float(dict_tank['initial_pressure'])
    #f_interpol = interpolate.interp1d(df_pvt_oil['pressure'], df_pvt_oil['oil_fvf'],
    #                                  bounds_error=False, fill_value="extrapolate")
    #Boi = np.interp(Pi, df_pvt_oil['pressure'], df_pvt_oil['oil_fvf'])
    #Bgi = np.interp(Pi, df_pvt_gas['pressure'], df_pvt_gas['gas_fvf']) / 1000
    #Boi = f_interpol(Pi)
    #f_interpol = interpolate.interp1d(df_pvt_oil['pressure'], df_pvt_gas['gas_fvf'],
    #                                  bounds_error=False, fill_value="extrapolate")
    #Bgi = f_interpol(Pi)
    Boi = dict_tank['Boi']
    Bgi = dict_tank['Bgi']
    Rsb = dict_pvtmaster['gor']
    Bti = mb.formation_total_volume_factor(Boi, Bgi, Rsb, Rsi)
    #####Water PVT
    Bw = 1.0  # dict_tank['Bw']
    Bwinj = 1.0
    #####Oil PVT
    pvt_oil_pressure = df_pvt_oil['pressure']
    pvt_oil_Bo = df_pvt_oil['oil_fvf']
    pvt_oil_Rs = df_pvt_oil['solution_gas']
    #####Gas PVT
    pvt_gas_pressure = df_pvt_gas['pressure']
    pvt_gas_Bg = df_pvt_gas['gas_fvf']
    pvt_gas_Bg = pvt_gas_Bg / 1000
    arr = np.array(pvt_oil_pressure)
    interpol = lambda P: np.interp(P, pvt_gas_pressure, pvt_gas_Bg)
    pvt_oil_Bg = interpol(arr)
    # pvt_oil_Bg = np.interp(P, pvt_gas_pressure, pvt_gas_Bg)

    aquifer_pres = [None] * len(Np)
    # Pres_calc_empty = [None] * len(Np)
    # DDI = [None] * len(Np)
    # SDI = [None] * len(Np)
    # WDI = [None] * len(Np)
    # CDI = [None] * len(Np)

    return Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
           Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
           ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type


def eval_mbal_input(dict):
    Pres_calc = []
    We = []
    aquifer_pres = []
    df_prod = dict['df_prod']
    Np = df_prod['np']
    dict_tank = dict['dict_tank']
    Pi = float(dict_tank['initial_pressure'])
    #Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
    #Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
    #ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type = mbal_setup(dict)


    for x in range(len(Np)):
        if x == 0:
            #Pres_calc.append(Pi)
            We.append(0.0)
            #aquifer_pres.append(Pi)
            Pres_calc.append(Pi)
        else:
            # data = (Np[x], Wp[x], Gp[x], N, x, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb,
            #         Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj[x], Bwinj, Ginj[x], We, Pres_calc, aquifer_pres, J,
            #         ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
            data = (dict, Pres_calc, We, aquifer_pres, x)
            Pres_calc.append(pressure_calculation(data, Pres_calc)[0])


            # prod_dict = {
            #     'np': Np[x],
            #     'wp': Wp[x],
            #     'gp': Gp[x],
            #     'winj': Winj[x],
            #     'ginj': Ginj[x],
            #     'pres_calc': Pres_calc[x]
            # }
            #dict_temp['prod_dict'] = prod_dict
    dict['Pres_calc'] = Pres_calc
    solution_set = drive_indices(dict)
    return solution_set


def drive_indices(dict):
    solution_set = pd.DataFrame()
    DDI = []
    SDI = []
    WDI = []
    CDI = []
    Pres_calc = dict['Pres_calc']
    Ncalc_array = []
    # Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
    # Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
    # ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type = mbal_setup(dict)
    aquifer_pres = []
    df_prod = dict['df_prod']
    Np = df_prod['np']
    Wp = df_prod['wp']
    dict_tank = dict['dict_tank']
    Pi = float(dict_tank['initial_pressure'])
    N = float(dict_tank['initial_inplace'])
    Boi = dict_tank['Boi']
    Bgi = dict_tank['Bgi']
    #aquifer_pres = [None] * len(Np)
    #aquifer_pres.append(Pi)
    dates = df_prod['datestamp']
    ts = pd.to_numeric(dates - dates.min()) / 864e11
    # Pres_calc_empty = [None] * len(Np)
    We = []
    We.append(0)




    for x in range(len(Np)):
        if x == 0:
            DDI.append(0)
            SDI.append(0)
            WDI.append(0)
            CDI.append(0)
            Ncalc_array.append(N)
        else:
            F, Eo, m, Eg, Efw, We, aquifer_pres, Bw, Bti, N = mbal_inner_calc(dict, Pres_calc[x], Pres_calc, We, aquifer_pres, x)
            # Wex = We[x]
            # P = Pres_calc[x]
            # Bo = np.interp(P, pvt_oil_pressure, pvt_oil_Bo)
            #
            # Bg = np.interp(P, pvt_oil_pressure, pvt_oil_Bg)
            # Bginj = Bg
            # Rs = np.interp(P, pvt_oil_pressure, pvt_oil_Rs)
            # Bt = mb.formation_total_volume_factor(Bo, Bg, Rsb, Rs)
            # Eo = mb.dissolved_oil_and_gas_expansion(Bt, Bti)
            # Eg2 = mb.gas_cap_expansion2(Bti, Bg, Bgi)
            # Eg = mb.gas_cap_expansion(Bti, Bg, Bgi)
            # dP = Pi - P
            # Efw = mb.pore_volume_reduction_connate_water_expansion(m, Boi, cw, Swi, cf, dP)
            # F, produced_oil_and_gas, produced_water, injected_gas, injected_water = \
            #     mb.production_injection_balance(Np[x],Bt, Rs, Rsi, Bg, Wp[x], Bw, Winj[x], Bwinj, Ginj[x], Bginj, Gp[x])
            # Wex, aq_pres = aquifer_influx(x, P, Wei, We, ts, Pres_calc, Pi, J, aquifer_pres)
            # We[x] = Wex
            # aquifer_pres[x] = aq_pres

            Ncalc = mb.oil_in_place(F, Eo, m, Eg, Efw, We[x], Bw, Bti)
            DDI.append(Ncalc * Eo / F)
            SDI.append(Ncalc * m * Eg * (Boi / Bgi) / F)
            WDI.append((We[x] * Bw - Wp[x] * Bw) / F)
            CDI.append(Ncalc * (1 + m) * Efw * Boi / F)
            Ncalc_array.append(Ncalc)
            #aquifer_pres.append(aq_pres)
            # of = (N - Ncalc)

    solution_set['ts'] = ts
    solution_set['pres_calc'] = Pres_calc
    solution_set['ddi'] = DDI
    solution_set['sdi'] = SDI
    solution_set['wdi'] = WDI
    solution_set['cdi'] = CDI
    solution_set['we'] = We
    solution_set['aquifer_pres'] = aquifer_pres
    solution_set['oip'] = Ncalc_array
    return solution_set
