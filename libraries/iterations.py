import libraries.matbal as mb
import numpy as np
from scipy.optimize import fsolve
import pandas as pd


def obj_funtion2(P, *data):
    P = P[0]
    Npx = data[0]
    Wpx = data[1]
    # Wex = data[2]
    Gpx = data[2]
    N = data[3]
    step = data[4]
    Wei = data[5]
    pvt_oil_pressure = data[6]
    pvt_oil_Bo = data[7]
    pvt_oil_Bg = data[8]
    pvt_oil_Rs = data[9]
    Rsb = data[10]
    Bti = data[11]
    Bgi = data[12]
    Pi = data[13]
    m = data[14]
    Boi = data[15]
    cw = data[16]
    Swi = data[17]
    cf = data[18]
    Rsi = data[19]
    Bw = data[20]
    Winj = data[21]
    Bwinj = data[22]
    Ginj = data[23]
    We = data[24]
    Pres_calc = data[25]
    aquifer_pres = data[26]
    J = data[27]
    ts = data[28]
    DDI = data[29]
    SDI = data[30]
    WDI = data[31]
    CDI = data[32]
    VEH_aq_type = data[33]
    td_array = data[34]
    VEH_dp_array = data[35]
    r = data[36]
    rr =  data[37]
    aq_type = data[38]
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
    Wex = 0.0
    Efw = mb.pore_volume_reduction_connate_water_expansion(m, Boi, cw, Swi, cf, dP)
    F, produced_oil_and_gas, produced_water, injected_gas, injected_water = mb.production_injection_balance(Npx, Bt, Rs,
                                                                                                         Rsi, Bg, Wpx,
                                                                                                         Bw, Winj,
                                                                                                         Bwinj, Ginj,
                                                                                                         Bginj, Gpx)
    if aq_type == 'VEH':
        VEH_avg_pressure = data[39]
        Wex, VEH_avg_pressure, VEH_dp_array  = mb.VEH_aquifer_influx(VEH_aq_type, step, ts, td_array, VEH_dp_array,
                                                                  r, rr, P, Pi, VEH_avg_pressure)
    if aq_type == 'Fetkovich':
        Wex, aq_pres = aquifer_influx(step, P, Wei, We, ts, Pres_calc, Pi, J, aquifer_pres)
        aquifer_pres[step] = aq_pres
    We[step] = Wex

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
    aquifer_pres[step] = aq_pres
    return aq_pres


def eval_mbal_input2(dict):
    df_prod = dict['df_prod']
    dict_tank = dict['dict_tank']
    dict_pvtmaster = dict['dict_pvtmaster']
    df_pvt_gas = dict['df_pvt_gas']
    df_pvt_oil = dict['df_pvt_oil']

    dates = df_prod['datestamp']
    ts = pd.to_numeric(dates - df_prod['datestamp'].min()) / 864e11
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
    Boi = np.interp(Pi, df_pvt_oil['pressure'], df_pvt_oil['oil_fvf'])
    Bgi = np.interp(Pi, df_pvt_gas['pressure'], df_pvt_gas['gas_fvf']) / 1000
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
    Pres_calc_empty = [None] * len(Np)
    DDI = [None] * len(Np)
    SDI = [None] * len(Np)
    WDI = [None] * len(Np)
    CDI = [None] * len(Np)
    Pres_calc = [None] * len(Np)

    for x in range(len(Np)):
        if x == 0:
            #Pres_calc.append(Pi)
            We[x] = 0.0
            aquifer_pres[0] = Pi
            Pres_calc[0] = Pi
        else:
            data = (Np[x], Wp[x], Gp[x], N, x, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb,
                    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj[x], Bwinj, Ginj[x], We, Pres_calc, aquifer_pres, J,
                    ts, DDI, SDI, WDI, CDI, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
            Pres_calc[x] = pressure_calculation(data, Pres_calc)[0]
    return Pres_calc, ts_obs, reservoir_pressure_obs, ts
