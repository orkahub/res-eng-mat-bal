import rematbal.matbal as mb
import numpy as np
from scipy.optimize import fsolve
import pandas as pd

def mbal_inner_calc(dict, P, Pres_calc, We, aquifer_pres, step):
    Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type = mbal_setup(dict)
    Bo = np.interp(P, pvt_oil_pressure, pvt_oil_Bo)
    Bg = np.interp(P, pvt_oil_pressure, pvt_oil_Bg)
    Bginj = Bg
    Rs = np.interp(P, pvt_oil_pressure, pvt_oil_Rs)
    Bt = mb.formation_total_volume_factor(Bo, Bg, Rsb, Rs)
    Eo = mb.dissolved_oil_and_gas_expansion(Bt, Bti)
    Eg = mb.gas_cap_expansion(Bti, Bg, Bgi)
    dP = Pi - P
    Efw = mb.pore_volume_reduction_connate_water_expansion(m, Boi, cw, Swi, cf, dP)
    Npx = Np[step]
    Wpx = Wp[step]
    Gpx = Gp[step]
    Winjx = Winj[step]
    Ginjx = Ginj[step]
    F, produced_oil_and_gas, produced_water, injected_gas, injected_water = mb.production_injection_balance(Npx, Bt, Rs,
                                                                                                            Rsi, Bg,
                                                                                                            Wpx,
                                                                                                            Bw, Winjx,
                                                                                                            Bwinj, Ginjx,
                                                                                                            Bginj, Gpx)
    if aq_type == 'VEH':
        Wex, VEH_avg_pressure, VEH_dp_array = mb.VEH_aquifer_influx(VEH_aq_type, step, ts, td_array, VEH_dp_array,
                                                                    r, rr, P, Pi, VEH_avg_pressure)
    if aq_type == 'Fetkovich':
        Wex, aq_pres = aquifer_influx(step, P, Wei, We, ts, Pres_calc, Pi, J, aquifer_pres)
        aquifer_pres[step] = aq_pres
    We[step] = Wex
    return F, Eo, m, Eg, Efw, We, aquifer_pres, Bw, Bti, N


def obj_funtion2(P, *data):
    dict = data[0]
    Pres_calc = data[1]
    We = data[2]
    aq_pres = data[3]
    step = data[4]
    F, Eo, m, Eg, Efw, We, aq_pres, Bw, Bti, N = mbal_inner_calc(dict, P, Pres_calc, We, aq_pres, step)
    Wex = We[step]
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
    return aq_pres


def mbal_setup(dict):
    df_prod = dict['df_prod']
    dict_tank = dict['dict_tank']
    dict_pvtmaster = dict['dict_pvtmaster']
    df_pvt_gas = dict['df_pvt_gas']
    df_pvt_oil = dict['df_pvt_oil']

    dates = df_prod['datestamp']
    ts = pd.to_numeric(dates - dates.min()) / 864e11
    Np = df_prod['np']
    Gp = df_prod['gp']
    Wp = df_prod['wp']
    N = float(dict_tank['initial_inplace'])
    Swi = float(dict_tank['swi'])
    cw = float(dict_tank['cw'])
    cf = float(dict_tank['cf'])

    m = float(dict_tank['initial_gascap'])
    Winj = df_prod['wi']
    Winj = Winj.fillna(0)
    Ginj = df_prod['gi']
    Ginj = Ginj.fillna(0)
    #####General PVT
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
    aquifer_pres = [None] * len(Np)
    return Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
           Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
           ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type


def eval_mbal_input(dict):
    Pres_calc = []
    df_prod = dict['df_prod']
    Np = df_prod['np']
    We = [None] * len(Np)
    aquifer_pres = [None] * len(Np)
    dict_tank = dict['dict_tank']
    Pi = float(dict_tank['initial_pressure'])

    for x in range(len(Np)):
        if x == 0:
            aquifer_pres[x] = Pi
            We[x] = 0.0
            Pres_calc.append(Pi)
        else:
            data = (dict, Pres_calc, We, aquifer_pres, x)
            Pres_calc.append(pressure_calculation(data, Pres_calc)[0])

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
    df_prod = dict['df_prod']
    Np = df_prod['np']
    We = [None] * len(Np)
    aquifer_pres = [None] * len(Np)
    Wp = df_prod['wp']
    dict_tank = dict['dict_tank']
    Pi = float(dict_tank['initial_pressure'])
    N = float(dict_tank['initial_inplace'])
    Boi = dict_tank['Boi']
    Bgi = dict_tank['Bgi']
    dates = df_prod['datestamp']
    ts = pd.to_numeric(dates - dates.min()) / 864e11

    for x in range(len(Np)):
        if x == 0:
            We[x] = 0.0
            aquifer_pres[x] = Pi
            DDI.append(0)
            SDI.append(0)
            WDI.append(0)
            CDI.append(0)
            Ncalc_array.append(N)
        else:
            F, Eo, m, Eg, Efw, We, aquifer_pres, Bw, Bti, N = mbal_inner_calc(dict, Pres_calc[x], Pres_calc, We, aquifer_pres, x)
            Ncalc = mb.oil_in_place(F, Eo, m, Eg, Efw, We[x], Bw, Bti)
            DDI.append(Ncalc * Eo / F)
            SDI.append(Ncalc * m * Eg * (Boi / Bgi) / F)
            WDI.append((We[x] * Bw - Wp[x] * Bw) / F)
            CDI.append(Ncalc * (1 + m) * Efw * Boi / F)
            Ncalc_array.append(Ncalc)

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
