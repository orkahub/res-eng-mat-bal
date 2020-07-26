import numpy as np
import rematbal.matbal as mb
import pandas as pd
import rematbal.iterations as itera


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
            Pres_calc[x] = itera.pressure_calculation(data, Pres_calc)[0]
    return Pres_calc, ts_obs, reservoir_pressure_obs, ts
