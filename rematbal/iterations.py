import rematbal.matbal as mb
import numpy as np
from scipy.optimize import fsolve
import pandas as pd

def mbal_inner_calc(dict, P, Pres_calc, We, aquifer_pres, step, dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo,
                    pvt_oil_Bg, pvt_oil_Rs, Rsb,
                    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J,
                    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type):

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
    dates = data[5]
    Np = data[6]
    Wp = data[7]
    Gp = data[8]
    N = data[9]
    Wei = data[10]
    pvt_oil_pressure = data[11]
    pvt_oil_Bo = data[12]
    pvt_oil_Bg = data[13]
    pvt_oil_Rs = data[14]
    Rsb = data[15]
    Bti = data[16]
    Bgi = data[17]
    Pi = data[18]
    m = data[19]
    Boi = data[20]
    cw = data[21]
    Swi = data[22]
    cf = data[23]
    Rsi = data[24]
    Bw = data[25]
    Winj = data[26]
    Bwinj = data[27]
    Ginj = data[28]
    J = data[29]
    ts = data[30]
    VEH_aq_type = data[31]
    td_array = data[32]
    VEH_dp_array = data[33]
    r = data[34]
    rr = data[35]
    aq_type = data[36]


    #= mbal_setup(dict)
    F, Eo, m, Eg, Efw, We, aq_pres, Bw, Bti, N = mbal_inner_calc('', P, Pres_calc, We, aq_pres, step, dates,
                                                Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, \
                                                pvt_oil_Bg, pvt_oil_Rs, Rsb, \
                                                Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, \
                                                Bwinj, Ginj, J, ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
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
    #df_prod = dict['df_prod']
    dict_tank = dict['dict_tank']
    dict_pvtmaster = dict['dict_pvtmaster']
    df_pvt_gas = dict['df_pvt_gas']
    df_pvt_oil = dict['df_pvt_oil']
    dict_prod = dict['df_prod']
    dates = dict_prod['datestamp']
    dates = dates.to_numpy()
    #ts = pd.to_numeric(dates - dates.min()) / 864e11
    ts = dates - dates.min()
    ts = [time.days for time in ts] # ts.days
    Np = dict_prod['np']
    Gp = dict_prod['gp']
    Wp = dict_prod['wp']
    N = float(dict_tank['initial_inplace'])
    Swi = float(dict_tank['swi'])
    cw = float(dict_tank['cw'])
    cf = float(dict_tank['cf'])

    m = float(dict_tank['initial_gascap'])
    Winj = dict_prod['wi'].to_numpy()
    #Winj = Winj.fillna(0)
    Ginj = dict_prod['gi'].to_numpy()
    #Ginj = Ginj.fillna(0)
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
    return dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
           Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
           ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type


def eval_mbal_input(dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg,
                                    pvt_oil_Rs, Rsb, Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
                                    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type):
    Pres_calc = np.array([])
    #df_prod = dict['df_prod']
    #Np = df_prod['np'].to_numpy()
    We = np.array([None] * len(Np))
    aquifer_pres = np.array([None] * len(Np))
    #dict_tank = dict['dict_tank']
    #Pi = float(dict_tank['initial_pressure'])
    #Np = df_prod['np'].to_numpy()
    #Wp = df_prod['wp'].to_numpy()
    #Gp = df_prod['wp'].to_numpy()
    #dates = df_prod['datestamp'].to_numpy()
    #Winj = df_prod['wi'].to_numpy()
    #Winj = Winj.fillna(0)
    #Ginj = df_prod['gi'].to_numpy()
    #Ginj = Ginj.fillna(0)
    dict_prod = {
        'np': Np,
        'wp': Wp,
        'gp': Gp,
        'wi': Winj,
        'gi': Ginj,
        'datestamp': dates
    }
    #dict['dict_prod'] = dict_prod

    for x in range(len(Np)):
        if x == 0:
            aquifer_pres[x] = Pi
            We[x] = 0.0
            Pres_calc = np.append(Pres_calc, [Pi])
            #Pres_calc.append(Pi)
        else:
            data = ('', Pres_calc, We, aquifer_pres, x, dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg,
                                    pvt_oil_Rs, Rsb, Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
                                    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
            Pres_calc = np.append(Pres_calc, [pressure_calculation(data, Pres_calc)[0]])
            #Pres_calc.append(pressure_calculation(data, Pres_calc)[0])

    #dict['Pres_calc'] = Pres_calc
    solution_set = drive_indices(Pres_calc, dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
           Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
           ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
    return solution_set


def drive_indices(Pres_calc, dates, Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
    Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
    ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type):
    solution_set = pd.DataFrame()
    DDI = []
    SDI = []
    WDI = []
    CDI = []
    #Pres_calc = dict['Pres_calc']
    Ncalc_array = []
    #df_prod = dict['df_prod']
    #Np = df_prod['np']
    We = [None] * len(Np)
    aquifer_pres = [None] * len(Np)
    #Wp = df_prod['wp']
    #dict_tank = dict['dict_tank']
    #Pi = float(dict_tank['initial_pressure'])
    #N = float(dict_tank['initial_inplace'])
    #Boi = dict_tank['Boi']
    #Bgi = dict_tank['Bgi']
    #dates = df_prod['datestamp']
    #ts = pd.to_numeric(dates - dates.min()) / 864e11
    ts = dates - dates.min()
    ts = [time.days for time in ts]
    #Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, pvt_oil_Bg, pvt_oil_Rs, Rsb, \
    #Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, Bwinj, Ginj, J, \
    #ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type = mbal_setup(dict)

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
            F, Eo, m, Eg, Efw, We, aquifer_pres, Bw, Bti, N = mbal_inner_calc('', Pres_calc[x], Pres_calc, We, aquifer_pres, x, dates,
                                                Np, Wp, Gp, N, Wei, pvt_oil_pressure, pvt_oil_Bo, \
                                                pvt_oil_Bg, pvt_oil_Rs, Rsb, \
                                                Bti, Bgi, Pi, m, Boi, cw, Swi, cf, Rsi, Bw, Winj, \
                                                Bwinj, Ginj, J, ts, VEH_aq_type, td_array, VEH_dp_array, r, rr, aq_type)
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
