"""
" Some of these functions were taken from
https://github.com/ESSS/kraken-macros/blob/master/src/macros/mbal/mbalcore/mbal_functions.py
"""
import json
import numpy as np
import pandas as pd
import plotly
from scipy.optimize import curve_fit
import math


class Mbal_resutls:
    def __init__(self, Pres_calc, We, aq_press):
        self.Pres_calc = Pres_calc
        self.We = We
        self.aq_press = aq_press


def VEH_td(VEH_aq_type, k, t, poro, visc, ct, rr, La):
    if VEH_aq_type == 'radial':
        td = 0.00633*k*t/(poro*visc*ct*rr*rr)
    if VEH_aq_type == 'linear':
        td = (2.309 * k * t) / (poro * visc * ct * La * La)
    return td # t in days


def VEH_Pd_constrate(k, h, pi, p, q, visc):
    return 0.00708*k*h*(pi-p)/(q*visc)


def VEH_Pd_constantpress(pi, p, pr):
    return (pi - p)/(pi - pr)


def VEH_rd(r, rr):
    return r/rr


def VEH_inifinite_aq_influx_const(poro, cr, h, rr, dp, Wd, angle):
    return 1.119*poro*cr*h*rr*rr*(angle/360)


def VEH_dimensionless_radialaq(td, t, red):
    tdstar = 0.4(red*red-1)
    Jstar = (math.pow(red,4)*math.log(red))/(red*red-1)+0.25*(1-3*red*red)
    if td > tdstar:
        Wd = 0.5*(red*red-1)*(1-math.exp(-2*td/Jstar))
    if td <= tdstar:
        if t <= 100:
            a7 = 4.8534e-12
            a6 = -1.8436e-9
            a5 = 2.8354e-7
            a4 = -2.2740e-5
            a3 = 1.0284e-3
            a2 = -2.7455e-2
            a1 = 8.5373e-1
            a0 = 8.1638e-1
            Wd = a7 * math.pow(td, 7) + a6 * math.pow(td, 6) + a5 * math.pow(td, 5) + a4 * math.pow(td, 4) \
                 + a3 * math.pow(td, 3) + a2 * math.pow(td, 2) + a1 * math.pow(td, 1) + a0
        if t > 100:
            Wd = 2 * td / math.log(td)
    return Wd


def VEH_linaq_La(VEH_linaq_type, Lr, h, Vpa, Vpr, poror, poroa):
    if VEH_linaq_type == 'edge':
        La = Lr*(Vpa/Vpr)(poror/poroa)
    if VEH_linaq_type == 'bottom':
        La = h*(Vpa/Vpr)(poror/poroa)
    return La


def VEH_dimensionless_linearaq(td):
    if td < 0.47:
        Wd = 2*math.sqrt(td/math.pi)
    if td >= 0.47:
        Wd= 1 - 0.065807*math.pow(td,-16512)
    return Wd


def VEH_dimensionless_influx(VEH_aq_type, dtd, t, r, rr):
    if VEH_aq_type == 'radial':
        red = VEH_rd(r, rr)
        VEH_dimensionless_radialaq(dtd, t, red)
    if VEH_aq_type == 'linear':
        VEH_dimensionless_linearaq(dtd)


def formation_total_volume_factor(Bo, Bg, Rsb, Rs):
    return np.where(Rs >= Rsb, Bo, Bo + Bg * (Rsb - Rs))


def formation_volume_factor_gas_from_z(z, Tres, Psc, P, Tsc):
    Bg = z * (1000 / 5.615) * (Tres + 460) * Psc / (P * (Tsc + 460))
    return Bg


def production_injection_balance(Np, Bt, Rs, Rsi, Bg, Wp, Bw, Winj, Bwinj, Ginj, Bginj, Gp):
    produced_oil_and_gas = (Np * (Bt + (Gp / Np - Rsi) * Bg))
    produced_water = Wp * Bw
    injected_water = Winj * Bwinj
    injected_gas = Ginj * Bginj

    F = (produced_oil_and_gas + produced_water - injected_water - injected_gas)

    return F, produced_oil_and_gas, produced_water, injected_gas, injected_water


def dissolved_oil_and_gas_expansion(Bt, Bti):
    Eo = (Bt - Bti)

    return Eo


def dissolved_oil_and_gas_expansion2(Bo, Boi, Rsi, Bg, Rs):
    Eo = (Bo - Boi + Bg * (Rsi - Rs))

    return Eo


def gas_cap_expansion(Bti, Bg, Bgi):
    Eg = ((Bg / Bgi) - 1)

    return Eg


def gas_cap_expansion2(Bti, Bg, Bgi):
    Eg2 = (Bg - Bgi)

    return Eg2


def deltaP(Pi, Pavg):
    deltaP = Pi - Pavg

    return deltaP


def pore_volume_reduction_connate_water_expansion(m, Boi, cw, Swi, cf, deltaP):
    Efw = ((cw * Swi + cf) / (1.0 - Swi)) * deltaP

    return Efw


def oil_in_place(F, Eo, m, Eg, Efw, We, Bw, Bti):
    oil_in_place = (F - We * Bw) / (Eo + (Bti) * m * Eg + (1.0 + m) * Bti * Efw)

    return oil_in_place


def oil_in_place_underg_withdrawal(F, Eo):
    oil_in_place_underg_withdrawal = F / Eo

    return oil_in_place_underg_withdrawal


def oil_in_place_gas_cap(F, Eo, m, Eg, Bti):
    oil_in_place_gas_cap = F / (Eo + (Bti) * m * Eg)


def oil_in_place_water_influx(F, We, Eo):
    oil_in_place_water_influx = (F - We) / Eo

    return oil_in_place_water_influx