import math

import numpy as np


def VEH_aquifer_influx(VEH_aq_type, step, ts, td_array, VEH_dp_array, r, rr, P, Pi, VEH_avg_pressure):
    return Wex, VEH_avg_pressure, VEH_dp_array

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
    tdstar = 0.4*(red*red-1)
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


def fet_aquifer_influx(step, P, Wei, We, ts, Pres_calc, Pi, J, aquifer_pres):
    We_prev = We[step - 1]
    ts_prev = ts[step - 1]
    tsx = ts[step]
    avg_pres = (Pres_calc[step - 1] + P) / 2
    aq_pres = fet_aquifer_pressure(step, Wei, We, aquifer_pres, Pi)
    # print(step,aq_pres)
    Wex = (Wei / Pi) * (aq_pres - avg_pres) * (1 - np.exp(-J * Pi * (tsx - ts_prev) / Wei))
    Wex = We_prev + Wex
    return Wex, aq_pres


def fet_aquifer_pressure(step, Wei, We, aquifer_pres, Pi):
    We_prev = We[step - 1]
    if step == 1:
        aq_pres = Pi
    else:
        aq_pres = Pi * (1 - We_prev / (Wei))
    return aq_pres
