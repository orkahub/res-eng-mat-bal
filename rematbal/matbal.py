"""
" Some of these functions were taken from
https://github.com/ESSS/kraken-macros/blob/master/src/macros/mbal/mbalcore/mbal_functions.py
"""
import numpy as np


# class Mbal_resutls:
#     def __init__(self, Pres_calc, We, aq_press):
#         self.Pres_calc = Pres_calc
#         self.We = We
#         self.aq_press = aq_press


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

