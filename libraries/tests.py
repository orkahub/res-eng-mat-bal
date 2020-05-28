from numpy.testing import assert_array_almost_equal
import numpy as np
import os
import pandas as pd
import datetime
import libraries.main as main


def test_calculated_pressure_at_different_timesteps():
    '''
    Test if Material Balance matched initial LSU Excel and its calculated pressure.
    :return:
    '''
    # Setup

    regress = False
    regress_config = None
    path = os.path.dirname(__name__)
    file = os.path.join(path, '../data/lsu_matbal/DONOTDELETE_TEST_OilMBEx.xls')
    keep_cols = ['Days', 'Np\nSTBO', 'Gp\nMCF', 'Wp\nSTBW', 'Meas P\npsia']
    df_prod = pd.read_excel(file, sheet_name='Calculations', skiprows=12, nrows=73)
    df_prod.rename(columns={'Np\nSTBO': 'np', 'Gp\nMCF': 'gp', 'Wp\nSTBW': 'wp', 'Meas P\npsia': 'pressure'},
                   inplace=True)
    startdate = datetime.date(2010, 1, 1)
    df_prod['datestamp'] = [datetime.timedelta(day) + startdate for day in df_prod['Days']]
    df_prod['wi'] = 0
    df_prod['gi'] = 0
    df_prod['gp'] = df_prod['gp'] * 1000.0
    keep_cols = ['Pressure\npsia', 'Bo\nrb/stbo', 'Rs\nscf/stbo', 'Bg rb/mcf', 'Bt rb/stbo']
    df_pvt_o = pd.read_excel(file, sheet_name='Oil PVT', skiprows=3, nrows=16)
    df_pvt_o.rename(columns={'Pressure\npsia': 'pressure', 'Bo\nrb/stbo': 'oil_fvf', 'Rs\nscf/stbo': 'solution_gas'},
                    inplace=True)
    keep_cols = ['Pressure psia', 'z', 'Bg\nrb/mcf']
    df_pvt_g = pd.read_excel(file, sheet_name='Z Factors', skiprows=3, nrows=12)
    df_pvt_g.rename(columns={'Pressure\npsia': 'pressure', 'Bg\nrb/mcf': 'gas_fvf'}, inplace=True)
    tank = {
        'initial_inplace': 13.7598436303087E6,
        'initial_gascap': 0,
        'initial_pressure': 10180,
        'wei': 141.815145478295e6,
        'J': 0.931250892112226,
        'swi': 0.2,
        'cw': 2.5e-6,
        'cf': 3e-5,
        'Boi': 1.7349,
        'Bgi': 0.6508
    }
    pvt_master = {
        'gor': 1720,
        'sat_press': 8227,
        'temperature': 219,
    }
    #df_prod['days'] = np.array([181, 577, 1338, 2160])
    #df_prod['measured_pressures'] = np.array([9643, 8778, 7202, 6877])

    # Exercise

    ts, Pres_calc, ts_obs, reservoir_pressure_obs, DDI, SDI, WDI, CDI, N, Wei, J = main.matbal_run(tank, df_prod,
                                                                                                   pvt_master, df_pvt_o,
                                                                                                   df_pvt_g, regress,
                                                                                                   regress_config)
    sum_indeces = np.array(DDI) + np.array(SDI) + np.array(WDI) + np.array(CDI)
    ones_array = np.repeat(1.0, len(DDI))

    # Verify

    #assert_array_almost_equal(df_prod['P\npsia'][1:10], Pres_calc[1:10], 0)
    assert_array_almost_equal(sum_indeces[50:70], ones_array[50:70], 3)
    assert_array_almost_equal(df_prod['P\npsia'][30:40], Pres_calc[30:40], 0)

    # Cleanup - none


def test_if_plots_are_ploted():
    pass