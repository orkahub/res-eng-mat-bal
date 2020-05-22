import sys
sys.path.insert(1, './')
import libraries.matbal as mb
from libraries.matbal import raisedError
import pytest
import numpy as np
import pandas as pd



#------------------------------------------------------
#--------------production injection balance START:
#------------------------------------------------------

def test_production_injection_balance_zero_div():
    """ Test to see if Np is zero""" 
    
    # Arrange
    Np = 0
    Bt = 1.735
    Rs = 1720
    Rsi = 1720
    Bg = 0.00064
    Wp = 0
    Bw = 1.0
    Winj = 0
    Bwinj = 1 
    Ginj = 0
    Bginj = 0.00064
    Gp = 3000000

    # Act & Assert
    with pytest.raises(raisedError) as context:
        f, produced_oil_and_gas, produced_water, injected_gas, injected_water = mb.production_injection_balance(Np,
                                                                                 Bt, Rs, Rsi, Bg, Wp, Bw, 
                                                                                 Winj, Bwinj, Ginj, Bginj, Gp)

def test_production_injection_balance_negatives():
    """Test for negative reservoir volume inputs""" 

    # Arrange
    Np = 4000
    Bt = 1.735
    Rs = 1720
    Rsi = 1720
    Bg = 0.00064
    Wp = 0
    Bw = 1.0
    Winj = 0
    Bwinj = 1 
    Ginj = 0
    Bginj = 0.00064
    Gp = -3000000

    # Act

    with pytest.raises(raisedError) as context:
        f, produced_oil_and_gas, produced_water, injected_gas, injected_water = mb.production_injection_balance(Np, 
                                                                                 Bt, Rs, Rsi, Bg, Wp, Bw, 
                                                                                 Winj, Bwinj, Ginj, Bginj, Gp)

#------------------------------------------------------
#--------------gas_cap_expansion START:
#------------------------------------------------------


def test_gas_cap_expansion_zeroDiv():
    """Test for negative division"""
    # Arrange
    Bti = None
    Bg = 1.2
    Bgi = 0
    
    # Assert act
    with pytest.raises(raisedError) as context:
        Eg = mb.gas_cap_expansion(Bti, Bg, Bgi)



#------------------------------------------------------
#--pore_volume_reduction_connate_water_expansion START:
#------------------------------------------------------
def test_pore_volume_reduction_connate_water_expansion_ZeroDiv():
    """Check that Swi != 1.0 && deltaP != 0"""

    # Arrange
    m = 1
    Boi = 1
    cw = 1
    Swi = 1
    cf = 1
    deltaP = 0

    # Act
    with pytest.raises(raisedError) as context:
        Efw = mb.pore_volume_reduction_connate_water_expansion(m, Boi, cw, Swi, cf, deltaP)
    




