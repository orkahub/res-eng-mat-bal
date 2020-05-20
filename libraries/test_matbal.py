import sys
sys.path.insert(1, './')
import libraries.matbal as mb
from libraries.matbal import raisedError
import pytest
import numpy as np
import pandas as pd

# Formation_total_volume factor START:
def test_formation_total_volume_factor_higher():
    Bo, Bg, Rsb, Rs = 1.735, 0.6508, 1720, 1720
    fvft = mb.formation_total_volume_factor(Bo, Bg, Rsb, Rs)
    assert fvft == Bo

def test_formation_total_volume_factor_lower():
    Bo, Bg, Rsb, Rs = 1.735, 0.6508, 1720, 1720
    test_val = Bo + Bg * (Rsb - Rs)
    fvft = mb.formation_total_volume_factor(Bo, Bg, Rsb, Rs)
    assert fvft == test_val
# Formation_total_volume factor END:

#------------------------------------------------------
#--------------production injection balance START:
#------------------------------------------------------

#------TEST DIV by ZERO------#
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

#------RESOLVED------#

#------TEST NEGATIVE VALUES------#
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

#------RESOLVED------#

#    assert produced_oil_and_gas >= 0
#    assert produced_water >= 0
#    assert injected_water >= 0
#    assert injected_gas >= 0 
#    assert produced_oil_and_gas == 4456.8

#------------------------------------------------------
#--------------Production injection balance END:
#------------------------------------------------------


#------------------------------------------------------
#--------------New Func START:
#------------------------------------------------------


#------------------------------------------------------
#--------------New Func END:
#------------------------------------------------------
