#!/usr/bin/env python

"""
@package ion_functions.test.sfl_functions
@file ion_functions/test/sfl_functions.py
@author Christopher Wingard
@brief Unit tests for sfl_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import sfl_functions as sflfunc

@attr('UNIT', group='func')
class TestSFLFunctionsUnit(BaseUnitTestCase):

    def test_sfl_trhph_vfltemp(self):
        """
        Test sfl_trhph_vfltemp function.

        Values based on those described in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_SPEC_TRHPHTE_OOI.pdf)
            
        Implemented by Christopher Wingard, April 2013
        """
        
        # calibration constants
        a = 1.98e-9
        b = -2.45e-6
        c = 9.28e-4
        d = -0.0888
        e = 0.731

        test_array = np.array([
            [1.506, 0.000, 12.01, 0.00, 12.0, -0.206, 11.8],
            [1.479, 0.015, 12.67, 3.68, 16.3, -0.483, 15.9],
            [1.926, 0.001, 2.47, 0.25, 2.7, 0.497, 3.2],
            [1.932, 0.274, 2.34, 67.12, 69.5, -1.735, 67.7],
            [1.927, 0.306, 2.45, 74.96, 77.4, -1.648, 75.8],
            [1.930, 0.393, 2.38, 96.27, 98.7, -1.162, 97.5],
            [1.929, 0.383, 2.40, 93.82, 96.2, -1.234, 95.0],
            [1.930, 0.388, 2.38, 95.05, 97.4, -1.199, 96.2],
            [1.930, 0.469, 2.38, 114.89, 117.3, -0.497, 116.8],
            [1.931, 1.077, 2.36, 263.83, 266.2, 6.579, 272.8],
            [1.926, 1.288, 2.47, 315.52, 318.0, 7.797, 325.8],
            [1.926, 1.305, 2.47, 319.69, 322.2, 7.847, 330.0],
            [1.928, 1.319, 2.43, 323.12, 325.5, 7.883, 333.4],
            [1.929, 1.318, 2.40, 322.87, 325.3, 7.880, 333.2]
        ])
        
        # set inputs
        V_s = test_array[:,0]
        V_c = test_array[:,1]

        # set known output
        T = test_array[:,6]
        
        # calculate the Vent Fluid Temperature
        Tout = sflfunc.sfl_trhph_vfltemp(V_s, V_c, a, b, c, d, e)
        
        # How'd we do?
        np.testing.assert_allclose(Tout, T, rtol=0.1, atol=0)


    def test_sfl_trhph_chloride(self):
        """
        Test sfl_trhph_chloride function.

        Values based on those described in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for Vent Fluid Chloride
            Concentration. Document Control Number 1341-00160.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_SPEC_TRHPHTE_OOI.pdf)
            
        Implemented by Christopher Wingard, April 2013
        """
        test_array = np.array([
            [0.906, 4.095, 4.095, 11.8, 4.530, 0.2208, -99999999., -99999999.],
            [0.890, 4.095, 4.095, 15.9, 4.450, 0.2247, -99999999., -99999999.],
            [0.891, 4.095, 4.095, 3.2, 4.455, 0.2245, -99999999., -99999999.],
            [0.184, 0.915, 4.064, 67.7, 0.915, 1.0929, -99999999., -99999999.],
            [0.198, 1.002, 4.095, 75.8, 1.002, 0.9980, -99999999., -99999999.],
            [0.172, 0.857, 4.082, 97.5, 0.857, 1.1669, -99999999., -99999999.],
            [0.183, 0.926, 4.076, 95.0, 0.926, 1.0799, -99999999., -99999999.],
            [0.233, 1.182, 4.072, 96.2, 1.182, 0.8460, -99999999., -99999999.],
            [0.146, 0.747, 3.634, 116.8, 0.727, 1.3759, 0.19507, 195],
            [0.134, 0.681, 3.405, 272.8, 0.681, 1.4684, 0.10893, 109],
            [0.131, 0.673, 3.293, 325.8, 0.659, 1.5184, 0.12813, 128],
            [0.133, 0.678, 3.396, 330.0, 0.679, 1.4723, 0.12664, 127],
            [0.135, 0.681, 3.409, 333.4, 0.682, 1.4667, 0.12878, 129],
            [0.135, 0.681, 3.426, 333.2, 0.685, 1.4594, 0.12802, 128]
        ])
        
        # set inputs
        V_R1 = test_array[:,0]
        V_R2 = test_array[:,1]
        V_R3 = test_array[:,2]
        T = test_array[:,3]

        # set known output
        Cl = test_array[:,7]
        
        # calculate the Vent Fluid Temperature
        Clout = sflfunc.sfl_trhph_chloride(V_R1, V_R2, V_R3, T)
        
        # How'd we do?
        np.testing.assert_allclose(Clout, Cl, rtol=0, atol=0)
