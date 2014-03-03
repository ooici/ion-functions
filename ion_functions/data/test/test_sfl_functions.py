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
from ion_functions.utils import fill_value


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

        #   [  V_s,   V_c,   T_s,  T_c,  T_u,  T_lc,     T]
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
        V_s = test_array[:, 0]
        V_c = test_array[:, 1]

        # set known outputs
        T = test_array[:, 6]

        # calculate the Vent Fluid Temperature
        Tout = sflfunc.sfl_trhph_vfltemp(V_s, V_c, a, b, c, d, e)

        # How'd we do?
        np.testing.assert_allclose(Tout, T, rtol=0.01, atol=0.0)

    def test_sfl_trhph_vfl_thermistor_temp(self):
        """
        Test sfl_trhph_vfl_thermistor_temp function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_SPEC_TRHPHTE_OOI.pdf)

        Implemented by Russell Desiderio, 28-Feb-2014
        """

        #   [  V_s,   V_c,   T_s,  T_c,  T_u,  T_lc,     T]
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
        V_s = test_array[:, 0]

        # set known outputs
        T_s = test_array[:, 2]

        # calculate the thermistor temperature
        Tout = sflfunc.sfl_trhph_vfl_thermistor_temp(V_s)

        # How'd we do?
        np.testing.assert_allclose(Tout, T_s, rtol=0.01, atol=0.0)

    def test_sfl_trhph_vflorp(self):
        """
        Test sfl_trhph_vflorp function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Vent Fluid Temperature from
            TRHPH. Document Control Number 1341-00150.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00150_Data_Product_SPEC_TRHPHTE_OOI.pdf)

        Implemented by Russell Desiderio, 28-Feb-2014
        """

        # test calibration constants
        offset = 2004.0
        gain = 4.0

        #   [  V      ORP[mV] ]
        test_array = np.array([
            [1.806,   -50],
            [1.541,  -116],
            [1.810,   -49],
            [0.735,  -317],
            [0.745,  -315],
            [0.715,  -322],
            [0.775,  -307],
            [0.799,  -301],
            [0.757,  -312],
            [0.542,  -366],
            [0.831,  -293],
            [0.867,  -284],
            [0.911,  -273]
        ])

        # set inputs
        V = test_array[:, 0]

        # set known outputs
        ORP = test_array[:, -1]

        # calculate the oxidation-reduction potential
        ORP_out = sflfunc.sfl_trhph_vflorp(V, offset, gain)
        ORP_out = np.round(ORP_out)

        # How'd we do?
        np.testing.assert_allclose(ORP_out, ORP, rtol=0.01, atol=0.0)

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
        Modified by Russell Desiderio, Feb. 28, 2014.
            Extended unit tests; see comments below.
        """
        # original unit test data from DPS:
        #    did not test v_r1 v_r2 v_r3 conditional.
        #    did not unambiguously test resistivity value out of range.
        #        (all of the test result fill_values could arise
        #         solely from temperature out-of-range)
        #test_array = np.array([
        #    [0.906, 4.095, 4.095, 11.8, 4.530, 0.2208, fill_value, fill_value],
        #    [0.890, 4.095, 4.095, 15.9, 4.450, 0.2247, fill_value, fill_value],
        #    [0.891, 4.095, 4.095, 3.2, 4.455, 0.2245, fill_value, fill_value],
        #    [0.184, 0.915, 4.064, 67.7, 0.915, 1.0929, fill_value, fill_value],
        #    [0.198, 1.002, 4.095, 75.8, 1.002, 0.9980, fill_value, fill_value],
        #    [0.172, 0.857, 4.082, 97.5, 0.857, 1.1669, fill_value, fill_value],
        #    [0.183, 0.926, 4.076, 95.0, 0.926, 1.0799, fill_value, fill_value],
        #    [0.233, 1.182, 4.072, 96.2, 1.182, 0.8460, fill_value, fill_value],
        #    [0.146, 0.747, 3.634, 116.8, 0.727, 1.3759, 0.19507, 195],
        #    [0.134, 0.681, 3.405, 272.8, 0.681, 1.4684, 0.10893, 109],
        #    [0.131, 0.673, 3.293, 325.8, 0.659, 1.5184, 0.12813, 128],
        #    [0.133, 0.678, 3.396, 330.0, 0.679, 1.4723, 0.12664, 127],
        #    [0.135, 0.681, 3.409, 333.4, 0.682, 1.4667, 0.12878, 129],
        #    [0.135, 0.681, 3.426, 333.2, 0.685, 1.4594, 0.12802, 128]
        #])

        # the first fill_value tests out-of-range resistivity value v_r1 which
        #     results in an out-of-range conductivity value given the curvature
        #     of the surface as a function of temperature; for other temperatures,
        #     this v_r1 value (0.380) is not out-of-range (cf 1st and 3rd rows).
        # the second fill_value tests an out-of-range temperature value (84.6).
        # the third fill_value tests an out_of_range resistivity value (v_r2=2.8).
        #
        #     v_r1    v_r2    v_r3   temp    chl [mmol/kg]
        test_array = np.array([
            [0.440,  4.095,  4.095,  105.4,    59.0],
            [0.380,  4.095,  4.095,  241.9,   fill_value],
            [0.320,  4.095,  4.095,  374.2,    60.0],
            [0.184,  0.915,  4.064,  105.4,   175.0],
            [0.198,  1.002,  4.095,  241.9,    71.0],
            [0.172,  0.857,  4.082,  374.2,   132.0],
            [0.183,  0.926,  4.076,   84.6,   fill_value],
            [0.233,  2.800,  4.072,  250.2,   fill_value],
            [0.146,  0.747,  3.634,  116.8,   195.0],
            [0.134,  0.681,  3.405,  272.8,   109.0],
            [0.131,  0.673,  3.293,  325.8,   128.0],
            [0.133,  0.678,  3.396,  330.0,   127.0],
            [0.135,  0.681,  2.000,  333.4,   239.0],
            [0.135,  0.681,  1.000,  333.2,   501.0]
        ])

        # set inputs
        V_R1 = test_array[:, 0]
        V_R2 = test_array[:, 1]
        V_R3 = test_array[:, 2]
        T = test_array[:, 3]

        # set known output
        Cl = test_array[:, -1]

        # calculate the Vent Fluid Temperature
        Clout = sflfunc.sfl_trhph_chloride(V_R1, V_R2, V_R3, T)

        ###########################################################################
        # How'd we do?
        np.testing.assert_allclose(Clout, Cl, rtol=0, atol=0)
        ###########################################################################

    def test_sfl_sflpres_l1(self):
                """
        Test the sfl_sflpres_l1 function.

        Value based on that described in DPS as available on Alfresco:

        OOI (2013). Data Product Specification for Seafloor Pressure from
        Sea-Bird SBE 26PLUS. Document Control Number 1341-00230.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00230_Data_Product_SPEC_SFLPRES_OOI.pdf)

        Implemented by Craig Risien, February 2014
        """

    # test inputs
    pressure_output = np.array([14.868])

    # expected outputs
    pressure_expected = np.array([10.2511])

    # compute chla and cdom values
    pressure_calc = np.zeros(1)
    for i in range(0, 1):
        pressure_calc[i] = sflfunc.sfl_sflpres_l1(pressure_output[i])

    # compare calculated results to expected
    np.testing.assert_allclose(pressure_calc, pressure_expected, rtol=0.0001, atol=0.0001)
