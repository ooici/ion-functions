#!/usr/bin/env python

"""
@package ion_functions.test.co2_functions
@file ion_functions/test/co2_functions.py
@author Christopher Wingard
@brief Unit tests for co2_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import co2_functions as co2func
from ion_functions.utils import fill_value


@attr('UNIT', group='func')
class Testpco2FunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        ###### Test data for PCO2W ######
        raw_strings = np.array([
            '*7E2705CBACEE7F007D007D0B2A00BF080500E00187034A008200790B2D00BE080600DE0C1406C98C',
            '*7E2704CBACEECB008000880B2900B2080600D300FB0263007F00890B2B00B4080700CE0C5106C884',
            '*7E2704CBACEF43007E00890B27014408070189045B0875007E00870B2B0140080201860C5506C601',
            '*7E2704CBACEFBB007E00820B2A042B0803051F16182785008000850B2A043C080405390C5506C5A9',
            '*7E2704CBACF033007F00840B28054D080606831CCC330A007F00850B290551080406800C5606C556',
            '*7E2704CBACF0AB007E00770B2804DE080605F31A672E9F008100790B2F04E0080705F70C5606C5EB',
            '*7E2704CBACF1230081007C0B2B00BF080800DB01BF0390007C00790B3000C7080B00EA0C5606C80A',
            '*7E2704CBACF19B008000750B2B01D20807023008310E94007E00730B2A01D2080502290C5706C0A2',
            '*7E2704CBACF213008000740B2A01D50808042F08501FC6007D00780B3201D8080A04350C5706C0A5',
            '*7E2704CBACF28B007F00710B2F0174080203570615189B008100730B2A0174080A03580C5706C125',
            '*7E2704CBACF303007E006B0B2C019B080803CC07001CBA007F006F0B300199080803D00C5706C4E3'
        ])

        # reagent constants (instrument and reagent bag specific)
        self.ea434 = np.ones(11) * 19706.
        self.ea620 = np.ones(11) * 34.
        self.eb434 = np.ones(11) * 3073.
        self.eb620 = np.ones(11) * 44327.
        self.calt = np.ones(11) * 16.5
        self.cala = np.ones(11) * 0.0459
        self.calb = np.ones(11) * 0.6257
        self.calc = np.ones(11) * -1.5406

        # expected outputs
        self.therm = np.array([18.8526, 18.8765, 18.9245, 18.9485,
                               18.9485, 18.9485, 18.8765, 19.0686,
                               19.0686, 19.0446, 18.9725])
        self.pco2 = np.array([fill_value, 294.1720, 311.3361, 319.0101,
                              319.8925, 319.8950, 305.8104, 317.9661,
                              284.3676, 280.2324, 280.0354])

        # parse the data strings
        self.mtype = np.zeros(11, dtype=np.int)
        self.light = np.zeros((11, 14), dtype=np.int)
        self.traw = np.zeros(11, dtype=np.int)
        for i in range(11):
            # parse the raw strings into subelements, such as the driver would
            # provide.
            s = raw_strings[i]
            self.mtype[i] = int((s[5:7]), 16)
            self.traw[i] = int((s[75:79]), 16)
            strt = 15
            step = 4
            for j in range(14):
                self.light[i, j] = int((s[strt:strt+step]), 16)
                strt += step

            if self.mtype[i] == 5:
                a434blnk = co2func.pco2_blank(self.light[i, 6])
                a620blnk = co2func.pco2_blank(self.light[i, 7])

        self.a434blnk = np.ones(11) * a434blnk
        self.a620blnk = np.ones(11) * a620blnk

    def test_pco2_pco2wat(self):
        """
        Test pco2_pco2wat function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00490.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """

        # compute the thermistor temperature in deg_C, derive blanks and then
        # calculate pco2.

        ### bulk case ###
        tout = co2func.pco2_thermistor(self.traw)
        pco2out = co2func.pco2_pco2wat(self.mtype, self.light, tout,
                                       self.ea434, self.eb434, self.ea620, self.eb620,
                                       self.calt, self.cala, self.calb, self.calc,
                                       self.a434blnk, self.a620blnk)

        np.testing.assert_allclose(pco2out, self.pco2, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(tout, self.therm, rtol=1e-4, atol=1e-4)

        ### single record case ###
        indx = 0
        for mtype in self.mtype:
            tout = co2func.pco2_thermistor(self.traw[indx])
            pco2out = co2func.pco2_pco2wat(mtype, self.light[indx, :], tout,
                                           self.ea434[indx], self.eb434[indx],
                                           self.ea620[indx], self.eb620[indx],
                                           self.calt[indx], self.cala[indx],
                                           self.calb[indx], self.calc[indx],
                                           self.a434blnk[indx], self.a620blnk[indx])

            np.testing.assert_allclose(pco2out, self.pco2[indx], rtol=1e-4, atol=1e-4)
            np.testing.assert_allclose(tout, self.therm[indx], rtol=1e-4, atol=1e-4)

            indx += 1

    def test_pco2_co2flux(self):
        """
        Test pco2_co2flux function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Flux of CO2 into the
            Atmosphere. Document Control Number 1341-00270.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00270_Data_Product_SPEC_CO2FLUX_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """
        test_data = np.array([
            [360, 390, 5, 0, 34, -2.063e-08],
            [360, 390, 5, 0, 35, -2.052e-08],
            [360, 390, 5, 10, 34, -1.942e-08],
            [360, 390, 5, 10, 35, -1.932e-08],
            [360, 390, 5, 20, 34, -1.869e-08],
            [360, 390, 5, 20, 35, -1.860e-08],
            [360, 390, 10, 0, 34, -8.250e-08],
            [360, 390, 10, 0, 35, -8.207e-08],
            [360, 390, 10, 10, 34, -7.767e-08],
            [360, 390, 10, 10, 35, -7.728e-08],
            [360, 390, 10, 20, 34, -7.475e-08],
            [360, 390, 10, 20, 35, -7.440e-08],
            [360, 390, 20, 0, 34, -3.300e-07],
            [360, 390, 20, 0, 35, -3.283e-07],
            [360, 390, 20, 10, 34, -3.107e-07],
            [360, 390, 20, 10, 35, -3.091e-07],
            [360, 390, 20, 20, 34, -2.990e-07],
            [360, 390, 20, 20, 35, -2.976e-07],
            [400, 390, 5, 0, 34, 6.875e-09],
            [400, 390, 5, 0, 35, 6.839e-09],
            [400, 390, 5, 10, 34, 6.472e-09],
            [400, 390, 5, 10, 35, 6.440e-09],
            [400, 390, 5, 20, 34, 6.229e-09],
            [400, 390, 5, 20, 35, 6.200e-09],
            [400, 390, 10, 0, 34, 2.750e-08],
            [400, 390, 10, 0, 35, 2.736e-08],
            [400, 390, 10, 10, 34, 2.589e-08],
            [400, 390, 10, 10, 35, 2.576e-08],
            [400, 390, 10, 20, 34, 2.492e-08],
            [400, 390, 10, 20, 35, 2.480e-08],
            [400, 390, 20, 0, 34, 1.100e-07],
            [400, 390, 20, 0, 35, 1.094e-07],
            [400, 390, 20, 10, 34, 1.036e-07],
            [400, 390, 20, 10, 35, 1.030e-07],
            [400, 390, 20, 20, 34, 9.966e-08],
            [400, 390, 20, 20, 35, 9.920e-08],
            [440, 390, 5, 0, 34, 3.438e-08],
            [440, 390, 5, 0, 35, 3.420e-08],
            [440, 390, 5, 10, 34, 3.236e-08],
            [440, 390, 5, 10, 35, 3.220e-08],
            [440, 390, 5, 20, 34, 3.114e-08],
            [440, 390, 5, 20, 35, 3.100e-08],
            [440, 390, 10, 0, 34, 1.375e-07],
            [440, 390, 10, 0, 35, 1.368e-07],
            [440, 390, 10, 10, 34, 1.294e-07],
            [440, 390, 10, 10, 35, 1.288e-07],
            [440, 390, 10, 20, 34, 1.246e-07],
            [440, 390, 10, 20, 35, 1.240e-07],
            [440, 390, 20, 0, 34, 5.500e-07],
            [440, 390, 20, 0, 35, 5.471e-07],
            [440, 390, 20, 10, 34, 5.178e-07],
            [440, 390, 20, 10, 35, 5.152e-07],
            [440, 390, 20, 20, 34, 4.983e-07],
            [440, 390, 20, 20, 35, 4.960e-07]
        ])

        # setup inputs and outputs
        pco2w = test_data[:, 0]
        pco2a = test_data[:, 1]
        u10 = test_data[:, 2]
        t = test_data[:, 3]
        s = test_data[:, 4]
        flux = test_data[:, 5]

        # compute the flux given the inputs
        out = co2func.pco2_co2flux(pco2w, pco2a, u10, t, s)

        # and compare the results
        self.assertTrue(np.allclose(out, flux, rtol=1e-9, atol=1e-9))
