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

    def test_co2_pco2wat(self):
        """
        Test co2_pco2wat function.

        Values based on those described in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00490.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
            
        Implemented by Christopher Wingard, April 2013
        """
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
        ea434 = 19706.
        ea620 = 34.
        eb434 = 3073.
        eb620 = 44327.
        calt = 16.5
        cala = 0.0459
        calb = 0.6257
        calc = -1.5406
        a434blnk = fill_value
        a620blnk = fill_value
        
        # expected outputs
        therm = np.array([18.8526, 18.8765, 18.9245, 18.9485,
                          18.9485, 18.9485, 18.8765, 19.0686,
                          19.0686, 19.0446, 18.9725])
        pco2 = np.array([fill_value, 294.1720, 311.3361, 319.0101,
                         319.8925, 319.8950, 305.8104, 317.9661,
                         284.3676, 280.2324, 280.0354
                         ])

        # parse the data strings
        light = np.zeros(14, dtype=np.int)
        pco2out = np.zeros(11)
        tout = np.zeros(11)
        vbout = np.zeros(11)
        for i in range(11):
            # parse the raw strings into subelements, such as the driver would
            # provide.
            s = raw_strings[i]
            mtype = int((s[5:7]), 16)
            traw = int((s[75:79]), 16)
            strt = 15; step = 4
            for j in range(14):
                light[j] = int((s[strt:strt+step]), 16)
                strt += step            
            
            # compute the thermistor temperature in deg_C, blanks and pco2
            tout[i] = co2func.pco2_thermistor(traw)
            a434blnk = co2func.pco2_abs434_blank(mtype, light, a434blnk)
            a620blnk = co2func.pco2_abs620_blank(mtype, light, a620blnk)
            pco2out[i] = co2func.pco2_pco2wat(mtype, light, tout[i], ea434,
                                               eb434, ea620, eb620, calt, cala,
                                               calb, calc, a434blnk, a620blnk)
        
        print pco2out
        self.assertTrue(np.allclose(pco2out, pco2, rtol=1e-4, atol=0))
        self.assertTrue(np.allclose(tout, therm, rtol=1e-4, atol=0))

