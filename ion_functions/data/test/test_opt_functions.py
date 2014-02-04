#!/usr/bin/env python

"""
@package ion_functions.test.opt_functions
@file ion_functions/test/opt_functions.py
@author Christopher Wingard
@brief Unit tests for opt_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import opt_functions as optfunc

@attr('UNIT', group='func')
class TestOptFunctionsUnit(BaseUnitTestCase):

    def test_opt_functions(self):
        """
        Test the functions in the opt_functions.py module.

        Values based on test data in DPSs and available on Alfresco.
        
        OOI (2013). Data Product Specification for Optical Beam Attenuation
            Coefficient. Document Control Number 1341-00690.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)
        
        OOI (2013). Data Product Specification for Optical Absorption
            Coefficient. Document Control Number 1341-00700.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)
            
        Implemented by Christopher Wingard, April 2013
        """

        # test inputs
        tbins = np.array([14.5036, 15.5200, 16.4706, 17.4833, 18.4831, 19.5196, 20.5565])
        tarr = np.array([
            [-0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, -0.004355, -0.004131, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, -0.004131, -0.003422, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, -0.003422, -0.002442]
        ])
        
        sig = np.array([50, 150, 250, 350, 450, 495])
        ref = np.array([500, 500, 500, 500, 500, 500])
        traw = np.array([48750, 48355, 47950, 47535, 47115, 46684])
        wlngth = np.array([500., 550., 600., 650., 700., 715.])        
        T = np.array([4., 8., 12., 16., 20., 24.])
        PS = np.array([10., 15., 20., 25., 30., 35])
        
        # expected outputs
        tint = np.array([15., 16., 17., 18., 19., 20.])
        deltaT = np.array([-0.0048, -0.0045, -0.0044, -0.0042, -0.0038, -0.0030])
        cpd = np.array([9.2251, 4.8304, 2.7870, 1.4409, 0.4352, 0.0532])
        apd = np.array([9.3151, 4.9204, 2.8770, 1.5309, 0.5252, 0.1432])
        cpd_ts = np.array([9.2260, 4.8312, 2.7955, 1.4412, 0.4406, 0.0447])
        apd_ts = np.array([9.3155, 4.9206, 2.8848, 1.5303, 0.5297, 0.1338])
        apd_ts_s = np.array([9.1811, 4.7864, 2.7507, 1.3965, 0.3959, 0.0000])
        
        # compute beam attenuation and optical absorption values
        dgC = np.zeros(6)
        dT = np.zeros(6)
        c = np.zeros(6)
        c_ts = np.zeros(6)
        a = np.zeros(6)
        a_ts = np.zeros(6)
        a_ts_s = np.zeros(6)
        for i in range(6):
            dgC[i] = optfunc.opt_internal_temp(traw[i])
            c[i], dT[i] = optfunc.opt_pd_calc(ref[i], sig[i], 0.01, dgC[i], tbins, tarr[i])
            c_ts[i] = optfunc.opt_tempsal_corr('c', c[i], wlngth[i], 20., T[i], PS[i])           
            a[i], foo = optfunc.opt_pd_calc(ref[i], sig[i], 0.10, dgC[i], tbins, tarr[i])
            a_ts[i] = optfunc.opt_tempsal_corr('a', a[i], wlngth[i], 20., T[i], PS[i])
            
        a_ts_s = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth)
        
        # compare calculated results to expected
        np.testing.assert_allclose(dgC, tint, rtol=0.1, atol=0.1)
        np.testing.assert_allclose(dT, deltaT, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(c, cpd, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a, apd, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a_ts, apd_ts, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a_ts_s, apd_ts_s, rtol=1e-4, atol=1e-4)

    def test_opt_par_satlantic(self):
                """
        Test the opt_par_satlantic function.

        Values based on that described in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY ACTIVE RADIATION (PAR)
        FROM SATLANTIC INSTRUMENT ON RSN SHALLOW PROFILER
        Document Control Number 1341-00720. https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00720_Data_Product_SPEC_OPTPARW_Satl_OOI.pdf)
            
        Implemented by Craig Risien, February 2014
        """

   # test inputs
    Im = 1.3589
    a0 = 2156849800.8
    a1 = 2.586852835e-006  
    
    par_count = np.array([2159403328, 2159400384, 2159396992, 2159400384, 2159407488, 2159399296, 2159400384,
                          2159404800, 2159403904, 2159402240, 2159403200, 2159409728, 2159409792, 2159408320,
                          2159407808, 2159402304, 2159402688, 2159407552, 2159404160, 2159403776, 2159402048,
                          2159404544])

    # expected outputs
    par_expected = np.array([8.976348585, 8.965999618, 8.954075807, 8.965999618, 8.990972126, 8.962174999,
                             8.965999618, 8.981523069, 8.978373383, 8.972523967, 8.97589863, 8.998846341,
                             8.999071318, 8.993896835, 8.992097014, 8.972748944, 8.97409881, 8.991197104,
                             8.979273293, 8.977923428, 8.971849034, 8.980623159])
        
    # compute par values
    par_calc = np.zeros(22)
    for i in range(0,22):
        par_calc[i] = optfunc.opt_par_satlantic(par_count[i], a0, a1, Im)
               
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected, rtol=0.000001, atol=0.000001)
    
    
    def test_opt_par_biospherical(self):
                """
        Test the opt_par_biospherical function.

        Values based on that described in DPS as available on Alfresco:
        
        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY ACTIVE RADIATION (PAR)
        FROM SATLANTIC INSTRUMENT ON RSN SHALLOW PROFILER
        Document Control Number 1341-00721. https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00721_Data_Product_SPEC_OPTPARW_Bios_OOI.pdf)
            
        Implemented by Craig Risien, February 2014
        """

   # test inputs
    offset = 0.0101
    scale = 5.897e-04
    
    par_volts = np.array([1.016793, 0.599800, 0.452400, 0.305000, 0.187000, 0.178900, 0.069100, 0.039600, 0.010100])

    # expected outputs
    par_expected = np.array([1707.13, 1000.00, 750.00, 500.00, 300.00, 286.25, 100.00, 50.00, 0.00])
        
    # compute par values
    par_calc = np.zeros(9)
    for i in range(0,9):
        par_calc[i] = optfunc.opt_par_biospherical(par_volts[i], offset, scale)
               
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected, rtol=0.01, atol=0.01)