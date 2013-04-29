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
        tbins = np.array([14.490159, 15.493818, 16.489434, 17.510408, 18.507083, 19.499111, 20.49766])
        tarr = np.array([
            [0.001607, 0.001476, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.000308, -0.000235, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.000180, 0.000294, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, -0.000357, -0.000206, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, -0.000630, -0.000514, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, -0.000507, -0.000326]
        ])
        
        sig = np.array([150, 200, 250, 300, 350, 325])
        ref = np.array([500, 500, 500, 500, 500, 500])
        wlngth = np.array([500., 550., 600., 650., 700., 715.])        
        tint = np.array([15., 16., 17., 18., 19., 20.])
        coff = np.array([0.8, 0.9, 0.8, 0.8, 0.6, 0.1])
        aoff = np.array([0.7, 1.0, 1.0, 1.0, 1.0, 0.5])
        T = np.array([4., 8., 12., 16., 20., 24.])
        PS = np.array([10., 15., 20., 25., 30., 35])
        
        # expected outputs
        deltaT = np.array([0.0015, -0.0003, 0.0002, -0.0003, -0.0006, -0.0004])
        cpd = np.array([5.6144, 4.5654, 3.5724, 2.8436, 2.0273, 1.8235])
        apd = np.array([5.5144, 4.6654, 3.7724, 3.0436, 2.4273, 2.2235])
        cpd_ts = np.array([5.6153, 4.5663, 3.5809, 2.8439, 2.0327, 1.8150])
        apd_ts = np.array([5.5148, 4.6656, 3.7802, 3.0430, 2.4318, 2.2141])
        apd_ts_s = np.array([6.0724, 4.1146, 2.6745, 1.9382, 0.2176, 0.0000])
        
        # compute beam attenuation and optical absorption values
        dT = np.zeros(6)
        c = np.zeros(6)
        c_ts = np.zeros(6)
        a = np.zeros(6)
        a_ts = np.zeros(6)
        a_ts_s = np.zeros(6)
        for i in range(6):            
            c[i], dT[i] = optfunc.opt_pd_calc(ref[i], sig[i], coff[i], tint[i], tbins, tarr[i])
            c_ts[i] = optfunc.opt_tempsal_corr('c', c[i], wlngth[i], 20., T[i], PS[i])           
            a[i], foo = optfunc.opt_pd_calc(ref[i], sig[i], aoff[i], tint[i], tbins, tarr[i])
            a_ts[i] = optfunc.opt_tempsal_corr('a', a[i], wlngth[i], 20., T[i], PS[i])
            
        a_ts_s = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth)
        
        # compare calculated results to expected
        np.testing.assert_allclose(dT, deltaT, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(c, cpd, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a, apd, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a_ts, apd_ts, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a_ts_s, apd_ts_s, rtol=1e-4, atol=1e-4)
