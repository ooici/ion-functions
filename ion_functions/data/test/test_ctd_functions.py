#!/usr/bin/env python

"""
@package ion_functions.test.ctd_functions
@file ion_functions/test/ctd_functions.py
@author Christopher Wingard
@brief Unit tests for ctd_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase
from ion_functions.data import ctd_functions as ctdfunc
import numpy as np
import inspect
from pygsw import vectors as gsw
#import sys


@attr('UNIT', group='func')
class TestCTDFunctionsUnit(BaseUnitTestCase):

    #def test_ctdfunc_isnan(self):
    #    """
    #    Test to ensure functions return a Nan, if inputs are NaN.
    #
    #    Initial code by Luke Campbell, 2013-05-10
    #    Implemented by Christopher Wingard, 2013-05-10
    #    Generalized function handling by Russell Desiderio, 2014-02-03
    #    Modified by Russell Desiderio, 2014-02-18:
    #        Could not import pyon.util.log function, so to keep stderr
    #        from writing to automatic builds I commented it out.
    #    """
    #    #sys.stderr.write('\n\ntest nan inputs to function:\n\n')
    #    functions = inspect.getmembers(ctdfunc, inspect.isfunction)
    #    for i in range(len(functions)):
    #        fname = functions[i][0]
    #        f = functions[i][1]
    #        argspec = inspect.getargspec(f)
    #        retval = f(*[np.nan for i in argspec.args])
    #        #stringout = fname + ':  ' + str(retval) + '\n'
    #        #sys.stderr.write(stringout)
    #        self.assertTrue(np.isnan(retval))

    def test_ctd_sbe37im_tempwat_instrument_recovered(self):
        """
        Test ctd_sbe37im_tempwat_instrument_recovered function.

        Values based on matlab script in the ctd folder in the matlab scripts folder:
        calc_sbe37im_instrument_recovered_unit_test.m

        Implemented by Russell Desiderio, June 16, 2016
        """
        # test inputs
        t0 = np.array([366964, 499888, 465784, 500403])
        a0 = -1.179278E-04
        a1 = 3.097942E-04
        a2 = -4.688854E-06
        a3 = 2.081274E-07
        output = ctdfunc.ctd_sbe37im_tempwat_instrument_recovered(t0, a0, a1, a2, a3)
        # check values are from SBE data processing
        check_values = np.array([10.9818000, 3.8488000, 5.4520000, 3.8255000])
        # check tolerance to 0.1 millidegree C, which is the precision of these SBE values.
        np.testing.assert_allclose(output, check_values, rtol=0, atol=1.e-4)

    def test_ctd_sbe16plus_tempwat(self):
        """
        Test ctd_sbe16plus_tempwat function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """
        # test inputs
        t0 = 248471
        a0 = 1.281651e-3
        a1 = 2.706002e-4
        a2 = -1.027561e-6
        a3 = 1.749446e-7
        tout = ctdfunc.ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3)
        np.testing.assert_allclose(tout, 22.544681, rtol=1e-6, atol=0)

    def test_ctd_sbe37im_tempwat(self):
        """
        Test ctd_sbe37im_tempwat function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)

        Implemented by Russell Desiderio, February 5, 2014
        """
        # test input
        t0 = 340357.0

        tout = ctdfunc.ctd_sbe37im_tempwat(t0)
        np.testing.assert_allclose(tout, 24.035700, rtol=1e-6, atol=0)

    def test_ctd_sbe52mp_tempwat(self):
        """
        Test ctd_sbe52mp_tempwat function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)

        Implemented by Russell Desiderio, February 17, 2014
        """
        # test input
        t0 = 200000.0

        tout = ctdfunc.ctd_sbe52mp_tempwat(t0)
        np.testing.assert_allclose(tout, 15.000000, rtol=1e-6, atol=0)

    def test_ctd_sbe37im_preswat_instrument_recovered(self):
        """
        Test ctd_sbe37im_preswat_instrument_recovered function.

        Values based on matlab script in the ctd folder in the matlab scripts folder:
        calc_sbe37im_instrument_recovered_unit_test.m

        Implemented by Russell Desiderio, June 16, 2016
        """
        p0 = np.array([533152, 571309, 632465, 828170])
        therm0 = np.array([1608, 1452, 1471, 1453])
        pa0 = 1.202594e-1
        pa1 = 4.514834e-3
        pa2 = -1.091899e-11
        ptempa0 = -6.953022e1
        ptempa1 = 5.115592e-2
        ptempa2 = -3.918145e-7
        ptca0 = 5.247204e5
        ptca1 = 9.617295e-1
        ptca2 = 6.296724e-3
        ptcb0 = 2.498163e1
        ptcb1 = -2.75e-4
        ptcb2 = 0

        output = ctdfunc.ctd_sbe37im_preswat_instrument_recovered(
                     p0, therm0, ptempa0, ptempa1, ptempa2, ptca0,
                     ptca1, ptca2, ptcb0, ptcb1, ptcb2, pa0, pa1, pa2)
        # check values are from SBE data processing
        check_values = np.array([16.1590, 134.9500, 325.2580, 933.8820])
        # check tolerance to about 5mm, which is about the precision of these SBE values.
        np.testing.assert_allclose(output, check_values, rtol=0, atol=0.005)

    def test_ctd_sbe16plus_preswat(self):
        """
        Test ctd_sbe16plus_preswat function.

        Values based on those defined in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """
        p0 = 528418
        therm0 = 24303
        ptempa0 = -6.87701000e+001
        ptempa1 = 5.05406200e+001
        ptempa2 = -2.15672900e-001
        ptca0 = 5.24965500e+005
        ptca1 = 7.23620100e+000
        ptca2 = -9.94485900e-002
        ptcb0 = 2.51220000e+001
        ptcb1 = -2.00000000e-004
        ptcb2 = 0.00000000e+000
        pa0 = 1.73472300e+000
        pa1 = 1.57475000e-002
        pa2 = -6.51927800e-010

        p = ctdfunc.ctd_sbe16plus_preswat(p0, therm0, ptempa0, ptempa1, ptempa2,
                                          ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                                          pa0, pa1, pa2)
        np.testing.assert_allclose(p, 27.282116, rtol=1e-6, atol=0)

    def test_ctd_sbe16digi_preswat(self):
        """
        Test ctd_sbe16digi_preswat function.

        Values based on those defined in the matlab script:
            construct_sbe_digiquartz_pressure_unit_test.m
            available at: ion-functions/ion_functions/data/matlab_scripts/ctd/

        Implemented by Russell Desiderio, February 2, 2014
        """
        p0 = 8833629.0
        t0 = 34336.0
        C1 = 991.3651
        C2 = 1.01360e-05
        C3 = -1.18210e-04
        D1 = 0.031072
        D2 = 0.0
        T1 = 27.67412
        T2 = -1.08033e-04
        T3 = 1.03670e-06
        T4 = 1.68749e-09
        T5 = 0.0

        p = ctdfunc.ctd_sbe16digi_preswat(p0, t0, C1, C2, C3, D1, D2, T1, T2, T3, T4, T5)
        np.testing.assert_allclose(p, 49.999967, rtol=1e-6, atol=0)

    def test_ctd_sbe37im_preswat(self):
        """
        Test ctd_sbe37im_preswat function.

        Values based on those defined in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)

        Implemented by Russell Desiderio, February 5, 2014
        """
        p0 = 2789.0
        p_range_psia = 1000.0

        p = ctdfunc.ctd_sbe37im_preswat(p0, p_range_psia)
        np.testing.assert_allclose(p, 0.04536611, rtol=1e-6, atol=0)

    def test_ctd_glider_preswat(self):
        """
        Test ctd_glider_preswat function.

        Implemented by Russell Desiderio, October 28, 2015
        """
        p0 = np.array([0.9, 5.6789, 11.07, 121.212, 1234.5678])
        xpctd = np.array([9.0, 56.789, 110.7, 1212.12, 12345.678])

        p = ctdfunc.ctd_glider_preswat(p0)
        np.testing.assert_allclose(p, xpctd, rtol=1e-8, atol=1e-8)

    def test_ctd_sbe52mp_preswat(self):
        """
        Test ctd_sbe52mp_preswat function.

        Values based on those defined in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)

        Implemented by Russell Desiderio, February 17, 2014
        """
        p0 = 201000.0

        p = ctdfunc.ctd_sbe52mp_preswat(p0)
        np.testing.assert_allclose(p, 2000.000000, rtol=1e-6, atol=0)

    def test_ctd_sbe37im_condwat_instrument_recovered(self):
        """
        Test ctd_sbe37im_condwat_instrument_recovered function.

        Values based on matlab script in the ctd folder in the matlab scripts folder:
        calc_sbe37im_instrument_recovered_unit_test.m

        Implemented by Russell Desiderio, June 16, 2016
        """
        c0 = np.array([1564991, 1457279, 1484332, 1462659])
        t1 = np.array([10.9818000, 3.8488000, 5.4520000, 3.8255000])
        p1 = np.array([16.1590, 134.9500, 325.2580, 933.8820])
        g = -9.899853E-01
        h = 1.314100E-01
        i = -4.181710E-04
        j = 4.723872E-05
        cpcor = -9.570000E-08
        ctcor = 3.250000E-06
        wbotc = 4.842900E-07

        output = ctdfunc.ctd_sbe37im_condwat_instrument_recovered(c0, t1, p1, g, h, i, j,
                                                                  cpcor, ctcor, wbotc)
        # check values are from SBE data processing
        check_values = np.array([3.891373, 3.240767, 3.399795, 3.272396])
        # check tolerance to the 6th decimal place, the precision of these SBE values.
        np.testing.assert_allclose(output, check_values, rtol=0, atol=1e-6)

    def test_ctd_sbe16plus_condwat(self):
        """
        Test ctd_sbe16plus_condwat function.

        Values based on those defined in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_PRACSAL_OOI.pdf)

        Implemented by Christopher Wingard, March 2013
        """
        c0 = 1673175
        t1 = 22.544681
        p1 = 27.282116
        g = -9.72193700e-001
        h = 1.38675900e-001
        i = -1.08398500e-004
        j = 2.63219300e-005
        cpcor = -9.57000000e-008
        ctcor = 3.2500e-006

        c = ctdfunc.ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor)
        np.testing.assert_allclose(c, 4.969069, rtol=1e-6, atol=0)

    def test_ctd_sbe37im_condwat(self):
        """
        Test ctd_sbe37im_condwat function.

        Values based on those defined in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_PRACSAL_OOI.pdf)

        Implemented by Russell Desiderio, February 5, 2014
        """
        c0 = 400000.0

        c = ctdfunc.ctd_sbe37im_condwat(c0)
        np.testing.assert_allclose(c, 3.500000, rtol=1e-6, atol=0)

    def test_ctd_sbe52mp_condwat(self):
        """
        Test ctd_sbe52mp_condwat function.

        Values based on those defined in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_PRACSAL_OOI.pdf)

        Implemented by Russell Desiderio, February 17, 2014
        """
        c0 = 305000.0

        c = ctdfunc.ctd_sbe52mp_condwat(c0)
        np.testing.assert_allclose(c, 3.000000, rtol=1e-6, atol=0)

    def test_ctd_pracsal(self):
        """
        Test ctd_pracsal function.

        Values based on those defined in DPS:

        OOI (2012). Data Product Specification for Salinty. Document Control
            Number 1341-00040. https://alfresco.oceanobservatories.org/ (See: 
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_PRACSAL_OOI.pdf)

        Implemented by Christopher Wingard, March 2013
        """

        c = np.array([5.407471, 5.407880, 5.041008, 3.463402, 3.272557, 3.273035])
        t = np.array([28., 28., 20., 6., 3., 2.])
        p = np.array([0., 10., 150., 800., 2500., 5000.])

        output = ctdfunc.ctd_pracsal(c, t, p)

        """
        Note, DPS rounds off output values to %.1f. For test to work, these were
        recalculated using the GSW Toolbox, Version 3.02 in Matlab R2013a and
        output using %.6f (see Matlab code snippet below). The DPS will be
        editted to correctly specify the higher precision.

        >> sprintf('%.6f\t',gsw_SP_from_C(c*10,t,p))
        ans =
        33.495229	33.495224	36.995774	34.898526	34.999244	34.999494
        """
        check_values = np.array([33.495229,
                                 33.495224,
                                 36.995774,
                                 34.898526,
                                 34.999244,
                                 34.999494])
        np.testing.assert_allclose(output, check_values, rtol=1e-6, atol=0)

    def test_ctd_density(self):
        """
        Test ctd_density function.

        Values based on those defined in DPS:

        OOI (2012). Data Product Specification for Density. Document Control
            Number 1341-00050. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf)

        Implemented by Christopher Wingard, March 2013
        """

        SP = np.array([33.5, 33.5, 37, 34.9, 35, 35])
        t = np.array([28., 28., 20., 6., 3., 2.])
        p = np.array([0., 10., 150., 800., 2500., 5000.])
        lat = 15.00
        lon = -55.00

        output = ctdfunc.ctd_density(SP, t, p, lat, lon)

        check_values = np.array([1021.26851,
                                 1021.31148,
                                 1026.94422,
                                 1031.13498,
                                 1039.28768,
                                 1050.30616])
        np.testing.assert_allclose(output, check_values, rtol=1e-6, atol=0)
