#!/usr/bin/env python

"""
@package ion_functions.test.opt_functions
@file ion_functions/test/opt_functions.py
@author Christopher Wingard, Russell Desiderio, Craig Risien
@brief Unit tests for opt_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import opt_functions as optfunc


@attr('UNIT', group='func')
class TestOptFunctionsUnit(BaseUnitTestCase):

    def test_opt_functions_OPTAA_sub_functions(self):
        """
        Test the OPTAA function subroutines in the opt_functions.py module.

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

        OOI (2014). OPTAA Unit Test. 1341-00700_OPTABSN Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00700_OPTABSN >>
            OPTAA_unit_test.xlsx)

        Implemented by Christopher Wingard, April 2013

        Modified by Russell Desiderio, 20-Feb-2014
            The opt_scatter_corr function was modified to trap out commonly
            occurring instances for when the scattering correction to absorption
            should not be applied (in these cases the additive correction term is
            set to 0). This requires that the unit test values for the abs and c
            data be physically reasonable (c >= abs). Therefore the unit test c
            data were changed by adding 1.0 to the c clear water offset.
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
        Tcal = 20.
        T = np.array([4., 8., 12., 16., 20., 24.])
        PS = np.array([10., 15., 20., 25., 30., 35])

        # clear water offsets
        c_off = 1.01  # previous tests used 0.01
        a_off = 0.1

        # also test case of user selected reference wavelength for abs scatter correction
        ref_wave = 695  # nearest wavelength (700nm) should become reference wavelength

        # expected outputs
        tint = np.array([15., 16., 17., 18., 19., 20.])
        deltaT = np.array([-0.0048, -0.0045, -0.0044, -0.0042, -0.0038, -0.0030])
        cpd = np.array([10.2251, 5.8304, 3.7870, 2.4409, 1.4352, 1.0532])
        apd = np.array([9.3151, 4.9204, 2.8770, 1.5309, 0.5252, 0.1432])
        cpd_ts = np.array([10.226025, 5.831245, 3.795494, 2.441203, 1.440651, 1.044652])
        apd_ts = np.array([9.3155, 4.9206, 2.8848, 1.5303, 0.5297, 0.1338])
        apd_ts_s_ref715 = np.array([9.181831, 4.786862, 2.751082, 1.396591, 0.396010, 0.000000])
        apd_ts_s_ref700 = np.array([8.785990, 4.390950, 2.355159, 1.000592, 0.000000, -0.396015])

        # compute beam attenuation and optical absorption values
        dgC = np.zeros(6)
        dT = np.zeros(6)
        c = np.zeros(6)
        c_ts = np.zeros(6)
        a = np.zeros(6)
        a_ts = np.zeros(6)
        a_ts_sdef = np.zeros(6)
        a_ts_s700 = np.zeros(6)
        for i in range(6):
            dgC[i] = optfunc.opt_internal_temp(traw[i])
            c[i], dT[i] = optfunc.opt_pd_calc(ref[i], sig[i], c_off, dgC[i], tbins, tarr[i])
            c_ts[i] = optfunc.opt_tempsal_corr('c', c[i], wlngth[i], Tcal, T[i], PS[i])
            a[i], foo = optfunc.opt_pd_calc(ref[i], sig[i], a_off, dgC[i], tbins, tarr[i])
            a_ts[i] = optfunc.opt_tempsal_corr('a', a[i], wlngth[i], Tcal, T[i], PS[i])

        # the scatter-correction-to-absorption check must be done outside the loop:
        #    for each iteration within the loop, abs and c for one wavelength are
        #    processed. because there is only one wavelength of data, that wavelength
        #    also becomes the reference wavelength, in which case the scatter-corrected
        #    abs value is calculated to be  a - (a/(c-a)) * (c-a) = identically 0 for
        #    each iteration wavelength.
        # case: default reference wavelength (715nm) used for scatter correction to abs.
        a_ts_sdef = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth)
        # case: user-selected reference wavelength (695nm --> 700nm is closest)
        a_ts_s700 = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth, ref_wave)

        # compare calculated results to expected
        np.testing.assert_allclose(dgC, tint, rtol=0.1, atol=0.1)
        np.testing.assert_allclose(dT, deltaT, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(c, cpd, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a, apd, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a_ts, apd_ts, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a_ts_sdef, apd_ts_s_ref715, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(a_ts_s700, apd_ts_s_ref700, rtol=1e-4, atol=1e-4)

        # now test the traps set for the cases in which unphysical scatter correction
        # calculations would be applied to the absorption data if not for the traps.
        a_ts_save = np.copy(a_ts)
        c_ts_save = np.copy(c_ts)

        # case: if a_ref < 0, do not apply the scatter correction to abs.
        # subcase: default reference wavelength
        a_ts[-1] = -0.01
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)
        # subcase: user selected reference wavelength
        a_ts = np.copy(a_ts_save)
        a_ts[-2] = -0.01
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth, ref_wave)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)

        # case: if a_ref > 0 but c_ref - a_ref < 0, do not apply the scatter correction to abs.
        # subcase: default reference wavelength
        a_ts = np.copy(a_ts_save)
        a_ts[-1] = 0.01
        c_ts[-1] = 0.005
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)
        # subcase: user selected reference wavelength
        a_ts = np.copy(a_ts_save)
        c_ts = np.copy(c_ts_save)
        a_ts[-2] = 0.01
        c_ts[-2] = 0.005
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth, ref_wave)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)

        # case: when both a_ref < 0 and c_ref - a_ref < 0, the scatter ratio does have
        # the correct sign; however, the scatter correction to absorption should not be
        # applied. this test is included because it is possible to code for the two
        # traps above without catching this case.
        # subcase: default reference wavelength
        a_ts = np.copy(a_ts_save)
        c_ts = np.copy(c_ts_save)
        a_ts[-1] = -0.01
        c_ts[-1] = -0.02
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)
        # subcase: user selected reference wavelength
        a_ts = np.copy(a_ts_save)
        c_ts = np.copy(c_ts_save)
        a_ts[-2] = -0.01
        c_ts[-2] = -0.02
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth, ref_wave)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)

        # case: when c_ref - a_ref = 0, do not apply the scatter correction to abs.
        # subcase: default reference wavelength
        a_ts = np.copy(a_ts_save)
        c_ts = np.copy(c_ts_save)
        a_ts[-1] = 0.01
        c_ts[-1] = 0.01
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)
        # subcase: user selected reference wavelength
        a_ts = np.copy(a_ts_save)
        c_ts = np.copy(c_ts_save)
        a_ts[-2] = 0.01
        c_ts[-2] = 0.01
        a_ts_out = np.zeros(6)
        a_ts_out = optfunc.opt_scatter_corr(a_ts, wlngth, c_ts, wlngth, ref_wave)
        np.testing.assert_allclose(a_ts_out, a_ts, rtol=1e-8, atol=1e-8)

    def test_opt_functions_OPTAA_wrapper_functions(self):
        """
        Test the OPTAA wrapper functions in the opt_functions.py module.
        Use realistically shaped ac-s data arrays in the calling arguments:
            offsets are dimensioned at the number of wavelengths;
            for internal temperature, L1 temperature, and salinity,
                one value per a-c dataset.

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

        OOI (2014). OPTAA Unit Test. 1341-00700_OPTABSN Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00700_OPTABSN >>
            OPTAA_unit_test.xlsx)

        Implemented by Russell Desiderio, 21-Feb-2014. Additional unit test data
            constructed as documented in the OPTAA Unit Test document artifact
            reference above.
        """

        # test inputs:
        # traw, T, and PS must be scalar.
        # a_off and c_off must have same dimensions as wlngth.
        tbins = np.array([14.5036, 15.5200, 16.4706, 17.4833, 18.4831, 19.5196, 20.5565])
        tarr = np.array([
            [0.0, -0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004355, -0.004131, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004131, -0.003422, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.003422, -0.002442, 0.0, 0.0, 0.0, 0.0]
        ])

        wlngth = np.array([500., 550., 600., 650., 700., 715.])
        c_sig = np.array([150., 225., 200., 350., 450., 495.])
        c_ref = np.array([550., 540., 530., 520., 510., 500.])
        c_off = np.array([1.35, 1.30, 1.25, 1.20, 1.15, 1.10])
        a_sig = np.array([250., 300., 210., 430., 470., 495.])
        a_ref = np.array([450., 460., 470., 480., 490., 500.])
        a_off = np.array([0.35, 0.30, 0.25, 0.20, 0.15, 0.10])
        traw = 48355.0
        Tcal = 20.0
        T = 12.0
        PS = 35.0

        # also test case of user selected reference wavelength for abs scatter correction
        ref_wave = 695  # nearest wavelength (700nm) should become reference wavelength

        # expected outputs
        cpd_ts = np.array([6.553646, 4.807949, 5.161655, 2.788220, 1.666362, 1.184530])
        apd_ts_s_ref715 = np.array([1.999989, 1.501766, 3.177048, 0.249944, 0.086438, 0.000000])
        apd_ts_s_ref700 = np.array([1.750859, 1.320885, 3.068470, 0.111075, 0.000000, -0.064806])

        # beam attenuation (beam c) wrapper test
        c_ts = np.zeros(6)
        c_ts = optfunc.opt_beam_attenuation(c_ref, c_sig, traw, wlngth, c_off, Tcal,
                                            tbins, tarr, T, PS)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-6, atol=1e-6)

        # absorption wrapper test
        # case: default reference wavelength for scatter correction
        a_ts_sdef = np.zeros(6)
        a_ts_sdef = optfunc.opt_optical_absorption(a_ref, a_sig, traw, wlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, wlngth, T, PS)
        np.testing.assert_allclose(a_ts_sdef, apd_ts_s_ref715, rtol=1e-6, atol=1e-6)
        # case: user selected reference wavelength for scatter correction
        a_ts_s700 = np.zeros(6)
        a_ts_s700 = optfunc.opt_optical_absorption(a_ref, a_sig, traw, wlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, wlngth, T, PS, ref_wave)
        np.testing.assert_allclose(a_ts_s700, apd_ts_s_ref700, rtol=1e-6, atol=1e-6)

    def test_opt_par_satlantic(self):
        """
        Test the opt_par_satlantic function.

        Values based on that described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY
        ACTIVE RADIATION (PAR) FROM SATLANTIC INSTRUMENT ON RSN SHALLOW
        PROFILER Document Control Number 1341-00720.
        https://alfresco.oceanobservatories.org/
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
    for i in range(0, 22):
        par_calc[i] = optfunc.opt_par_satlantic(par_count[i], a0, a1, Im)

    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected, rtol=0.000001, atol=0.000001)

    def opt_par_biospherical_mobile(self):
                """
        Test the opt_par_biospherical_mobile function.

        Values based on that described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY
        ACTIVE RADIATION (PAR) FROM BIOSPHERICAL INSTRUMENT QSP 2100
        ON CGSN MOBILE ASSETS Document Control Number 1341-00721.
        https://alfresco.oceanobservatories.org/
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
    for i in range(0, 9):
        par_calc[i] = optfunc.opt_par_biospherical_mobile(par_volts[i], offset, scale)

    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected, rtol=0.01, atol=0.01)

    def opt_par_biospherical_wfp(self):
                """
        Test the opt_par_biospherical_wfp function.

        Values based on that described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for PHOTOSYNTHETICALLY
        ACTIVE RADIATION (PAR) FROM BIOSPHERICAL INSTRUMENT QSP 2200
        ON CGSN PROFILERS Document Control Number 1341-00721.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00721_Data_Product_SPEC_OPTPARW_Bios_OOI.pdf)

        Implemented by Craig Risien, March 2014
        """

    # test inputs
    #offset units = mvolts
    offset = 10.1
    scale = 5.897e-04 / (6.02 * 10**13)

    par_mvolts = np.array([1016.793, 599.800, 452.400, 305.000, 187.000, 178.900, 69.100, 39.600, 10.100])

    # expected outputs
    par_expected = np.array([1707.13, 1000.00, 750.00, 500.00, 300.00, 286.25, 100.00, 50.00, 0.00])

    # compute par values
    par_calc = np.zeros(9)
    for i in range(0, 9):
        par_calc[i] = optfunc.opt_par_biospherical_wfp(par_mvolts[i], offset, scale)

    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected, rtol=0.01, atol=0.01)
