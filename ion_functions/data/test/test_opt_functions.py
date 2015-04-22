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
        2015-04-17: Russell Desiderio. Changed signal and reference count input from float to fix.
        """

        # test inputs:
        tbins = np.array([14.5036, 15.5200, 16.4706, 17.4833, 18.4831, 19.5196, 20.5565])
        # nrows of tarr = length(wlngth); ncols of tarr = length(tbins)
        tarr = np.array([
            [0.0, -0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004355, -0.004131, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004131, -0.003422, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.003422, -0.002442, 0.0, 0.0, 0.0, 0.0]
        ])

        # a_off and c_off must have same dimensions as wlngth.
        wlngth = np.array([500., 550., 600., 650., 700., 715.])
        c_sig = np.array([150, 225, 200, 350, 450, 495])
        c_ref = np.array([550, 540, 530, 520, 510, 500])
        c_off = np.array([1.35, 1.30, 1.25, 1.20, 1.15, 1.10])
        a_sig = np.array([250, 300, 210, 430, 470, 495])
        a_ref = np.array([450, 460, 470, 480, 490, 500])
        a_off = np.array([0.35, 0.30, 0.25, 0.20, 0.15, 0.10])

        # traw, T, and PS must be scalar (before vectorization).
        traw = 48355
        Tcal = 20.0
        T = 12.0
        PS = 35.0

        # also test case of user selected reference wavelength for abs scatter correction
        ref_wave = 695  # nearest wavelength (700nm) should become reference wavelength

        # expected outputs
        cpd_ts = np.array([6.553646, 4.807949, 5.161655, 2.788220, 1.666362, 1.184530], ndmin=2)
        apd_ts_s_ref715 = np.array([1.999989, 1.501766, 3.177048, 0.249944, 0.086438, 0.000000],
                                   ndmin=2)
        apd_ts_s_ref700 = np.array([1.750859, 1.320885, 3.068470, 0.111075, 0.000000, -0.064806],
                                   ndmin=2)

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

    def test_opt_functions_OPTAA_a_and_c_wavelength_sets(self):
        """
        Test the OPTAA wrapper functions in the opt_functions.py module.
            Test a-c dataset with different wavelength values for the absorption versus
            beam c optical channels. Included to test that the c optical values are correctly
            interpolated onto the abs wavelengths in the scatter correction algorithm.

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

        Implemented by Russell Desiderio, 26-Feb-2014. Additional unit test data
            constructed as documented in the OPTAA Unit Test document artifact
            reference above.
        2015-04-17: Russell Desiderio. Changed signal and reference count input from float to fix.
        """

        # test inputs:
        tbins = np.array([14.5036, 15.5200, 16.4706, 17.4833, 18.4831, 19.5196, 20.5565])
        # nrows of tarr = length(wlngth); ncols of tarr = length(tbins)
        tarr = np.array([
            [0.0, -0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004355, -0.004131, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004131, -0.003422, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.003422, -0.002442, 0.0, 0.0, 0.0, 0.0]
        ])

        # a_off and c_off must have same dimensions as wlngth.
        cwlngth = np.array([510., 540., 580., 630., 670., 710.])
        awlngth = np.array([500., 550., 600., 650., 700., 715.])
        c_sig = np.array([150, 225, 200, 350, 450, 495])
        c_ref = np.array([550, 540, 530, 520, 510, 500])
        c_off = np.array([1.35, 1.30, 1.25, 1.20, 1.15, 1.10])
        a_sig = np.array([250, 300, 210, 430, 470, 495])
        a_ref = np.array([450, 460, 470, 480, 490, 500])
        a_off = np.array([0.35, 0.30, 0.25, 0.20, 0.15, 0.10])

        # traw, T, and PS must be scalar (before vectorization).
        traw = 48355
        Tcal = 20.0
        T = 12.0
        PS = 35.0

        # also test case of user selected reference wavelength for abs scatter correction
        ref_wave = 695  # nearest abs wavelength (700nm) should become reference wavelength

        # expected outputs
        cpd_ts = np.array([6.553771, 4.807914, 5.156010, 2.788715, 1.655607, 1.171965], ndmin=2)
        apd_ts_s_ref715 = np.array([1.990992, 1.479089, 3.350109, 0.350108, 0.152712, 0.000000],
                                   ndmin=2)
        apd_ts_s_ref700 = np.array([1.379858, 1.021574, 3.235057, 0.099367, 0.000000, -0.156972],
                                   ndmin=2)

        # beam attenuation (beam c) wrapper test
        c_ts = np.zeros(6)
        c_ts = optfunc.opt_beam_attenuation(c_ref, c_sig, traw, cwlngth, c_off, Tcal,
                                            tbins, tarr, T, PS)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-6, atol=1e-6)

        # absorption wrapper test
        # case: default reference wavelength for scatter correction
        a_ts_sdef = np.zeros(6)
        a_ts_sdef = optfunc.opt_optical_absorption(a_ref, a_sig, traw, awlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, cwlngth, T, PS)
        np.testing.assert_allclose(a_ts_sdef, apd_ts_s_ref715, rtol=1e-6, atol=1e-6)
        # case: user selected reference wavelength for scatter correction
        a_ts_s700 = np.zeros(6)
        a_ts_s700 = optfunc.opt_optical_absorption(a_ref, a_sig, traw, awlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, cwlngth, T, PS, ref_wave)
        np.testing.assert_allclose(a_ts_s700, apd_ts_s_ref700, rtol=1e-6, atol=1e-6)

    def test_opt_functions_OPTAA_vectorization(self):
        """
        Test the vectorization of the OPTAA wrapper functions in the opt_functions.py module.

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

        Implemented by Russell Desiderio, 05-Mar-2014. Additional unit test data
            constructed as documented in the OPTAA Unit Test document artifact
            reference above.
        2015-04-17: Russell Desiderio. Changed signal and reference count input from float to fix.
                                       Replaced 'exec' statements.
        2015-04-21: Russell Desiderio. Corrected input data shapes to current CI implementation.
                                       Added scalar time record case.
        """
        ### set test inputs: scalar (that is, 1) time record case
        # native scalars will come into the DPA with shape (1,)
        traw = np.array([48355])
        Tcal = np.array([20.0])
        T = np.array([12.0])
        PS = np.array([35.0])
        # native 1D arrays will come in as 2D row vectors
        cwlngth = np.array([[510., 540., 580., 630., 670., 710.]])
        awlngth = np.array([[500., 550., 600., 650., 700., 715.]])
        c_sig = np.array([[150, 225, 200, 350, 450, 495]])
        c_ref = np.array([[550, 540, 530, 520, 510, 500]])
        c_off = np.array([[1.35, 1.30, 1.25, 1.20, 1.15, 1.10]])
        a_sig = np.array([[250, 300, 210, 430, 470, 495]])
        a_ref = np.array([[450, 460, 470, 480, 490, 500]])
        a_off = np.array([[0.35, 0.30, 0.25, 0.20, 0.15, 0.10]])
        tbins = np.array([[14.5036, 15.5200, 16.4706, 17.4833, 18.4831, 19.5196, 20.5565]])
        # native 2D arrays will come in as 3D arrays, for this single time record case;
        # a singleton dimension will be prepended to the native 2D array.
        # native dimensions of tarr are (nwavelengths, ntempbins),
        # so that tarr.shape = (1,6,7)
        tarr = np.array([[
            [0.0, -0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004355, -0.004131, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004131, -0.003422, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.003422, -0.002442, 0.0, 0.0, 0.0, 0.0]
        ]])

        ## expected outputs to scalar time record test: 2D row vectors
        # beam attenuation (beam c) test
        # expected:
        cpd_ts = np.array([[6.553771, 4.807914, 5.156010, 2.788715, 1.655607, 1.171965]])
        # calculated:
        c_ts = optfunc.opt_beam_attenuation(c_ref, c_sig, traw, cwlngth, c_off, Tcal,
                                            tbins, tarr, T, PS)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-6, atol=1e-6)

        # absorption test
        # case: default reference wavelength for scatter correction
        apd_ts_s_ref715 = np.array([[1.990992, 1.479089, 3.350109, 0.350108, 0.152712, 0.000000]])
        a_ts_sdef = optfunc.opt_optical_absorption(a_ref, a_sig, traw, awlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, cwlngth, T, PS)
        np.testing.assert_allclose(a_ts_sdef, apd_ts_s_ref715, rtol=1e-6, atol=1e-6)

        # absorption test
        # case: user selected reference wavelength for scatter correction to abs values
        ref_wave = 695  # nearest abs wavelength (700nm) should become reference wavelength

        apd_ts_s_ref700 = np.array([[1.379858, 1.021574, 3.235057, 0.099367, 0.000000, -0.156972]])
        a_ts_s700 = optfunc.opt_optical_absorption(a_ref, a_sig, traw, awlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, cwlngth, T, PS, ref_wave)
        np.testing.assert_allclose(a_ts_s700, apd_ts_s_ref700, rtol=1e-6, atol=1e-6)

        ### multiple time records case
        # replicate the inputs to represent 3 data packets
        npackets = 3
        [traw, Tcal, T, PS] = [np.tile(xarray, npackets) for xarray in [traw, Tcal, T, PS]]
        [cwlngth, awlngth, c_sig, c_ref, c_off, a_sig, a_ref, a_off, tbins] = [
            np.tile(xarray, (npackets, 1)) for xarray in [
                cwlngth, awlngth, c_sig, c_ref, c_off, a_sig, a_ref, a_off, tbins]]
        tarr = np.tile(tarr, (npackets, 1, 1))

        # replicate the expected output data products
        cpd_ts = np.tile(cpd_ts, (npackets, 1))
        apd_ts_s_ref715 = np.tile(apd_ts_s_ref715, (npackets, 1))
        apd_ts_s_ref700 = np.tile(apd_ts_s_ref700, (npackets, 1))

        # beam attenuation (beam c) test
        c_ts = optfunc.opt_beam_attenuation(c_ref, c_sig, traw, cwlngth, c_off, Tcal,
                                            tbins, tarr, T, PS)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-6, atol=1e-6)

        # absorption test
        # case: default reference wavelength for scatter correction
        a_ts_sdef = optfunc.opt_optical_absorption(a_ref, a_sig, traw, awlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, cwlngth, T, PS)
        np.testing.assert_allclose(a_ts_sdef, apd_ts_s_ref715, rtol=1e-6, atol=1e-6)

        # case: user selected reference wavelength for scatter correction
        a_ts_s700 = optfunc.opt_optical_absorption(a_ref, a_sig, traw, awlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, cwlngth, T, PS, ref_wave)
        np.testing.assert_allclose(a_ts_s700, apd_ts_s_ref700, rtol=1e-6, atol=1e-6)

    def test_opt_functions_OPTAA_vectorization_with_tscor_nans(self):
        """
        Test the vectorized OPTAA wrapper functions in the opt_functions.py module;
        include "a" and "c" test wavelengths outside of the (original) range of the
        wavelength keys in the tscor.py file, which contains the dictionary of the
        temperature and salinity correction coefficients. The dictionary keys have
        been extended from [400.0 755.0] to [380.0 775.0] with entries of np.nan.

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

        Initial code by Russell Desiderio, 29-Apr-2014. Added tests for when OPTAA
            wavelengths are outside the range of the empirically derived temperature
            and salinity correction wavelength keys in the original tscor.py dictionary
            by modifying def test_opt_functions_OPTAA_vectorization test data.
        2015-04-17: Russell Desiderio. Use np.nan instead of fill_value.
                                       Changed signal and reference count input from float to fix.
                                       Replaced 'exec' statements.
        2015-04-21: Russell Desiderio. Cosmetic statement re-ordering and shaping of input arrays.
        """
        # test inputs:
        # replicate the inputs to represent 3 data packets
        npackets = 3

        ## natively scalar inputs
        traw = np.array([48355])
        Tcal = np.array([20.0])
        T = np.array([12.0])
        PS = np.array([35.0])
        # time-vectorize
        [traw, Tcal, T, PS] = [np.tile(xarray, npackets) for xarray in [traw, Tcal, T, PS]]

        ## natively 1D inputs
        # valid tscor dictionary entries (those without nan values) are
        # keyed at [400.0 755.0] nm. The dictionary has been extended to
        # to [380.0 775.0] using np.nan as fill values. test:
        cwlngth = np.array([[398.5, 540., 580., 630., 670., 710.]])
        awlngth = np.array([[500., 550., 600., 650., 700., 761.2]])
        c_sig = np.array([[150, 225, 200, 350, 450, 495]])
        c_ref = np.array([[550, 540, 530, 520, 510, 500]])
        c_off = np.array([[1.35, 1.30, 1.25, 1.20, 1.15, 1.10]])
        a_sig = np.array([[250, 300, 210, 430, 470, 495]])
        a_ref = np.array([[450, 460, 470, 480, 490, 500]])
        a_off = np.array([[0.35, 0.30, 0.25, 0.20, 0.15, 0.10]])
        tbins = np.array([[14.5036, 15.5200, 16.4706, 17.4833, 18.4831, 19.5196, 20.5565]])
        # time-vectorize
        [cwlngth, awlngth, c_sig, c_ref, c_off, a_sig, a_ref, a_off, tbins] = [
            np.tile(xarray, (npackets, 1)) for xarray in [
                cwlngth, awlngth, c_sig, c_ref, c_off, a_sig, a_ref, a_off, tbins]]

        # nrows of tarr = length(nwlngth); ncols of tarr = length(tbins)
        tarr = np.array([[
            [0.0, -0.004929, -0.004611, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004611, -0.004418, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004418, -0.004355, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004355, -0.004131, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.004131, -0.003422, 0.0, 0.0, 0.0, 0.0],
            [0.0, -0.003422, -0.002442, 0.0, 0.0, 0.0, 0.0]
        ]])
        # time-vectorize
        # dimensions will be (npackets, nwavelengths, ntempbins)
        tarr = np.tile(tarr, (npackets, 1, 1))

        # expected outputs:
        # 1st c wavelength is out-of-range of the empirical TS corrections, so the DPA
        # should return nan.
        cpd_ts = np.array([[np.nan, 4.807914, 5.156010, 2.788715, 1.655607, 1.171965]])
        # to scatter-correct the 1st abs wavelength value, the c values are interpolated;
        # because the 1st c value is out-of-TScor-range, the 1st abs value will also be
        # a nan value. and, the last abs value will also be a nan value because 761.2>755.0.
        apd_ts_s_ref700 = np.array([[np.nan, 1.021574, 3.235057, 0.099367, 0.000000, np.nan]])

        cpd_ts = np.tile(cpd_ts, (npackets, 1))
        apd_ts_s_ref700 = np.tile(apd_ts_s_ref700, (npackets, 1))

        # beam attenuation (beam c)  test
        c_ts = optfunc.opt_beam_attenuation(c_ref, c_sig, traw, cwlngth, c_off, Tcal,
                                            tbins, tarr, T, PS)
        np.testing.assert_allclose(c_ts, cpd_ts, rtol=1e-6, atol=1e-6)

        # absorption test
        a_ts_s700 = optfunc.opt_optical_absorption(a_ref, a_sig, traw, awlngth, a_off, Tcal,
                                                   tbins, tarr, cpd_ts, cwlngth, T, PS)
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
        Modified by Russell Desiderio, April 09, 2015.
            check function operation without iterating on the input data using for loop.
            check both cal coeff cases: (1) calcoeffs are scalars
                                        (2) calcoeffs are 1d arrays in time.
        Modified by Russell Desiderio, April 21, 2015.
            removed test using 1D data input array with scalar calcoeffs.
            added scalar data with scalar calcoeff test.
        """

    # test inputs
    Im = np.array([1.3589])
    a0 = np.array([2156849800.8])
    a1 = np.array([2.586852835e-006])

    par_count = np.array([2159403328, 2159400384, 2159396992, 2159400384, 2159407488, 2159399296, 2159400384,
                          2159404800, 2159403904, 2159402240, 2159403200, 2159409728, 2159409792, 2159408320,
                          2159407808, 2159402304, 2159402688, 2159407552, 2159404160, 2159403776, 2159402048,
                          2159404544])

    # expected outputs
    par_expected = np.array([8.976348585, 8.965999618, 8.954075807, 8.965999618, 8.990972126, 8.962174999,
                             8.965999618, 8.981523069, 8.978373383, 8.972523967, 8.97589863, 8.998846341,
                             8.999071318, 8.993896835, 8.992097014, 8.972748944, 8.97409881, 8.991197104,
                             8.979273293, 8.977923428, 8.971849034, 8.980623159])

    # compute par value for 'scalar' in time case.
    par_calc = optfunc.opt_par_satlantic(par_count[0], a0, a1, Im)
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected[0], rtol=0.000001, atol=0.000001)

    # compute par values: time-vectorized cal coeff case
    # all inputs and output are 1D arrays of shape (tval,)
    tval = par_count.shape[0]
    a0_tiled = np.tile(a0, tval)
    a1_tiled = np.tile(a1, tval)
    Im_tiled = np.tile(Im, tval)
    par_calc = optfunc.opt_par_satlantic(par_count, a0_tiled, a1_tiled, Im_tiled)
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected, rtol=0.000001, atol=0.000001)

    def test_opt_par_wetlabs(self):
        """
        Test the opt_par_wetlabs function.

        Values based on that described in DPS as available on Alfresco:

        OOI (2014). Data Product Specification for DATA PRODUCT SPECIFICATION
        FOR PHOTOSYNTHETICALLY ACTIVE RADIATION (PAR) FROM WET LABS INSTRUMENT
        ON COASTAL SURFACE PIERCING PROFILER Document Control Number 1341-00722.
        https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00722_Data_Product_SPEC_OPTPARW_WETLabs_OOI.pdf)

        Implemented by Craig Risien, December 2014
        Modified by Russell Desiderio, April 09, 2015.
            check function operation without iterating on the input data using for loop.
            check both cal coeff cases: (1) calcoeffs are scalars
                                        (2) calcoeffs are 1d arrays in time.
        Modified by Russell Desiderio, April 21, 2015.
            removed test using 1D data input array with scalar calcoeffs.
            added scalar data with scalar calcoeff test.
        """

    # test inputs
    Im = np.array([1.3589])
    a0 = np.array([4381])
    a1 = np.array([2904])

    par_count = np.array([4975, 4999, 5077, 5230, 5413, 5483, 5427, 5344, 5285, 5284,
                          5344, 5470, 5649, 5835, 6007, 11824, 4975, 4999, 5077, 5230,
                          5413, 11824])

    # expected outputs
    par_expected = np.array([2.176371115, 2.218183222, 2.359700693, 2.664053061, 3.08006158,
                             3.255847708, 3.114442636, 2.916077547, 2.782801729, 2.780596116,
                             2.916077547, 3.222459729, 3.713869311, 4.304039059, 4.932928879,
                             496.8256708, 2.176371115, 2.218183222, 2.359700693, 2.664053061,
                             3.08006158, 496.8256708])

    # compute par value for 'scalar' in time case.
    par_calc = optfunc.opt_par_wetlabs(par_count[0], a0, a1, Im)
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected[0], rtol=0.000001, atol=0.000001)

    # compute par values: time-vectorized cal coeff case
    tval = par_count.shape[0]
    a0_tiled = np.tile(a0, tval)
    a1_tiled = np.tile(a1, tval)
    Im_tiled = np.tile(Im, tval)
    par_calc = optfunc.opt_par_wetlabs(par_count, a0_tiled, a1_tiled, Im_tiled)
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
        Modified by Russell Desiderio, April 09, 2015.
            check function operation without iterating on the input data using for loop.
            check both cal coeff cases: (1) calcoeffs are scalars
                                        (2) calcoeffs are 1d arrays in time.
        Modified by Russell Desiderio, April 21, 2015.
            removed test using 1D data input array with scalar calcoeffs.
            added scalar data with scalar calcoeff test.
        """

    # test inputs
    offset = np.array([0.0101])
    scale = np.array([5.897e-04])

    par_volts = np.array([1.016793, 0.599800, 0.452400, 0.305000, 0.187000, 0.178900, 0.069100, 0.039600, 0.010100])

    par_expected = np.array([1707.13, 1000.00, 750.00, 500.00, 300.00, 286.25, 100.00, 50.00, 0.00])

    # scalar time case
    par_calc = optfunc.opt_par_biospherical_mobile(par_volts[0], offset, scale)
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected[0], rtol=0.01, atol=0.01)

    # compute par values: time-vectorized cal coeff case
    tval = par_volts.shape[0]
    offset_tiled = np.tile(offset, tval)
    scale_tiled = np.tile(scale, tval)
    par_calc = optfunc.opt_par_biospherical_mobile(par_volts, offset_tiled, scale_tiled)
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
        Modified by Russell Desiderio, April 09, 2015.
            check function operation without iterating on the input data using for loop.
            check both cal coeff cases: (1) calcoeffs are scalars
                                        (2) calcoeffs are 1d arrays in time.
        Modified by Russell Desiderio, April 21, 2015.
            removed test using 1D data input array with scalar calcoeffs.
            added scalar data with scalar calcoeff test.
        """

    # test inputs
    #offset units = mvolts
    offset = np.array([10.1])
    scale = np.array([5.897e-04 / (6.02 * 10**13)])

    par_mvolts = np.array([1016.793, 599.800, 452.400, 305.000, 187.000, 178.900, 69.100, 39.600, 10.100])

    par_expected = np.array([1707.13, 1000.00, 750.00, 500.00, 300.00, 286.25, 100.00, 50.00, 0.00])

    # scalar case
    par_calc = optfunc.opt_par_biospherical_wfp(par_mvolts[0], offset, scale)
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected[0], rtol=0.01, atol=0.01)

    # compute par values: time-vectorized cal coeff case
    tval = par_mvolts.shape[0]
    offset_tiled = np.tile(offset, tval)
    scale_tiled = np.tile(scale, tval)
    par_calc = optfunc.opt_par_biospherical_wfp(par_mvolts, offset_tiled, scale_tiled)
    # compare calculated results to expected
    np.testing.assert_allclose(par_calc, par_expected, rtol=0.01, atol=0.01)

    def test_opt_ocr507_irradiance(self):
        """
        Test opt_ocr507_irradiance function.

        Values calculated using data product algorithm in the DPS in an excel spreadsheet
        as available on Alfresco:

        OOI (2014). Data Product Specification for Downwelling Spectral Irradiance.
            Document Control Number 1341-00730. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00730__???.pdf)

        OOI (2014). SPECTIR Unit Test. 1341-00730_SPECTIR Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00730_SPECTIR >>
            SPKIR_SPECTIR_unit_test.xlsx)


        Implemented by Russell Desiderio, March 14, 2014.
        Modified by Russell Desiderio, March 25, 2014.
            Transposed test array.
            Set known output to be 2D row vector.
        Modified by Russell Desiderio, April 09, 2015.
            Added unit tests for multiple data packets and exception checking
            Set known output for 1 data packet case to be 1D array.
        Modified by Russell Desiderio, April 21, 2015.
            Conditioned inputs for 1 data packet case to be 2D row vectors.
            Set known output for 1 data packet case back to 2D row vector.
        """
#              counts         offset          scale        mrsn     Ed
        test_array = np.transpose(np.array([
            [2148370944, 2148377867.8, 2.09023117662E-07, 1.368, -0.00198],
            [2200000000, 2148218092.4, 2.06543624674E-07, 1.410, 15.080],
            [2300000000, 2147607229.7, 2.12484770952E-07, 1.365, 44.200],
            [2400000000, 2147789959.1, 2.07241106309E-07, 1.354, 70.771],
            [2500000000, 2148047456.7, 1.99358530187E-07, 1.372, 96.266],
            [2600000000, 2147335412.8, 2.06033896796E-07, 1.404, 130.943],
            [2700000000, 2146998228.4, 2.14806273478E-07, 1.347, 160.008]
        ]))

        ## one data packet case
        # set inputs
        counts = test_array[[0], :]
        offset = test_array[[1], :]
        scale = test_array[[2], :]
        immersion_factor = test_array[[3], :]

        # set known output to be a 2D row vector
        Ed = test_array[[-1], :]
        #print Ed.shape

        # calculate the downwelling irradiance
        Ed_out = optfunc.opt_ocr507_irradiance(counts, offset, scale, immersion_factor)

        ###########################################################################
        # The DPS specifies a precision of 0.25 uW/cm^2/nm
        np.testing.assert_allclose(Ed_out, Ed, rtol=0.0, atol=0.1)
        ###########################################################################

        ## multiple data packets case
        ## cal coeffs will also be "time-vectorized" by CI
        # create 10 data packets
        z_value = 10
        zcounts = np.tile(counts, (z_value, 1))
        zoffset = np.tile(offset, (z_value, 1))
        zscale = np.tile(scale, (z_value, 1))
        zimmersion_factor = np.tile(immersion_factor, (z_value, 1))

        # set known output
        zEd = np.tile(Ed, (z_value, 1))
        #print Ed.shape

        # calculate the downwelling irradiance
        zEd_out = optfunc.opt_ocr507_irradiance(zcounts, zoffset, zscale, zimmersion_factor)

        ###########################################################################
        # The DPS specifies a precision of 0.25 uW/cm^2/nm
        np.testing.assert_allclose(zEd_out, zEd, rtol=0.0, atol=0.1)
        ###########################################################################

        ## test error-checking
        # data array counts does not have 7 elements
        counts_wrongshape = test_array[0, 0:-1]
        np.testing.assert_raises(ValueError, optfunc.opt_ocr507_irradiance, counts_wrongshape,
                                 offset, scale, immersion_factor)

        # cal coeff scale does not have 7 elements
        scale_wrongshape = test_array[2, 0:-1]
        np.testing.assert_raises(ValueError, optfunc.opt_ocr507_irradiance, counts, offset,
                                 scale_wrongshape, immersion_factor)

