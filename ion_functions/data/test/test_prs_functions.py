#!/usr/bin/env python

"""
@package ion_functions.test.prs_functions
@file ion_functions/test/test_prs_functions.py
@author Russell Desiderio, Chris Wingard
@brief Unit tests for prs_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
import scipy.io as sio
import datetime as dt
import matplotlib.dates as mdates
from ion_functions.data import prs_functions as prsfunc
from ion_functions.utils import fill_value


@attr('UNIT', group='func')
class TestPRSFunctionsUnit(BaseUnitTestCase):

    # BOTTILT TEST DATA
    # setup up a list with the bottilt test data used by the three test functions
    #      xtilt ,  ytilt , scmp, sernum , ccmp, tmag, angle, tdir
    lily = [
        [-150.000, 110.000, 0.15, 'N9651', 183, 186.0, -36.3, 129],
        [-140.000, 100.000, 71.99, 'N9651', 170, 172.0, -35.5, 116],
        [-130.000, 90.000, 144.21, 'N9651', 92, 158.1, -34.7, 37],
        [-120.000, 80.000, 215.53, 'N9651', 317, 144.2, -33.7, 261],
        [-110.000, 70.000, 288.12, 'N9651', 204, 130.4, -32.5, 146],
        [-100.000, 60.000, 359.76, 'N9651', 183, 116.6, -31.0, 124],
        [-90.000, 50.000, 0.15, 'N9676', 200, 103.0, -29.1, 139],
        [-80.000, 40.000, 71.99, 'N9676', 132, 89.4, -26.6, 69],
        [-70.000, 30.000, 144.21, 'N9676', 53, 76.2, -23.2, 346],
        [0.000, 20.000, 215.53, 'N9676', 332, 20.0, 90.0, 332],
        [0.000, 10.000, 288.12, 'N9676', 263, 10.0, 90.0, 263],
        [0.000, 0.000, 359.76, 'N9676', 200, 0.0, 0.0, 290],
        [0.000, -10.000, 0.15, 'N9656', 166, 10.0, -90.0, 346],
        [0.000, -20.000, 71.99, 'N9656', 130, 20.0, -90.0, 310],
        [10.000, -30.000, 144.21, 'N9656', 87, 31.6, -71.6, 249],
        [20.000, -40.000, 215.53, 'N9656', 351, 44.7, -63.4, 144],
        [30.000, -50.000, 288.12, 'N9656', 240, 58.3, -59.0, 29],
        [40.000, -60.000, 359.76, 'N9656', 166, 72.1, -56.3, 312],
        [50.000, -70.000, 0.15, 'N9652', 173, 86.0, -54.5, 317],
        [60.000, -80.000, 71.99, 'N9652', 133, 100.0, -53.1, 276],
        [70.000, -90.000, 144.21, 'N9652', 82, 114.0, -52.1, 224],
        [80.000, -100.000, 215.53, 'N9652', 347, 128.1, -51.3, 128],
        [90.000, -110.000, 288.12, 'N9652', 243, 142.1, -50.7, 24],
        [100.000, -120.000, 359.76, 'N9652', 173, 156.2, -50.2, 313],
        [110.000, -130.000, 0.15, 'N9655', 173, 170.3, -49.8, 313],
        [120.000, 160.000, 71.99, 'N9655', 152, 200.0, 53.1, 189],
        [130.000, 150.000, 144.21, 'N9655', 98, 198.5, 49.1, 139],
        [140.000, 140.000, 215.53, 'N9655', 337, 198.0, 45.0, 22],
        [150.000, 130.000, 288.12, 'N9655', 222, 198.5, 40.9, 271],
        [160.000, 120.000, 359.76, 'N9655', 173, 200.0, 36.9, 226],
        [-10.000, -120.000, 0.15, 'N9651', 183, 120.4, 85.2, 8],
        [-60.000, -220.000, 71.99, 'N9651', 170, 228.0, 74.7, 5],
        [-110.000, -5.000, 144.21, 'N9651', 92, 110.1, 2.6, 359],
        [-160.000, -150.000, 215.53, 'N9651', 317, 219.3, 43.2, 184],
        [-240.000, -260.000, 288.12, 'N9651', 204, 353.8, 47.3, 67],
        [-310.000, -10.000, 359.76, 'N9651', 183, 310.2, 1.8, 91]
    ]

    # BOTSFLU TEST DATA
    # (1) construct test values for event warnings.
    # generate 3 rows of random numbers in the interval [0.0, 1.0)
    event_values = np.random.random_sample((3, 15))
    # set up 3 rows of data:
    #    row 0: no |values| >= 1.0
    #    row 1: one value >= 1.0
    #    row 2: one value <= =1.0
    # also put in a nan in each row to test robustness
    event_values[:, 4] = np.nan
    neg_elem = np.array([1, 3, 11])
    event_values[:, neg_elem] = -event_values[:, neg_elem]
    event_values[1, 6] = 1.1
    event_values[2, 12] = -1.1

    """
    BOTSFLU unit test data are 20 Hz pressure data reverse engineered from
    15-sec binned depth data from Feb-Apr 2011.

    Implemented by Russell Desiderio, January 2015.

    Originally the unit test data were imported into ionfunctions;
    however, there is a github filesize limit of 100 MB, so that
    the binary testdata files could not be uploaded in github.

    This code is retained but commented out, with the matpath set
    to a directory outside of the virtual environment where the
    large files are stored, in case I need to re-run it. The unit
    test check values were calculated using these data and comparing
    their calculated data products with those calculated by a matlab
    program I also wrote.
    """
    ## (2) import test data as unsigned 4 byte integers in v7 mat files.
    #matpath = '/media/sf_OOI/'
    ## import seconds_since_1970, most significant part [units = ksec]
    #dict_ss1970_ksec = sio.loadmat(matpath + 'botsflu_ss1970_ksec.mat')
    #ss1970_ksec = dict_ss1970_ksec['ss1970_ksec']
    ## import seconds_since_1970, least significant part [units = msec]
    #dict_ss1970_msec = sio.loadmat(matpath + 'botsflu_ss1970_msec.mat')
    #ss1970_msec = dict_ss1970_msec['ss1970_msec']
    ## construct OOI timestamps [sec since 01-01-1900]
    ## to convert between unix [sec since 01-01-1970] and OOI epochs, use
    #delta_epoch = 2208988800.0  # [sec]
    #ss1900 = (1000 * ss1970_ksec.reshape((-1)) +
    #          0.001 * ss1970_msec.reshape((-1))) + delta_epoch
    ## bottom pressure [units = 0.0001 psi]
    #dict_psi_0001 = sio.loadmat(matpath + 'botsflu_psi_0001.mat')
    #botpres = 0.0001 * dict_psi_0001['psi_0001']
    #botpres = botpres.reshape((-1))

    """
    Because the binary data files in the above code could not be uploaded in
    github, the reverse engineering of 20 hz pressure data from 15sec binned
    depth data is done below; the size of the binary depth file is much less
    than 100 MB.

    BEGIN generating BOTSFLU 20 hz unit test pressure data
    """
    # depth [units = 0.00001 m]
    matpath = 'ion_functions/data/matlab_scripts/botpt/'
    dict_depth = sio.loadmat(matpath + 'botsflu_15secbin_depth_00001.mat')
    depth15s = 0.00001 * dict_depth['depth_00001']
    depth15s = depth15s.reshape((-1))

    # first, convert missing pressure values from 0 to nan
    depth15s[depth15s == 0] = np.nan

    # generate 20 Hz data by ramping 300 data points around 0.
    # use a 0.002m maximum variablity.
    ramp = np.linspace(-0.002, 0.002, 300)
    # randomly permute them
    vvar = ramp[np.random.permutation(300)]
    # use the same permutation for each depth value;
    # add each depth as a scalar to its own row vector copy of vvar
    # matlab way of doing this: botpres = bsxfun(@plus, depth15s, vvar)
    botpres = depth15s[:, np.newaxis] + vvar
    # botpres is 2D; flatten into a vector
    botpres = botpres.flatten()  # default is row-major, as desired
    # convert to pressure by inverting the DPS method, and ignoring atm pressure.
    botpres = botpres / -0.67
    # pressure values were formerly imported as unsigned integers in units of
    # 0.0001 psi. if the rounding operation below is omitted, one of the unit
    # tests (10minrate) will fail. True in Jan 2015 and May 2017.
    botpres = np.around(botpres, decimals=4)

    # generate OOI 20 Hz time stamps.
    # hard code first timestamp:
    # to convert between unix [sec since 01-01-1970] and OOI epochs, use
    delta_epoch = 2208988800.0  # [sec]
    starttime = 1296518392.525 + delta_epoch  # seconds since 1900-01-01
    ss1900 = np.arange(botpres.size) * 0.05  # each step is 1/20 sec
    ss1900 = ss1900 + starttime

    #### May 2017
    #### Input Nans are now trapped out by the DPA, so the next 3 executable
    #### statements are commented out. Nevertheless, DO NOT DELETE THEM.
    ## find where the nan values are, and delete in both the
    ## pressure record and in the timestamps.
    #nan_position_mask = np.isnan(botpres)
    #botpres = botpres[~nan_position_mask]
    #ss1900 = ss1900[~nan_position_mask]

    """
    END generating 20 hz unit test data for BOTSFLU
    """

    # indices to check for BOTSFLU data products associated with 15sec timestamps
    size_15s_minus_1 = 501119.0
    idx_15s = np.around(np.linspace(0.0, size_15s_minus_1, 20)).astype(int)
    # indices to check for data products associated with 24hr timestamps
    size_24h_minus_1 = 89.0
    idx_24h = np.around(np.linspace(0.0, size_24h_minus_1, 16)).astype(int)

    ### ########################## ###################################################
    ### BOTTILT DATA PRODUCT TESTS ###################################################
    ### ########################## ###################################################
    def test_prs_bottilt_ccmp(self):
        """
        Test prs_bottilt_ccmp function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)

        Note, DPS does not specify a test for this function. Values from CCMP
        lookup table in DPS used to create this test by the function and
        test_function author.

        Implemented by Christopher Wingard, July 2013
        """
        # set known inputs
        scmp = np.atleast_1d([row[2] for row in self.lily])
        snum = np.atleast_1d([row[3] for row in self.lily])

        # set known output for the corrected compass direction
        ccmp = np.atleast_1d([row[4] for row in self.lily])

        # calculate the corrected compass direction
        out = prsfunc.prs_bottilt_ccmp(scmp, snum)

        # How'd we do?
        np.testing.assert_array_equal(out, ccmp)

    def test_prs_bottilt_tmag(self):
        """
        Test prs_bottilt_tmag function.

        Note, DPS specifies a test data set that does not vary the X-tilt and
        Y-tilt inputs to this function. They are 330 and -330 microradians,
        respectively, for the entirety of the dataset. Test values below were
        created by function and test_function author to ensure function
        operates as expected over a range of potential values.

        OOI (2012). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)

        Implemented by Christopher Wingard, July 2013
        """
        # set inputs
        xtilt = np.atleast_1d([row[0] for row in self.lily])
        ytilt = np.atleast_1d([row[1] for row in self.lily])

        # set known output for the tilt magnitude
        tmag = np.atleast_1d([row[5] for row in self.lily])

        # calculate the tilt magnitude
        out = prsfunc.prs_bottilt_tmag(xtilt, ytilt)

        # How'd we do?
        np.testing.assert_allclose(out, tmag, rtol=0.1, atol=0.1)

    def test_prs_bottilt_tdir(self):
        """
        Test prs_bottilt_tdir function.

        Note, DPS specifies a test data set that does not vary the X-tilt and
        Y-tilt inputs to this function. They are 330 and -330 microradians,
        respectively, for the entirety of the dataset. Test values below were
        created by function and test_function author to ensure function
        operates as expected over a range of potential values.

        OOI (2012). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)

        Implemented by:

            2013-07: Christopher Wingard. Initial code.
            2014-03-20: Russell Desiderio. Added 3rd quadrant (xtilt, ytilt) test values
                        to check atan2 implementation.

        """
        # set inputs
        xtilt = np.atleast_1d([row[0] for row in self.lily])
        ytilt = np.atleast_1d([row[1] for row in self.lily])
        ccmp = np.atleast_1d([row[4] for row in self.lily])

        # set known output
        tdir = np.atleast_1d([row[7] for row in self.lily])

        # calculate the tilt direction
        out = prsfunc.prs_bottilt_tdir(xtilt, ytilt, ccmp)

        # How'd we do?
        np.testing.assert_array_equal(out, tdir)

    ### ########################## ###################################################
    ### BOTSFLU DATA PRODUCT TESTS ###################################################
    ### ########################## ###################################################
    def test_anchor_bin_detided_data_to_24h(self):
        """
        Test the function anchor_bin_detided_data_to_24h.
            (1) the 24hr timestamps now include values for empty bins (bin values are nan).
            (2) the 90% threshold coverage criterion is tested.

        Even though there are only 2 24-hour gaps of missing data, 4 consecutive 24-hr bins
        will have no data because (a) there is a 90% coverage threshold and (b) the missing
        data are not aligned with the boundaries of the midnight-centered bins (which would
        be noons).

        The pattern of the data gaps is taken directly from actual data (February 2011,
        Axial Seamount) supplied by RSN.

        Russell Desiderio, May 15, 2017. Original code.
        """
        # generate 8 days of 15sec timestamps with the following gaps:
        #    24 hrs of missing data from day 3 at 03:00:00 to day 4 at 02:59:45
        #    24 hrs of missing data from day 4 at 21:00:00 to day 5 at 20:59:45
        time15s = np.hstack((1 * 86400.0 + np.arange(0, 86386, 15),
                             2 * 86400.0 + np.arange(0, 86386, 15),
                             3 * 86400.0 + np.arange(0, 10786, 15),
                             4 * 86400.0 + np.arange(10800, 75586, 15),
                             5 * 86400.0 + np.arange(75600, 86386, 15),
                             6 * 86400.0 + np.arange(0, 86386, 15),
                             7 * 86400.0 + np.arange(0, 86386, 15),
                             8 * 86400.0 + np.arange(0, 86386, 15)))
        data = time15s * 0 + 1  # therefore all the non-Nan bin values will be 1.0
        dday_coverage = 0.90
        # expected values
        xpctd_timestamps = np.arange(1.0, 10.0) * 86400  # 9 days centered at midnight
        nan = np.nan
        xpctd_data = np.array([nan, 1.0, nan, nan, nan, nan, 1.0, 1.0, nan])
        # there are 86400/15 = 5760 15-sec bins in 1 day
        xpctd_bincount = np.array([2880, 5760, 3600, 2160, 2160, 3600, 5760, 5760, 2880])
        # calc is a tuple of three elements
        calc = prsfunc.anchor_bin_detided_data_to_24h(time15s, data, dday_coverage)
        # test
        np.testing.assert_array_equal(calc[0], xpctd_timestamps)
        np.testing.assert_array_equal(calc[1], xpctd_data)
        np.testing.assert_array_equal(calc[2], xpctd_bincount)

        # also test potentially pathological case mimicked with an unphysical coverage:
        dday_coverage = 1.1  # 110%
        xpctd_timestamps = np.arange(1.0, 10.0) * 86400  # same as above
        xpctd_data = np.zeros(9.0) + np.nan
        xpctd_bincount = np.array([2880, 5760, 3600, 2160, 2160, 3600, 5760, 5760, 2880])  # as above
        # calc is a tuple of three elements
        calc = prsfunc.anchor_bin_detided_data_to_24h(time15s, data, dday_coverage)
        # test
        np.testing.assert_array_equal(calc[0], xpctd_timestamps)
        np.testing.assert_array_equal(calc[1], xpctd_data)
        np.testing.assert_array_equal(calc[2], xpctd_bincount)

    def test_calculate_sliding_means(self):
        """
        Test the calculation of window-centered rolling means.

        Matlab and numpy have different conventions for convolutions involving
        even window sizes; specifically, make sure edge effects are nan'd out
        appropriately.

        Russell Desiderio, January 2015. Original code.
        May 2017: Added case of odd window-size for use in sliding_slope algorithms
                  implementing threshold coverage criteria.
        """
        xpctd = np.array([np.nan, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, np.nan, np.nan])
        data = np.arange(10.0)
        window_size = 4
        calc = prsfunc.calculate_sliding_means(data, window_size)
        np.testing.assert_array_equal(calc, xpctd)

        xpctd = np.array([np.nan, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, np.nan])
        data = np.arange(10.0)
        window_size = 3
        calc = prsfunc.calculate_sliding_means(data, window_size)
        np.testing.assert_array_equal(calc, xpctd)

    def test_calculate_sliding_slopes(self):
        """
        Test the calculation of backwards-looking rolling slopes.

        Verify the placement of nans in the output array.

        Russell Desiderio, January 2015. Original code.
        Russell Desiderio, May 2017: updated function call by adding rate_coverage argument
        """
        # xpctd and calc are not equal, presumably due to roundoff error
        xpctd = np.array([np.nan, np.nan, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        data = np.arange(10.0)
        window_size = 3
        rate_coverage = 1.0
        calc = prsfunc.calculate_sliding_slopes(data, window_size, rate_coverage)
        np.testing.assert_allclose(calc, xpctd, rtol=0.0, atol=1.e-12)

        # Verify nan placement after routine has been altered to accept a 75% fill criterion.
        # set up test data with slope=2
        data = 2 * np.arange(150.0)
        # nan out middle section
        data[50:100] = np.nan

        # when the sliding slope routine is applied, the result will be sets of contiguous nans
        # and contiguous non-Nan values whose indices depend on the window size and on the
        # location of nans in the test data.
        #    29 point window (4weekrate case)
        xpctd_w29 = np.ones(150) + 1.0
        xpctd_w29[0:21] = np.nan
        xpctd_w29[57:121] = np.nan
        window_size = 29
        rate_coverage = 0.75
        calc = prsfunc.calculate_sliding_slopes(data, window_size, rate_coverage)
        np.testing.assert_array_equal(calc, xpctd_w29)
        #    57 point window (8weekrate case)
        xpctd_w57 = np.ones(150) + 1.0
        xpctd_w57[0:42] = np.nan
        xpctd_w57[64:142] = np.nan
        window_size = 57
        rate_coverage = 0.75
        calc = prsfunc.calculate_sliding_slopes(data, window_size, rate_coverage)
        np.testing.assert_array_equal(calc, xpctd_w57)

        # make sure that routine doesn't crash if no windows satisfy the
        # rate_coverage criterion
        rate_coverage = 1.1  # 110%; that ought to do it
        xpctd = np.ones(150) + np.nan
        calc = prsfunc.calculate_sliding_slopes(data, window_size, rate_coverage)
        np.testing.assert_array_equal(calc, xpctd)

        # one last potentially pathological case
        rate_coverage = 0.0
        xpctd = np.ones(150) + np.nan
        allnans = np.copy(xpctd)
        calc = prsfunc.calculate_sliding_slopes(allnans, window_size, rate_coverage)
        np.testing.assert_array_equal(calc, xpctd)

    def test_prs_botsflu_auxiliary_timestamps(self):
        """
        Test the calculation of the auxiliary timestamp data products
        TIME15S and TIME24H.

        Uses matplotlib.dates as mdates.

        Implemented by Russell Desiderio, January 2015.
        No change after May 2017 modifications.
        """
        xpctd_time15s = ['2011 02 01 00 00 00',
                         '2011 02 05 13 53 45',
                         '2011 02 10 03 47 15',
                         '2011 02 16 17 41 00',
                         '2011 02 21 07 34 45',
                         '2011 02 25 21 28 15',
                         '2011 03 02 11 22 00',
                         '2011 03 07 01 15 45',
                         '2011 03 11 15 09 15',
                         '2011 03 16 05 03 00',
                         '2011 03 20 18 56 45',
                         '2011 03 25 08 50 30',
                         '2011 03 29 22 44 00',
                         '2011 04 03 12 37 45',
                         '2011 04 08 02 31 30',
                         '2011 04 12 16 25 00',
                         '2011 04 17 06 18 45',
                         '2011 04 21 20 12 30',
                         '2011 04 26 10 06 00',
                         '2011 04 30 23 59 45']

        reftime = mdates.date2num(dt.datetime(1900, 1, 1))
        time15s = prsfunc.prs_botsflu_time15s(self.ss1900, self.botpres)
        # convert to datetime objects
        t15s_dto = np.array(mdates.num2date(time15s/86400.0+reftime))
        # convert 20 values to readable date-times
        calc = [dto.strftime('%Y %m %d %H %M %S') for dto in t15s_dto[self.idx_15s]]
        np.testing.assert_equal(calc, xpctd_time15s)

        xpctd_time24h = ['2011 02 01 00 00 00',
                         '2011 02 07 00 00 00',
                         '2011 02 13 00 00 00',
                         '2011 02 19 00 00 00',
                         '2011 02 25 00 00 00',
                         '2011 03 03 00 00 00',
                         '2011 03 09 00 00 00',
                         '2011 03 15 00 00 00',
                         '2011 03 20 00 00 00',
                         '2011 03 26 00 00 00',
                         '2011 04 01 00 00 00',
                         '2011 04 07 00 00 00',
                         '2011 04 13 00 00 00',
                         '2011 04 19 00 00 00',
                         '2011 04 25 00 00 00',
                         '2011 05 01 00 00 00']

        time24h = prsfunc.prs_botsflu_time24h(time15s)
        # convert to datetime objects
        t24h_dto = np.array(mdates.num2date(time24h/86400.0+reftime))
        # convert 16 values to readable date-times
        calc = [dto.strftime('%Y %m %d %H %M %S') for dto in t24h_dto[self.idx_24h]]
        np.testing.assert_equal(calc, xpctd_time24h)

    def test_prs_botsflu_meanpres(self):
        """
        Test the calculation of MEANPRES.

        Implemented by Russell Desiderio, January 2015.
        """
        xpctd_meanpres = np.array([2253.2836, 2254.6806, 2254.7552, 2256.6851, 2255.5060,
                                   2253.4896, 2254.6896, 2254.4776, 2255.0015, 2255.7582,
                                   2255.7612, 2255.4791, 2253.7328, 2254.0881, 2257.8000,
                                   2259.0492, 2260.3119, 2257.8776, 2258.0000, 2257.5761])

        meanpres = prsfunc.prs_botsflu_meanpres(self.ss1900, self.botpres)
        calc = meanpres[self.idx_15s]
        np.testing.assert_allclose(calc, xpctd_meanpres, rtol=0, atol=1.e-4)

    def test_prs_botsflu_predtide(self):
        """
        Test the look-up of PREDTIDE values in the 2014-2019 matfile binary tide table:
        ('ion_functions/data/prs_functions_tides_2014_thru_2019.mat').

        The unit test data cover Feb 2011 to Apr 2011, and the tides for this time
        period are in the binary file 'ion_functions/data/matlab_scripts/botpt/
        tides_15sec_2011_for_unit_tests.mat'. These values are effectively tested
        in the unit tests for MEANPRES and the data products downstream from it.

        Implemented by Russell Desiderio, January 2015.

        Matlab code to calculate tides using TPXO7.2 global model:
            http://polaris.esr.org/ptm_index.html
        Further documentation for the TPXO7.2 global tide model:
            http://volkov.oce.orst.edu/tides/global.html
        """
        test_times = [(2014, 1, 1, 0, 0, 0),
                      (2014, 8, 30, 11, 23, 45),
                      (2015, 12, 31, 10, 0, 0),
                      (2016, 2, 28, 23, 59, 45),
                      (2016, 2, 29, 0, 0, 0),
                      (2016, 2, 29, 12, 0, 0),
                      (2016, 11, 9, 18, 44, 15),
                      (2017, 4, 2, 15, 15, 15),
                      (2018, 10, 17, 4, 23, 30),
                      (2018, 6, 12, 12, 39, 45),
                      (2019, 5, 23, 17, 9, 30),
                      (2020, 1, 1, 0, 0, 0)]

        # the tide values were obtained by running a tide calculation program written in
        # matlab by converting the above date/times to matlab's serial date numbers.
        xpctd_tides = [-1.245763, 0.628805, 0.165821, 0.467970,
                       0.468104, 0.822229, 0.245643, -0.024056,
                       0.041149, -1.458176, -1.262106, 0.665962]

        # construct ooi timestamps
        reftime = mdates.date2num(dt.datetime(1900, 1, 1))
        times = mdates.date2num([dt.datetime(*q) for q in test_times])
        ooi_timestamps = (times-reftime)*86400
        predtide = prsfunc.prs_botsflu_predtide(ooi_timestamps)
        np.testing.assert_allclose(predtide, xpctd_tides, rtol=0, atol=1.e-6)

    def test_prs_botsflu_meandepth(self):
        """
        Test the calculation of MEANDEPTH.

        Implemented by Russell Desiderio, January 2015.
        """
        xpctd_meandepth = np.array([-1510.8648, -1510.9475, -1510.8527, -1510.8415, -1510.8779,
                                    -1510.8516, -1510.8708, -1510.8519, -1510.8164, -1510.8517,
                                    -1510.8154, -1510.8588, -1510.8562, -1510.9655, -1512.7381,
                                    -1513.2261, -1513.1438, -1513.1528, -1513.0736, -1513.0669])

        meandepth = prsfunc.prs_botsflu_meandepth(self.ss1900, self.botpres)
        calc = meandepth[self.idx_15s]
        np.testing.assert_allclose(calc, xpctd_meandepth, rtol=0, atol=1.e-4)

    def test_prs_botsflu_5minrate(self):
        """
        Test the calculation of 5MINRATE.

        Implemented by Russell Desiderio, January 2015.
        """
        xpctd_5minrate = np.array([np.nan,  -0.27994,  0.00367, -0.12633, -0.03283,
                                   -0.00276, -0.02017, -0.02958,  0.32096, -0.03186,
                                   0.05080,  0.01794,  0.04610, -0.00838, -0.06690,
                                   -0.02708, -0.06918,  0.03318, -0.01891,  0.04748])

        b_5minrate = prsfunc.prs_botsflu_5minrate(self.ss1900, self.botpres)
        calc = b_5minrate[self.idx_15s]
        np.testing.assert_allclose(calc, xpctd_5minrate, rtol=0, atol=1.e-4)

    def test_prs_botsflu_10minrate(self):
        """
        Test the calculation of 10MINRATE.

        Implemented by Russell Desiderio, January 2015.
        """
        xpctd_10minrate = np.array([np.nan,  -1.96881,  0.62240,  0.35972,  1.75263,
                                    -0.13044,  0.69612, -2.05217, 25.82545,  1.04359,
                                    3.08025,  0.57874,  1.32764, -1.17460, -1.70137,
                                    -1.36129,  0.04627,  1.13735, -0.30540,  np.nan])

        b_10minrate = prsfunc.prs_botsflu_10minrate(self.ss1900, self.botpres)
        calc = b_10minrate[self.idx_15s]
        np.testing.assert_allclose(calc, xpctd_10minrate, rtol=0, atol=1.e-4)

    def test_prs_botsflu_daydepth(self):
        """
        Test the calculation of DAYDEPTH.

        Russell Desiderio, January 2015. Original code.
        Russell Desiderio, May 12, 2017. Added 0.9 coverage and default coverage test.
        """
        # case(1) January 2015 unit test values
        xpctd_daydepth = np.array([-1510.8705, -1510.9055, -1510.9119, -1510.8491,
                                   -1510.8596, -1510.8354, -1510.8791, -1510.8807,
                                   -1510.8378, -1510.8279, -1510.8530, -1512.2859,
                                   -1513.2018, -1513.1660, -1513.1128, -1513.0478])
        dday_coverage = 0.0  # if a bin contains at least one non-Nan value, bin value!=nan
        daydepth = prsfunc.prs_botsflu_daydepth(self.ss1900, self.botpres, dday_coverage)
        calc = daydepth[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_daydepth, rtol=0, atol=1.e-4)

         # case(2) May 2017 unit test values, coverage = 0.9
        xpctd_daydepth = np.array([np.nan, -1510.9055, np.nan, -1510.8491,
                                   -1510.8596, -1510.8354, -1510.8791, -1510.8807,
                                   -1510.8378, -1510.8279, -1510.8530, -1512.2859,
                                   -1513.2018, -1513.1660, -1513.1128, np.nan])
        dday_coverage = 0.9  # binvalue is nan unless >= 90% of the bin values are good.
        daydepth = prsfunc.prs_botsflu_daydepth(self.ss1900, self.botpres, dday_coverage)
        calc = daydepth[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_daydepth, rtol=0, atol=1.e-4)

         # case(3) May 2017 default, no coverage specified in argument list
         #     default coverage value is 0.90, so expected values are same as for case (2)
        xpctd_daydepth = np.array([np.nan, -1510.9055, np.nan, -1510.8491,
                                   -1510.8596, -1510.8354, -1510.8791, -1510.8807,
                                   -1510.8378, -1510.8279, -1510.8530, -1512.2859,
                                   -1513.2018, -1513.1660, -1513.1128, np.nan])
        daydepth = prsfunc.prs_botsflu_daydepth(self.ss1900, self.botpres)
        calc = daydepth[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_daydepth, rtol=0, atol=1.e-4)

    def test_prs_botsflu_4wkrate(self):
        """
        Test the calculation of 4WKRATE.

        Implemented by Russell Desiderio, January 2015.
        RDesiderio. May 12, 2015. 3 additional tests for coverage parameters; 4 total
            (1) dday_coverage =  0.0;    rate_coverage = 1.0  (January 2015 unit test values)
            (2) dday_coverage =  0.0;    rate_coverage = 0.75 (intermediate case)
            (3) dday_coverage =  0.9;    rate_coverage = 0.75
            (4) default (no argument)    default (no argument)  (should be same as (3))
        """
        # case (1)
        dday_coverage = 0.0
        rate_coverage = 1.0
        xpctd_4wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  np.nan, 70.5000, 28.8000, -23.6372,
                                  7.8316, 38.7091, 57.6564, -407.2349,
                                  -3024.1664, -4291.1577, -3982.1142, -2118.0963])
        b_4wkrate = prsfunc.prs_botsflu_4wkrate(
            self.ss1900, self.botpres, dday_coverage, rate_coverage)
        calc = b_4wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_4wkrate, rtol=0, atol=1.e-4)

        # case (2)
        dday_coverage = 0.0
        rate_coverage = 0.75
        xpctd_4wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  64.1039, 70.5000, 28.8000, -23.6372,
                                  7.8316, 38.7091, 57.6564, -407.2349,
                                  -3024.1664, -4291.1577, -3982.1142, -2118.0963])
        b_4wkrate = prsfunc.prs_botsflu_4wkrate(
            self.ss1900, self.botpres, dday_coverage, rate_coverage)
        calc = b_4wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_4wkrate, rtol=0, atol=1.e-4)

        # case (3)
        dday_coverage = 0.9
        rate_coverage = 0.75
        xpctd_4wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  np.nan, 63.2291, 0.4971, -23.6372,
                                  7.8316, 38.7091, 57.6564, -407.2349,
                                  -3024.1664, -4291.1577, -3982.1142, -2276.7628])
        b_4wkrate = prsfunc.prs_botsflu_4wkrate(
            self.ss1900, self.botpres, dday_coverage, rate_coverage)
        calc = b_4wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_4wkrate, rtol=0, atol=1.e-4)

        # case (4)
        xpctd_4wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  np.nan, 63.2291, 0.4971, -23.6372,
                                  7.8316, 38.7091, 57.6564, -407.2349,
                                  -3024.1664, -4291.1577, -3982.1142, -2276.7628])
        b_4wkrate = prsfunc.prs_botsflu_4wkrate(self.ss1900, self.botpres)
        calc = b_4wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_4wkrate, rtol=0, atol=1.e-4)

    def test_prs_botsflu_8wkrate(self):
        """
        Test the calculation of 8WKRATE.

        Implemented by Russell Desiderio, January 2015.
        RDesiderio. May 12, 2015. 3 additional tests for coverage parameters; 4 total
            (1) dday_coverage =  0.0;    rate_coverage = 1.0  (January 2015 unit test values)
            (2) dday_coverage =  0.0;    rate_coverage = 0.75 (intermediate case)
            (3) dday_coverage =  0.9;    rate_coverage = 0.75
            (4) default (no argument)    default (no argument)  (should be same as (3))
        """
        # case (1)
        dday_coverage = 0.0
        rate_coverage = 1.0
        xpctd_8wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, 35.7686, -90.1254,
                                  -882.7773, -1506.2385, -1921.8774, -2120.7657])
        b_8wkrate = prsfunc.prs_botsflu_8wkrate(
            self.ss1900, self.botpres, dday_coverage, rate_coverage)
        calc = b_8wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_8wkrate, rtol=0, atol=1.e-4)

        # case (2)
        dday_coverage = 0.0
        rate_coverage = 0.75
        xpctd_8wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, np.nan, 24.5993,
                                  28.7563, 32.6131, 35.7686, -90.1254,
                                  -882.7773, -1506.2385, -1921.8774, -2120.7657])
        b_8wkrate = prsfunc.prs_botsflu_8wkrate(
            self.ss1900, self.botpres, dday_coverage, rate_coverage)
        calc = b_8wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_8wkrate, rtol=0, atol=1.e-4)

        # case (3)
        dday_coverage = 0.9
        rate_coverage = 0.75
        xpctd_8wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, np.nan, np.nan,
                                  26.1114, 30.4322, 31.6141, -113.8871,
                                  -882.7773, -1506.2385, -1921.8774, -2150.3933])
        b_8wkrate = prsfunc.prs_botsflu_8wkrate(
            self.ss1900, self.botpres, dday_coverage, rate_coverage)
        calc = b_8wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_8wkrate, rtol=0, atol=1.e-4)

        # case (4)
        xpctd_8wkrate = np.array([np.nan, np.nan, np.nan, np.nan,
                                  np.nan, np.nan, np.nan, np.nan,
                                  26.1114, 30.4322, 31.6141, -113.8871,
                                  -882.7773, -1506.2385, -1921.8774, -2150.3933])
        b_8wkrate = prsfunc.prs_botsflu_8wkrate(self.ss1900, self.botpres)
        calc = b_8wkrate[self.idx_24h]
        np.testing.assert_allclose(calc, xpctd_8wkrate, rtol=0, atol=1.e-4)

    def test_prs_tsunami_detection(self):
        """
        Test prs_tsunami_detection.

        Implemented by Russell Desiderio, January 2015.
        """
        # threshold is 1.0 cm/min
        data = self.event_values
        xpctd = np.array([False, True, True])

        tf = np.zeros(3) * np.nan
        for ii in range(3):
            tf[ii] = prsfunc.prs_tsunami_detection(data[ii, :])
        np.testing.assert_array_equal(tf, xpctd)

    def test_prs_eruption_imminent(self):
        """
        Test prs_eruption_imminent.

        Implemented by Russell Desiderio, January 2015.
        """
        # threshold is 5.0 cm/hr
        data = self.event_values * 5.0
        xpctd = np.array([False, True, False])

        tf = np.zeros(3) * np.nan
        for ii in range(3):
            tf[ii] = prsfunc.prs_eruption_imminent(data[ii, :])
        np.testing.assert_array_equal(tf, xpctd)

    def test_prs_eruption_occurred(self):
        """
        Test prs_eruption_occurred.

        Implemented by Russell Desiderio, January 2015.
        """
        # threshold is -5.0 cm/hr
        data = self.event_values * 5.0
        xpctd = np.array([False, False, True])

        tf = np.zeros(3) * np.nan
        for ii in range(3):
            tf[ii] = prsfunc.prs_eruption_occurred(data[ii, :])
        np.testing.assert_array_equal(tf, xpctd)
