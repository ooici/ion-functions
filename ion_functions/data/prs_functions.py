#!/usr/bin/env python
"""
@package ion_functions.data.prs_functions
@file ion_functions/data/prs_functions.py
@author Russell Desiderio, Chris Wingard
@brief Module containing calculations related to instruments in the Seafloor
    Pressure family.
"""
import pkg_resources
import numexpr as ne
import numpy as np
import scipy.io as spio
from scipy import signal


"""
    Listing of functions, in order encountered.

    Functions calculating data products.

      BOTTILT:

        prs_bottilt_ccmp -- computes the BOTTILT-CCMP_L1 data product
        prs_bottilt_tmag -- computes the BOTTILT-TMAG_L1 data product
        prs_bottilt_tdir -- computes the BOTTILT-TDIR_L1 data product

      BOTSFLU:

        prs_botsflu_time15s -- computes the TIME15S-AUX auxiliary data product
        prs_botsflu_meanpres -- computes the BOTSFLU-MEANPRES_L2 data product
        prs_botsflu_predtide -- computes the BOTSFLU-PREDTIDE_L2 data product
        prs_botsflu_meandepth -- computes the BOTSFLU-MEANDEPTH_L2 data product
        prs_botsflu_5minrate -- computes the BOTSFLU-5MINRATE_L2 data product
        prs_botsflu_10minrate -- computes the BOTSFLU-10MINRATE_L2 data product
        prs_botsflu_time24h -- computes the TIME24H-AUX auxiliary data product
        prs_botsflu_daydepth -- computes the BOTSFLU-DAYDEPTH_L2 data product
        prs_botsflu_4wkrate -- computes the BOTSFLU-4WKRATE_L2 data product
        prs_botsflu_8wkrate -- computes the BOTSFLU-8WKRATE_L2 data product

    Worker functions called by functions calculating data products.

      BOTSFLU:

        anchor_bin
        calc_meandepth_plus
        calculate_sliding_means
        calculate_sliding_slopes

    Functions calculating event notifications; they return either True or False.

      BOTSFLU:

        prs_tsunami_detection -- event notification specified by DPS
        prs_eruption_imminent -- event notification specified by DPS
        prs_eruption_occurred -- event notification specified by DPS
        calculate_sliding_slopes
"""


def prs_bottilt_ccmp(scmp, sn):
    """
    Description:

        OOI Level 1 Seafloor High-Resolution tilt (BOTTILT) core data product,
        derived from data output by the Applied Geomechanics LILY tilt sensor
        on board the Bottom Pressure Tilt (BOTPT) instruments on the Regional
        Scale Nodes (RSN) at Axial Seamount. This function computes
        BOTTILT-CCMP_L1.

    Implemented by:

        2013-06-10: Christopher Wingard. Initial code.
        2014-03-20: Russell Desiderio. Alternate code: faster, but less direct.

    Usage:

        ccmp = prs_bottilt_ccmp(scmp, sn)

            where

        ccmp = Corrected compass direction (BOTTILT-CCMP_L1) [degrees]
        scmp = Uncorrected sensor compass direction (BOTTILT-SCMP_L0) [degrees]
        sn = LILY sensor serial number [unitless]

    References:

        OOI (2013). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)
    """

    """
        Currently, there are two coded algorithms:
            (1) the straightforward original, which uses a two-element keyed dictionary;
            (2) a faster version, which uses serial number keys to the dictionary.

        Since each algorithm uses its own dictionary, the corresponding import statements
            are TEMPORARILY placed at the beginning of their respective code sections
            instead of at module top.
    """
    ###  Original coding, using a dictionary constructed with 2-element keys.

    # load the corrected compass directions table [(sn, scmp) keys]
    from ion_functions.data.prs_functions_ccmp import cmp_lookup

    # use the lookup table to get the ccmp
    ccmp = np.zeros(len(scmp))

    for i in range(len(scmp)):
        ccmp[i] = cmp_lookup[(sn[i], int(round(scmp[i])))]
    return ccmp


    ####  Faster coding, using a dictionary constructed with 1-element keys.
    #
    ## load the corrected compass directions table [sn keys]
    #from ion_functions.data.prs_functions_ccmp_lily_compass_cals import cmp_cal
    #
    ## initialize output array for vectorized masking operations. this will 'break'
    ##    the code if an invalid serial number is specified in the argument list.
    #ccmp = np.zeros(len(scmp)) + np.nan
    #
    ## round the uncorrected compass values to the nearest integer as specified in the DPS,
    ##    which uses a lookup table consisting of integral values to do the correction.
    #scmp = np.round(scmp)
    #
    ## find the supported tilt sensor serial numbers, which are keys in the dictionary
    #sernum = cmp_cal.keys()
    #
    #for ii in range(len(sernum)):
    #    # get the cal coeffs as a function of the iterated serial number;
    #    #    x is the raw, uncorrected reading (scmp)
    #    #    y is the corrected reading (ccmp)
    #    [x, y] = cmp_cal[sernum[ii]]
    #
    #    # the boolean mask has 'true' entries where the elements of input vector sn
    #    #    agree with the iterated serial number.
    #    # np.core.defchararray.equal handles vector string comparisons.
    #    mask = np.core.defchararray.equal(sn, sernum[ii])
    #
    #    ## np.interp is used to do the 'lookup' for performance reasons (vectorized)
    #    ccmp[mask] = np.interp(scmp[mask], x, y)
    #
    ## round to make sure we get an integral value (but not int type)
    #return np.round(ccmp)


def prs_bottilt_tmag(x_tilt, y_tilt):
    """
    Description:

        OOI Level 1 Seafloor High-Resolution Tilt (BOTTILT) core data product,
        derived from data output by the Applied Geomechanics LILY tilt sensor
        on board the Bottom Pressure Tilt (BOTPT) instruments on the Regional
        Scale Nodes (RSN) at Axial Seamount. This function computes
        BOTTILT-TMAG_L1.

    Implemented by:

        2013-06-10: Christopher Wingard. Initial code.

    Usage:

        tmag = prs_bottilt(x_tilt, y_tilt)

            where

        tmag = Resultant tilt magnitude (BOTTILT-TMAG_L1) [microradians]
        x_tilt = Sensor X_tilt (BOTTILT-XTLT_L0) [microradians]
        y_tilt = Sensor Y_tilt (BOTTILT-YTLT_L0) [microradians]

    References:

        OOI (2013). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)
     """
    tmag = ne.evaluate('sqrt(x_tilt**2 + y_tilt**2)')
    return tmag


def prs_bottilt_tdir(x_tilt, y_tilt, ccmp):
    """
    Description:

        OOI Level 1 Seafloor High-Resolution Tilt (BOTTILT) core data product,
        derived from data output by the Applied Geomechanics LILY tilt sensor
        on board the Bottom Pressure Tilt (BOTPT) instruments on the Regional
        Scale Nodes (RSN) at Axial Seamount. This function computes
        BOTTILT-TDIR_L1.

    Implemented by:

        2013-06-10: Christopher Wingard. Initial code.
        2014-03-20: Russell Desiderio. Replaced initial code with arctan2 implementation.

    Usage:

        tdir = prs_bottilt(x_tilt, y_tilt, ccmp)

            where

        tdir = Resultant tilt direction (BOTTILT-TDIR_L1) [degrees]
        x_tilt = Sensor X_tilt (BOTTILT-XTLT_L0) [microradians]
        y_tilt = Sensor Y_tilt (BOTTILT-YTLT_L0) [microradians]
        ccmp = Corrected compass direction (BOTTILT-CCMP_L1) [degrees]

    References:

        OOI (2013). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-00060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00060_Data_Product_SPEC_BOTTILT_OOI.pdf)
     """
    ### As originally coded, according to the algorithm specified in the DPS:

    ## Calculate the angle to use in the tilt direction formula
    ## default angle calculation -- in degrees
    #angle = ne.evaluate('arctan(y_tilt / x_tilt)')
    #angle = np.degrees(angle)
    #
    ## if X-Tilt == 0 and Y-Tilt > 0
    #mask = np.logical_and(x_tilt == 0, y_tilt > 0)
    #angle[mask] = 90.0
    #
    ## if X-Tilt == 0 and Y-Tilt < 0
    #mask = np.logical_and(x_tilt == 0, y_tilt < 0)
    #angle[mask] = -90.0
    #
    ## if Y-Tilt == 0
    #mask = np.equal(y_tilt, np.zeros(len(y_tilt)))
    #angle[mask] = 0.0
    #
    ### Calculate the tilt direction, using the X-Tilt to set the equation
    ## default tilt direction equation
    #tdir = ne.evaluate('(270 - angle + ccmp) % 360')
    #
    ## if X-Tilt >= 0
    #tmp = ne.evaluate('(90 - angle + ccmp) % 360')
    #mask = np.greater_equal(x_tilt, np.zeros(len(x_tilt)))
    #tdir[mask] = tmp[mask]
    #
    #return np.round(tdir)

    # The preceding calculation is faster and simpler if the arctan2 function is used.
    # Use 450 as an addend in the first argument to the mod function to make sure the result is positive.
    return np.round(np.mod(450 - np.degrees(np.arctan2(y_tilt, x_tilt)) + ccmp, 360))


#**********************************************************************
#.. BOTSFLU: Functions calculating data products
#**********************************************************************
def prs_botsflu_time15s(timestamp, botpres):
    """
    Description:

        Calculates the auxiliary BOTSFLU data product TIME15S-AUX. These are timestamps
        anchored at multiples of 15 seconds past the minute which correspond to the time
        base for the BOTSFLU data products which are binned on 15 seconds.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.

    Usage

        time15s = prs_botsflu_time15s(timestamp)

            where

        time15s = BOTSFLU-TIME15S-AUX [sec since 01-01-1900]
        timestamp = OOI system timestamps [sec since 01-01-1900]

    Notes:

        The BOTSFLU data products associated with this timebase are:
        MEANPRES
        PREDTIDE
        MEANDEPTH
        5MINRATE
        10MINRATE

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    ###### [in development]: the second calling argument is a placeholder
    time15s, _, _ = anchor_bin_rawdata_to_15s(timestamp, botpres)

    return time15s


def prs_botsflu_meanpres(timestamp, botpres):
    """
    Description:

        Calculates the BOTSFLU data product MEANPRES_L1.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.

    Usage

        meanpres = prs_botsflu_meanpres(timestamp, botpres)

            where

        meanpres = BOTSFLU-MEANPRES_L2 [psi]
        timestamp = OOI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The timebase data product associated with this data product is TIME15S.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    _, meanpres, _ = anchor_bin_rawdata_to_15s(timestamp, botpres)

    return meanpres


def prs_botsflu_predtide(time):
    """
    Description:

        Assigns tide values for the 3 BOTPT instrument sites about 500 km west of Astoria.

        When the input argument is the data product TIME15S, the output of this function
        will be the BOTSFLU data product PREDTIDE.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.

    Usage:

        PREDTIDE = prs_botsflu_predtide(TIME15S)

            where

        PREDTIDE = BOTSFLU-PREDTIDE data product [m]
        TIME15S = BOTSFLU-TIME15S data product [sec since 01-01-1900].

    Notes:

        Lookup table in binary file: 'ion_functions/data/prs_functions_tides_2014_thru_2019.mat'

        The lookup table contains tide values every 15 seconds from 2014-01-01 to 2020-01-01
        at lat = 45.95547 lon = -130.00957 calculated by the Tide Model Driver software
        written in Matlab (Mathworks, Natick, MA) using the TPXO7.2 global model. The tides
        corresponding to time are determined by positional indexing (the first value is for
        2014-01-01 00:00:00, the second is for 2014-01-01 00:00:15, etc). The 3 BOTPT sites
        are close enough together that the caldera center location can be used for all, as
        above: lat = 45.95547 lon = -130.00957.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.

        Matlab code to calculate tides using TPXO7.2 global model:
            http://polaris.esr.org/ptm_index.html
        Further documentation for the TPXO7.2 global tide model:
            http://volkov.oce.orst.edu/tides/global.html
    """
    time0 = 3597523200.0  # midnight, 2014-01-01
    time_interval = 15.0  # seconds

    # for unit test data, only, feb-apr 2011
    if time[0] < time0:
        time0 = 3502828800.0  # midnight, 2011-01-01
        matpath = 'data/matlab_scripts/botpt/tides_15sec_2011_for_unit_tests.mat'
    else:
        # else, OOI data from 2014 onwards
        matpath = 'data/prs_functions_tides_2014_thru_2019.mat'

    matstream = pkg_resources.resource_stream('ion_functions', matpath)
    dict_tides = spio.loadmat(matstream)
    # tide values are signed 4 byte integers, units [0.001mm]
    tidevector = 0.000001 * dict_tides['tides_mat']
    tidevector = tidevector.reshape((-1))
    # calculate tide vector index as a function of timestamp
    idx = np.around((time - time0) / time_interval)
    tide = tidevector[idx.astype(int)]
    return tide


def prs_botsflu_meandepth(timestamp, botpres):
    """
    Description:

        Calculates the BOTSFLU data product MEANDEPTH_L2, de-tided bottom depth
        as a function of time (15sec bins).

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage

        meandepth = prs_botsflu_meandepth(timestamp, botpres)

            where

        meandepth = BOTSFLU-MEANDEPTH_L2 [m]
        timestamp = OOI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The timebase data product associated with this data product is TIME15S.

        The DPS specifies that atmospheric pressure not be subtracted from the
        L1 pressure data even though its units are [psia].

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    _, meandepth, _ = calc_meandepth_plus(timestamp, botpres)

    return meandepth


def prs_botsflu_5minrate(timestamp, botpres):
    """
    Description:

        Calculates the BOTSFLU data product 5MINRATE_L2, the instantaneous rate of
        depth change using 5 minute backwards-looking meandepth data.

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage

        botsflu_5minrate = pprs_botsflu_5minrate(timestamp, botpres)

            where

        botsflu_5minrate = BOTSFLU-5MINRATE_L2 [cm/min]
        timestamp = CI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The timebase data product associated with this data product is TIME15S.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # calculate de-tided depth and the positions of non-zero bins in the original data.
    _, meandepth, mask_nonzero = calc_meandepth_plus(timestamp, botpres)

    # initialize data product including elements representing data gap positions
    botsflu_5minrate = np.zeros(mask_nonzero.size) + np.nan

    # re-constitute the original data, with data gaps represented by nans.
    data_w_gaps = np.copy(botsflu_5minrate)
    data_w_gaps[mask_nonzero] = meandepth

    # for 15s binned data, 5 minutes comes out to (5 minutes)/(0.25 min) = 20 intervals
    shift = 20
    # units of the subtraction are meter/5min; to convert to cm/min,
    # multiply by 100cm/m and divide by 5 = 20.
    botsflu_5minrate[shift:] = 20.0 * (data_w_gaps[shift:] - data_w_gaps[:-shift])

    # this rate product now has potentially two sources of nans;
    # definitely those at the start of the data record, and any that might
    # have been propagated into the calculation because of the presence of
    # data gaps. remove those only at the data dropout positions (if present)
    # so that this data product will have a 1:1 correspondence with
    # its associated timestamp variable (TIME15S).
    botsflu_5minrate = botsflu_5minrate[mask_nonzero]

    return botsflu_5minrate


def prs_botsflu_10minrate(timestamp, botpres):
    """
    Description:

        Calculates the BOTSFLU data product 10MINRATE_L2, the mean seafloor uplift rate
        calculated using 10 minute backwards-looking 10 minute running mean depth data.

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage

        botsflu_10minrate = pprs_botsflu_10minrate(timestamp, botpres)

            where

        botsflu_10minrate = BOTSFLU-10MINRATE_L2 [cm/hr]
        timestamp = OOI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The timebase data product associated with this data product is TIME15S.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # calculate de-tided depth and the positions of non-zero bins in the original data.
    _, meandepth, mask_nonzero = calc_meandepth_plus(timestamp, botpres)

    # initialize data product including elements representing data gap positions
    botsflu_10minrate = np.zeros(mask_nonzero.size) + np.nan

    # re-constitute the original data, with data gaps represented by nans.
    data_w_gaps = np.copy(botsflu_10minrate)
    data_w_gaps[mask_nonzero] = meandepth

    # now calculate sliding 10 minute means.
    # the mean of the 1st 40 values will be located at timestamp position 20
    # (python index 19).
    window_size = 40  # 10min averages on 0.25min binned data
    means = calculate_sliding_means(data_w_gaps, window_size)

    # as above, 10 minutes = 40 intervals for 15sec binned data.
    shift = 40
    # units of the subtraction are meter/10min; to convert to cm/hr,
    # multiply by 100cm/m and multiply by 6 = 600.
    botsflu_10minrate[shift:] = 600.0 * (means[shift:] - means[:-shift])

    # this rate product now has potentially two sources of nans;
    # definitely those at the start of the data record, and any that might
    # have been propagated into the calculation because of the presence of
    # data gaps. remove those only at the data dropout positions (if present)
    # so that this data product will have a 1:1 correspondence with
    # its associated timestamp variable (TIME15S).
    botsflu_10minrate = botsflu_10minrate[mask_nonzero]

    return botsflu_10minrate


def prs_botsflu_time24h(time15s):
    """
    Description:

        Calculates the auxiliary BOTSFLU data product TIME24H-AUX. These are
        timestamps anchored at midnight which correspond to the time base for
        the BOTSFLU data products which are binned on a day's worth of data.

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.
        2017-05-05: Russell Desiderio. Changed time24h time base to span the entire
                                       dataset including data gaps. This change is
                                       made in the function anchor_bin_detided_to_24h.

    Usage

        time24h = prs_botsflu_time24h(time15s)

            where

        time24h = BOTSFLU-TIME24H-AUX [sec since 01-01-1900]
        time15s = BOTSFLU-TIME15S-AUX [sec since 01-01-1900]

    Notes:

        The BOTSFLU data products associated with this timebase are:
        DAYDEPTH
        4WKRATE
        8WKRATE

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # the second calling argument is a placeholder
    time24h, _ = anchor_bin_detided_to_24h(time15s, None, None)

    return time24h


def prs_botsflu_daydepth(timestamp, botpres, dday_coverage=0.90):
    """
    Description:

        Calculates the BOTSFLU data product DAYDEPTH_L2, de-tided bottom depth
        as a function of time (1 day bins).

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.
        2017-05-05: Russell Desiderio. Changed time24h time base to span the entire
                                       dataset including data gaps.


    Usage

        daydepth = prs_botsflu_daydepth(timestamp, botpres)

            where

        daydepth = BOTSFLU-DAYDEPTH_L2 [m]
        timestamp = OOI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The timebase data product associated with this data product is TIME24H.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # calculate 15sec bin timestamps and de-tided depth.
    time15s, meandepth, _ = calc_meandepth_plus(timestamp, botpres)

    # bin the 15sec data into 24 hour bins so that the timestamps are at midnight.
    # to calculate daydepth, don't need the time24h timestamps.

    _, daydepth = anchor_bin_detided_to_24h(time15s, meandepth, dday_coverage)

    # downstream data products no longer require the mask_nonzero variable
    return daydepth

    #daydepth = calc_daydepth_plus(timestamp, botpres)
    #
    #return daydepth


def prs_botsflu_4wkrate(timestamp, botpres, dday_coverage=0.9, rate_coverage=0.75):
    """
    Description:

        Calculates the BOTSFLU data product 4WKRATE_L2, the mean rate of seafloor
        change as calculated by 4-week backwards-looking linear regressions.

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.
        2017-05-05: Russell Desiderio. Changed time24h time base to span the entire
                                       dataset including data gaps. Therefore removed
                                       the last masking operation that removed the
                                       values at bins that had zero data.

    Usage

        botsflu_4wkrate = pprs_botsflu_4wkrate(timestamp, botpres)

            where

        botsflu_4wkrate = BOTSFLU-4WKRATE_L2 [cm/yr]
        timestamp = CI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The timebase data product associated with this data product is TIME24H.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # calculate daydepth
    daydepth = prs_botsflu_daydepth(timestamp, botpres, dday_coverage)

    # 4 weeks of data
    window_size = 29
    botsflu_4wkrate = calculate_sliding_slopes(daydepth, window_size, rate_coverage)
    #  convert units:
    #    the units of the slopes are [y]/[x] = meters/day;
    #    to get units of cm/yr, multiply by 100cm/m * 365 days/yr
    botsflu_4wkrate = 100.0 * 365.0 * botsflu_4wkrate

    return botsflu_4wkrate


def prs_botsflu_8wkrate(timestamp, botpres, dday_coverage=0.9, rate_coverage=0.75):
    """
    Description:

        Calculates the BOTSFLU data product 8WKRATE_L2, the mean rate of seafloor
        change as calculated by 8-week backwards-looking linear regressions.

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.
        2017-05-05: Russell Desiderio. Changed time24h time base to span the entire
                                       dataset including data gaps. Therefore removed
                                       the last masking operation that removed the
                                       values at bins that had zero data.

    Usage

        botsflu_8wkrate = pprs_botsflu_8wkrate(timestamp, botpres)

            where

        botsflu_8wkrate = BOTSFLU-8WKRATE_L2 [cm/yr]
        timestamp = OOI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The timebase data product associated with this data product is TIME24H.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # calculate daydepth
    daydepth = prs_botsflu_daydepth(timestamp, botpres, dday_coverage)

   # 8 weeks of data
    window_size = 57
    botsflu_8wkrate = calculate_sliding_slopes(daydepth, window_size, rate_coverage)
    #  convert units:
    #    the units of the slopes are [y]/[x] = meters/day;
    #    to get units of cm/yr, multiply by 100cm/m * 365 days/yr
    botsflu_8wkrate = 100.0 * 365.0 * botsflu_8wkrate

    return botsflu_8wkrate


#**********************************************************************
#.. BOTSFLU: Worker functions called by the data product functions
#**********************************************************************
def anchor_bin_rawdata_to_15s(time, data):
    """
    Description:

        Calculates 'anchored' timestamps (see Notes) and binned data based on timestamps
        in units of seconds since midnight. Written explicitly for the BOTSFLU DPA which
        requires two stages of binning: 20hz data on 15 seconds, then the 15sec data on 24 hours.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.
        2015-01-14: Russell Desiderio. Changed output arguments and incorporated conditionals
                                       to improve program efficiency.

    Usage (1):

        bin_timestamps = anchor_bin(time, None, bin_duration, 'time')

            where

        bin_timestamps = 1D array of centered timestamps for non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        None = not used; python placeholder object
        bin_duration = size of bin [s]
        mode = the string 'time'

    Usage (2):

        binned_data, mask_nonzero = anchor_bin(time, data, bin_duration, 'data')

            where

        binned_data = 1D array of binned data; no empty bins are represented
        mask_nonzero = boolean where True values represent locations of non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        data = data to be binned
        bin_duration = size of bin [s]
        mode = the string 'data'

    Usage (3):

        bin_timestamps, binned_data, mask_nonzero = anchor_bin(time, data, bin_duration, 'both')

            where

        bin_timestamps = 1D array of centered timestamps for non-empty bins
        binned_data = 1D array of binned data; no empty bins are represented
        mask_nonzero = boolean where True values represent locations of non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        data = data to be binned
        bin_duration = size of bin [s]
        mode = the string 'both'

    Notes:

        The conditional construction is used so that only necessary statements are executed;
        when multiple years' worth of 20 Hz data is operated on, each np.bincount operation
        may take multiple tens of seconds to execute.

        The np.bincount routine is used in the same way accumarray in matlab is used
        to bin data. The key to the routine is to convert the timestamps into elapsed
        time in units of bin_duration and to construct bins based on the floored
        bin_duration times. The summing is then carried out by using the weighting
        feature of the np.bincount function, as described in the example in the
        numpy.bincount documentation as listed in the References.

        The BOTSFLU data products require binning at two stages. Bin results both with
        and without empty bins are required. The output arguments have been selected to
        provide this flexibility (in particular mask_nonzero).

        This routine has been constructed to supply 'anchored' timestamps. For example,
        if the bin_duration is 86400 (the number of seconds in a day) then the start time
        will be half a bin earlier than the first day of data (at noon) and all timestamps
        will be 'anchored' at midnight. Similarly, if the bin_duration is 15 sec, all
        timestamps will be at 00, 15, 30, and 45 seconds past the minute.

    References:

        http://docs.scipy.org/doc/numpy-1.8.1/reference/generated/numpy.bincount.html.
    """
    bin_duration = 15.0  # seconds
    half_bin = bin_duration/2.0

    # anchor time-centered bins by determining the start time to be half a bin
    # before the first 'anchor timestamp', which will be an integral number of
    # bin_durations after midnight.
    start_time = np.floor((time[0] - half_bin)/bin_duration) * bin_duration + half_bin
    # calculate elapsed time from start in units of bin_duration.
    time_elapsed = (time - start_time)/bin_duration
    # assign each timestamp a bin number index based on its elapsed time.
    bin_number = np.floor(time_elapsed).astype(int)
    # the number of elements in each bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    # create a logical mask of non-zero bin_count values
    mask_nonzero = (bin_count != 0)

    # directly calculate bin timestamp, units of [sec]:
    # the midpoint of the data interval is used.
    bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)
    # keep only the bins with values
    bin_timestamps = bin_timestamps[mask_nonzero]

    # sum the values in each time bin, and put into the variable binned_data
    binned_data = np.bincount(bin_number, data)
    # divide the values in non-empty bins by the number of values in each bin
    binned_data = binned_data[mask_nonzero]/bin_count[mask_nonzero]

    return bin_timestamps, binned_data, mask_nonzero


def anchor_bin_detided_to_24h(time, data, dday_coverage):
    """
    Description:

        Calculates 'anchored' timestamps (see Notes) and binned data based on timestamps
        in units of seconds since midnight. Written explicitly for the BOTSFLU DPA which
        requires two stages of binning: 20hz data on 15 seconds, then the 15sec data on 24 hours.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.
        2015-01-14: Russell Desiderio. Changed output arguments and incorporated conditionals
                                       to improve program efficiency.

    Usage (1):

        bin_timestamps = anchor_bin(time, None, bin_duration, 'time')

            where

        bin_timestamps = 1D array of centered timestamps for non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        None = not used; python placeholder object
        bin_duration = size of bin [s]
        mode = the string 'time'

    Usage (2):

        binned_data, mask_nonzero = anchor_bin(time, data, bin_duration, 'data')

            where

        binned_data = 1D array of binned data; no empty bins are represented
        mask_nonzero = boolean where True values represent locations of non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        data = data to be binned
        bin_duration = size of bin [s]
        mode = the string 'data'

    Usage (3):

        bin_timestamps, binned_data, mask_nonzero = anchor_bin(time, data, bin_duration, 'both')

            where

        bin_timestamps = 1D array of centered timestamps for non-empty bins
        binned_data = 1D array of binned data; no empty bins are represented
        mask_nonzero = boolean where True values represent locations of non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        data = data to be binned
        bin_duration = size of bin [s]
        mode = the string 'both'

    Notes:

        The conditional construction is used so that only necessary statements are executed;
        when multiple years' worth of 20 Hz data is operated on, each np.bincount operation
        may take multiple tens of seconds to execute.

        The np.bincount routine is used in the same way accumarray in matlab is used
        to bin data. The key to the routine is to convert the timestamps into elapsed
        time in units of bin_duration and to construct bins based on the floored
        bin_duration times. The summing is then carried out by using the weighting
        feature of the np.bincount function, as described in the example in the
        numpy.bincount documentation as listed in the References.

        The BOTSFLU data products require binning at two stages. Bin results both with
        and without empty bins are required. The output arguments have been selected to
        provide this flexibility (in particular mask_nonzero).

        This routine has been constructed to supply 'anchored' timestamps. For example,
        if the bin_duration is 86400 (the number of seconds in a day) then the start time
        will be half a bin earlier than the first day of data (at noon) and all timestamps
        will be 'anchored' at midnight. Similarly, if the bin_duration is 15 sec, all
        timestamps will be at 00, 15, 30, and 45 seconds past the minute.

    References:

        http://docs.scipy.org/doc/numpy-1.8.1/reference/generated/numpy.bincount.html.
    """
    bin_duration = 86400.0  # number of seconds in a day
    half_bin = bin_duration/2.0
    max_count = 86400.0/15.0  # maximum number of values in a day's bin

    # anchor time-centered bins by determining the start time to be half a bin
    # before the first 'anchor timestamp', which will an integral number of
    # bin_durations after midnight.
    start_time = np.floor((time[0] - half_bin)/bin_duration) * bin_duration + half_bin
    # calculate elapsed time from start in units of bin_duration.
    time_elapsed = (time - start_time)/bin_duration
    # assign each timestamp a bin number index based on its elapsed time.
    bin_number = np.floor(time_elapsed).astype(int)
    # the number of elements in each bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    # bin_count is used as a divisor to calculate mean values at each bin
    #    bins with bincounts below the threshold value will have a nan value
    bin_count[bin_count/max_count < dday_coverage] = np.nan
    #    replace 0 with nan
    bin_count[bin_count == 0] = np.nan

    # directly calculate bin timestamp, units of [sec]:
    # the midpoint of the data interval is used.
    bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)

    # sum the values in each time bin, and put into the variable binned_data
    binned_data = np.bincount(bin_number, data)
    # divide the values in each bin by the number of values in each bin
    daydepth = binned_data/bin_count

    return bin_timestamps, daydepth


def calc_meandepth_plus(timestamp, botpres):
    """
    Description:

        Worker function to calculate the botsflu data product meandepth plus
        additional variables required to calculate other botsflu data products
        downstream from meandepth.

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage

        time15s, meandepth, mask_nonzero = calc_meandepth_plus(timestamp, botpres)

            where

        time15s = TIME15S [sec since 01-01-1900]
        meandepth = BOTSFLU-MEANDEPTH_L2 [m]
        mask_nonzero = boolean of positions of non-empty bins in the original data
        timestamp = OOI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

        The DPS specifies that atmospheric pressure not be subtracted from the
        L1 pressure data even though its units are [psia].

        The DPS convention is that depths are negative, so that to detide the
        pressure record, the predicted tide is added to the negative depths.

        This function was written as a way to eliminate the execution of time
        consuming duplicate calculations in the botsflu coding within the
        OOI CI architecture constraints.

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # The pressure values do have units of psia. However, historically at these sites
    # atmospheric pressure has *not* been subtracted when converting the pressure data
    # to depth. Therefore the DPS authors do not want atmospheric pressure subtracted
    # in the DPA. To emphasize this, I have created the variable atm_press_psi and set
    # it to 0.
    atm_press_psi = 0.0
    psi_2_depth = -0.67  # psi to depth in meters
    bin_duration = 15.0  # seconds

    time15s, meanpres, mask_nonzero = anchor_bin_rawdata_to_15s(timestamp, botpres)
    # look up tide data
    tide = prs_botsflu_predtide(time15s)
    # de-tide
    meandepth = ((meanpres - atm_press_psi) * psi_2_depth) + tide

    # downstream data products require the time15s and mask_nonzero variables,
    # so pass these as output arguments so that they won't have to be recalculated.
    return time15s, meandepth, mask_nonzero


def calculate_sliding_means(data, window_size):
    """
    Description:

        Calculates time-centered means using digital convolution for the
        BOTSFLU data product 10MINRATE.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.
        2017-05-06: Russell Desiderio. Updated to make work with odd window_size.

    Usage

        means = calculate_sliding_means(data, window_size)

            where

        means = 1D array of sliding means
        data = 1D array of data
        window_size = window size, integer data type

    Notes

        The unit test values were calculated in Matlab, so that the python convolution
        result for even sized windows is shifted by 1 element to match the matlab result.

    """
    kk = np.ones(window_size) / window_size
    means = np.convolve(data, kk, 'same')
    # nan out data with boundary effects at edges.
    # integer arithmetic will 'truncate' 5/2 to 2 and -5/2 to -3.
    means[:window_size/2] = np.nan
    means[-((window_size-1)/2):] = np.nan

    # matlab and numpy behave differently for even window sizes, so roll the python
    # result to mimic matlab
    means = np.roll(means, -np.mod(window_size+1, 2))  # roll only if window_size is even

    return means


def calculate_sliding_slopes(data, window_size, coverage_threshold):
    """
    Description:

        Calculates backwards-looking sliding slopes using the normal linear
        regression equations rewritten to be less susceptible to round-off
        error; required for the BOTSFLU data products 4WKRATE and 8WKRATE.

    Implemented by:

        2017-05-03: Russell Desiderio. Initial code. Replaces the much faster Moore-Penrose
                                       pseudo-inverse method so that nan-masking can be
                                       incorporated.
        2017-05-08: Russell Desiderio. Added the 70% coverage criterion: if the number of non-Nan
                                       data points in a window is greater than or equal to 70% of
                                       the maximum possible window points calculate the slope.

    Usage

        slopes = calculate_sliding_slopes(data, window_size)

            where

        slopes = 1D array of sliding slopes
        data = 1D array of data
        window_size = integer
        coverage = fractional window fill threshold for calculation of slope values

    Notes

        The robust regression equations are taken from equations 14.2.15-14.2.17 in the Numerical
        Recipes reference below.

        Before the May 2017 modifications, just one Nan within a window would result in a Nan value
        for the slope. The routine now will calculate slopes by ignoring Nans, and if the coverage
        is 70% or better, a calculated value will result. If the coverage is less than 70%, then the
        output value is Nan.

        The data vector is padded so that bins within a window of the beginning of the data record
        that satisfy the coverage criterion will have non-Nan data product values.

    References:

        Press, Flannery, Teukolsky and Vetterling. Numerical Recipes, 1986; 1987 reprint.
        Cambridge University Press. page 507.
    """
    # ODD WINDOW SIZES are expected; if even, increment by one
    window_size = window_size + 1 - np.mod(window_size, 2)
    half_window = window_size / 2  # this will 'floor' when window_size is odd as desired

    # pad front end of data with Nans (because for loop is backwards-looking)
    npts = 2 * half_window + data.size
    padded_data = np.zeros(npts) + np.nan
    padded_data[window_size-1:] = data

    # first calculate values for all sliding windows
    slopes = np.zeros(npts) + np.nan
    abscissa = np.arange(npts).astype('float')
    for ii in range(window_size-1, npts):
        x = abscissa[(ii-2*half_window):(ii+1)]  # rather than np.arange-ing in each iteration
        y = padded_data[(ii-2*half_window):(ii+1)]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        if x.size < 2:
            continue  # trap out cases of either 0 or 1 valid datapoint in window
        tti = x - x.mean(axis=0)
        stt = np.sum(tti * tti)
        slopes[ii] = np.sum(tti * y) / stt
    # get rid of padding
    slopes = slopes[2*half_window:]

    # now determine the fraction of good values per window (vectorized):
    # first change non-nan values to 1, then nan values to 0, and take the average
    padded_data[~np.isnan(padded_data)] = 1.0
    padded_data[np.isnan(padded_data)] = 0.0
    fraction_good = calculate_sliding_means(padded_data, window_size)
    # convert to backwards-looking
    fraction_good = np.roll(fraction_good, half_window)
    # get rid of padding
    fraction_good = fraction_good[2*half_window:]

    # nan out fractional values (means) in sliding windows less than coverage threshold
    # there can be issues involving roundoff error when considering 100% coverage, so:
    machine_epsilon = np.finfo(np.float).eps
    fraction_good = fraction_good + 100 * machine_epsilon
    # avoid a python warning message by trapping out nans in the conditional
    fraction_good[np.isnan(fraction_good)] = -999.0  # any negative number will work as intended
    slopes[fraction_good < coverage_threshold] = np.nan

    return slopes


#**********************************************************************
#.. EVENT NOTIFICATION: tsunami detection
#**********************************************************************
def prs_tsunami_detection(botsflu_5minrate, tsunami_detection_threshold=1.0):
    """
    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage:

        TF = prs_tsunami_detection(BOTSFLU-5MINRATE_L2)

            where

            TF = True or False; whether a tsunami event has been detected.

    WARNING: This function and its data product input argument were coded as instructed
             in the DPS using the pseudocode specified. The robustness of this code has
             not been checked with actual data.
    """
    # units of variable and threshold are [cm/min]
    boolean_tsunami_detection = False
    # get rid of runtime warnings if nans are present
    botsflu_5minrate[np.isnan(botsflu_5minrate)] = 0.0
    if np.any(np.abs(botsflu_5minrate) >= tsunami_detection_threshold):
        boolean_tsunami_detection = True
    return boolean_tsunami_detection


#**********************************************************************
#.. EVENT NOTIFICATION: eruption imminent
#**********************************************************************
def prs_eruption_imminent(botsflu_10minrate, eruption_imminent_threshold=5.0):
    """
    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage:

        TF = prs_eruption_imminent(BOTSFLU-10MINRATE_L2)

            where

            TF = True or False; whether an eruption event is imminent.

    WARNING: This function and its data product input argument were coded as instructed
             in the DPS using the pseudocode specified. The robustness of this code has
             not been checked with actual data.
    """
    # units of variable and threshold are [cm/hr]
    boolean_eruption_imminent = False
    # get rid of runtime warnings if nans are present
    botsflu_10minrate[np.isnan(botsflu_10minrate)] = 0.0
    if np.any(botsflu_10minrate >= eruption_imminent_threshold):
        boolean_eruption_imminent = True
    return boolean_eruption_imminent


#**********************************************************************
#.. EVENT NOTIFICATION: eruption occurred
#**********************************************************************
def prs_eruption_occurred(botsflu_10minrate, eruption_occurred_threshold=-5.0):
    """
    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage:

        TF = prs_eruption_occurred(BOTSFLU-10MINRATE_L2)

            where

            TF = True or False; whether an eruption event has occurred.

    WARNING: This function and its data product input argument were coded as instructed
             in the DPS using the pseudocode specified. The robustness of this code has
             not been checked with actual data.
    """
    # units of variable and threshold are [cm/hr]
    boolean_eruption_occurred = False
    # get rid of runtime warnings if nans are present
    botsflu_10minrate[np.isnan(botsflu_10minrate)] = 0.0
    if np.any(botsflu_10minrate <= eruption_occurred_threshold):
        boolean_eruption_occurred = True
    return boolean_eruption_occurred


#**********************************************************************
#.. BOTSFLU functions deprecated May 2017 but retained.
#**********************************************************************


def anchor_bin_DEPRECATED(time, data, bin_duration, mode):
    """
    Description:

        Calculates 'anchored' timestamps (see Notes) and binned data based on timestamps
        in units of seconds since midnight. Written explicitly for the BOTSFLU DPA which
        requires two stages of binning: 20hz data on 15 seconds, then the 15sec data on 24 hours.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.
        2015-01-14: Russell Desiderio. Changed output arguments and incorporated conditionals
                                       to improve program efficiency.
        2017-05-05: Russell Desiderio. Deprecated because the new code requires different
                                       modifications to the rawdata and detided data binning:
                                       (1) bad value check in the rawdata
                                       (2) 'extended' 24hr timestamp records to incorporate
                                           non-Nan coverage thresholds for the detided data.

    Usage (1):

        bin_timestamps = anchor_bin(time, None, bin_duration, 'time')

            where

        bin_timestamps = 1D array of centered timestamps for non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        None = not used; python placeholder object
        bin_duration = size of bin [s]
        mode = the string 'time'

    Usage (2):

        binned_data, mask_nonzero = anchor_bin(time, data, bin_duration, 'data')

            where

        binned_data = 1D array of binned data; no empty bins are represented
        mask_nonzero = boolean where True values represent locations of non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        data = data to be binned
        bin_duration = size of bin [s]
        mode = the string 'data'

    Usage (3):

        bin_timestamps, binned_data, mask_nonzero = anchor_bin(time, data, bin_duration, 'both')

            where

        bin_timestamps = 1D array of centered timestamps for non-empty bins
        binned_data = 1D array of binned data; no empty bins are represented
        mask_nonzero = boolean where True values represent locations of non-empty bins
        time = 1D array of timestamps, units of sec since 01-01-1900
        data = data to be binned
        bin_duration = size of bin [s]
        mode = the string 'both'

    Notes:

        The conditional construction is used so that only necessary statements are executed;
        when multiple years' worth of 20 Hz data is operated on, each np.bincount operation
        may take multiple tens of seconds to execute.

        The np.bincount routine is used in the same way accumarray in matlab is used
        to bin data. The key to the routine is to convert the timestamps into elapsed
        time in units of bin_duration and to construct bins based on the floored
        bin_duration times. The summing is then carried out by using the weighting
        feature of the np.bincount function, as described in the example in the
        numpy.bincount documentation as listed in the References.

        The BOTSFLU data products require binning at two stages. Bin results both with
        and without empty bins are required. The output arguments have been selected to
        provide this flexibility (in particular mask_nonzero).

        This routine has been constructed to supply 'anchored' timestamps. For example,
        if the bin_duration is 86400 (the number of seconds in a day) then the start time
        will be half a bin earlier than the first day of data (at noon) and all timestamps
        will be 'anchored' at midnight. Similarly, if the bin_duration is 15 sec, all
        timestamps will be at 00, 15, 30, and 45 seconds past the minute.

    References:

        http://docs.scipy.org/doc/numpy-1.8.1/reference/generated/numpy.bincount.html.
    """
    half_bin = bin_duration/2.0
    # anchor time-centered bins by determining the start time to be half a bin
    # before the first 'anchor timestamp', which will an integral number of
    # bin_durations after midnight.
    start_time = np.floor((time[0] - half_bin)/bin_duration) * bin_duration + half_bin
    # calculate elapsed time from start in units of bin_duration.
    time_elapsed = (time - start_time)/bin_duration
    # assign each timestamp a bin number index based on its elapsed time.
    bin_number = np.floor(time_elapsed).astype(int)
    # the number of elements in each bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    # create a logical mask of non-zero bin_count values
    mask_nonzero = (bin_count != 0)

    # to calculate timestamps and to get tides, without also binning data.
    # mask_nonzero is not needed.
    if mode == 'time':
        # directly calculate bin timestamp, units of [sec]:
        # the midpoint of the data interval is used.
        bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)
        # keep only the bins with values
        bin_timestamps = bin_timestamps[mask_nonzero]
        return bin_timestamps

    # for binning data when the resultant timestamps are not explicitly required.
    # daydepth_plus also requires mask_nonzero for downstream products 4wkrate and 8wkrate.
    elif mode == 'data':
        # sum the values in each time bin, and put into the variable binned_data
        binned_data = np.bincount(bin_number, data)
        # divide the values in non-empty bins by the number of values in each bin
        binned_data = binned_data[mask_nonzero]/bin_count[mask_nonzero]
        return binned_data, mask_nonzero

    # for when both timestamps and binned data are required.
    elif mode == 'both':
        bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)
        bin_timestamps = bin_timestamps[mask_nonzero]
        binned_data = np.bincount(bin_number, data)
        binned_data = binned_data[mask_nonzero]/bin_count[mask_nonzero]
        return bin_timestamps, binned_data, mask_nonzero


def calc_daydepth_plus_deprecated(timestamp, botpres):
    """
    Description:

        Worker function to calculate the botsflu data product daydepth plus an
        additional boolean mask required to calculate other botsflu data products
        downstream from daydepth.

    Implemented by:

        2015-01-14: Russell Desiderio. Initial code.

    Usage

        daydepth, mask_nonzero = calc_daydepth_plus(timestamp, botpres)

            where

        daydepth = BOTSFLU-DAYDEPTH_L2 [m]
        mask_nonzero = boolean of positions of non-empty 24 hr bins
        timestamp = OOI system timestamps [sec since 01-01-1900]
        botpres = BOTPRES_L1 [psia]

    Notes:

    References:

        OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence
            (BOTSFLU) from the BOTPT instrument. Document Control Number 1341-00080.
    """
    # calculate 15sec bin timestamps and de-tided depth.
    time15s, meandepth, _ = calc_meandepth_plus(timestamp, botpres)

    # bin the 15sec data into 24 hour bins so that the timestamps are at midnight.
    # to calculate daydepth, don't need the time24h timestamps.

    _, daydepth = anchor_bin_detided_to_24h(time15s, meandepth)

    # downstream data products no longer require the mask_nonzero variable
    return daydepth


def calculate_sliding_slopes__MoorePenrose(data, window_size):
    # DEPRECATED because nan-masking cannot be implemented with this algorithm #
    """
    Description:

        Calculates backwards-looking sliding slopes using Moore_Penrose
        pseudo-inverse matrices; required for the BOTSFLU data products
        4WKRATE and 8WKRATE.

    Implemented by:

        2015-01-13: Russell Desiderio. Initial code.
        2017-05-03: Russell Desiderio. Deprecated because nan-masking cannot
                                       be implemented with this algorithm.

    Usage
 
        slopes = calculate_sliding_slopes(data, window_size)

            where

        slopes = 1D array of sliding slopes
        data = 1D array of data
        window_size = integer

    Notes

        Code lifted from John D'Errico's response on Matlab Central (thread 49181)
        to a query on how to calculate vectorized rolling regressions. For a more
        generalized application of the pinv\filter method, see D'Errico's 2007 code
        for movingslope.m on Matlab Central's file exchange (16997).

        The slopes are backwards-looking, not centered. The first non-nan value occurs
        at (python) index window_size, and is the slope of a regression of the first
        window_size points.
    """
    column1 = np.ones((window_size, 1))
    column2 = -np.arange(float(window_size)).reshape(-1, 1)
    X = np.hstack((column1, column2))
    filtercoef = np.linalg.pinv(X)
    slopes = signal.lfilter(filtercoef[1, :], 1, data)
    slopes[0:window_size-1] = np.nan

    # if time-centered slopes are desired, circularly shift this by half a window.

    return slopes
