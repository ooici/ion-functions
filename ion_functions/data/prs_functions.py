#!/usr/bin/env python
"""
@package ion_functions.data.prs_functions
@file ion_functions/data/prs_functions.py
@author Christopher Wingard
@brief Module containing calculations related to instruments in the Seafloor
    Pressure family.
"""

import numexpr as ne
import numpy as np


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
            Tilt. Document Control Number 1341-000060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-000060_Data_Product_SPEC_BOTTILT_OOI.pdf)
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
            Tilt. Document Control Number 1341-000060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-000060_Data_Product_SPEC_BOTTILT_OOI.pdf)
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
            Tilt. Document Control Number 1341-000060.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-000060_Data_Product_SPEC_BOTTILT_OOI.pdf)
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
