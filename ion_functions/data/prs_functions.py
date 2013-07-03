#!/usr/bin/env python
"""
@package ion_functions.data.prs_functions
@file ion_functions/data/prs_functions.py
@author Christopher Wingard
@brief Module containing Seafloor Pressure related calculations.
"""

import numexpr as ne
import numpy as np


def prs_bottilt_ccmp(scmp, sn):
    """
    Description:

        OOI Level 1 Seafloor high-resolution tilt (BOTTILT) core data product,
        derived from data output by the Applied Geomechanincs LILY tilt sensor
        on board the Bottom Pressure Tilt (BOTPT) instruments on the Regional
        Scale Nodes (RSN) at Axial Seamount. This function computes
        BOTTILT-CCMP_L1.

    Implemented by:

        2013-06-10: Christopher Wingard. Initial code.

    Usage:

        ccmp = prs_bottilt_ccmp(scmp, sn)

            where

        ccmp = Corrected compass direction (BOTTILT-CCMP_L1) [degrees]
        scmp = Sensor compass direction (BOTTILT-SCMP_L0) [degrees]
        sn = LILY sensor serial number [unitless]

    References:

        OOI (2013). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-000060. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-000060_Data_Product_SPEC_BOTTILT_OOI.pdf)
    """
    # load the corrected compass directions table
    from ion_functions.data.prs_functions_ccmp import cmp_lookup

    # use the lookup table to get the ccmp
    ccmp = np.zeros(len(scmp))

    for i in range(len(scmp)):
        ccmp[i] = cmp_lookup[(sn[i], int(round(scmp[i])))]

    return ccmp


def prs_bottilt_tmag(x_tilt, y_tilt):
    """
    Description:

        OOI Level 1 Seafloor high-resolution tilt (BOTTILT) core data product,
        derived from data output by the Applied Geomechanincs LILY tilt sensor
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
        y_tilt = Sensor X_tilt (BOTTILT-YTLT_L0) [microradians]

    References:

        OOI (2013). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-000060. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-000060_Data_Product_SPEC_BOTTILT_OOI.pdf)
     """
    tmag = ne.evaluate('sqrt(x_tilt**2 + y_tilt**2)')
    return tmag


def prs_bottilt_tdir(x_tilt, y_tilt, ccmp):
    """
    Description:

        OOI Level 1 Seafloor high-resolution tilt (BOTTILT) core data product,
        derived from data output by the Applied Geomechanincs LILY tilt sensor
        on board the Bottom Pressure Tilt (BOTPT) instruments on the Regional
        Scale Nodes (RSN) at Axial Seamount.

    Implemented by:

        2013-06-10: Christopher Wingard. Initial code.

    Usage:

        tdir = prs_bottilt(x_tilt, y_tilt, ccmp)

            where

        tdir = Resultant tilt direction (BOTTILT-TDIR_L1) [degrees]
        x_tilt = Sensor X_tilt (BOTTILT-XTLT_L0) [microradians]
        y_tilt = Sensor X_tilt (BOTTILT-YTLT_L0) [microradians]
        ccmp = Sensor compass direction (BOTTILT-SCMP_L0) [degrees]

    References:

        OOI (2013). Data Product Specification for Seafloor High-Resolution
            Tilt. Document Control Number 1341-000060. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-000060_Data_Product_SPEC_BOTTILT_OOI.pdf)
     """
    ### Calculate the angle to use in the tilt direction formula
    # default angle calculation -- in degrees
    angle = ne.evaluate('arctan(y_tilt / x_tilt)')
    angle = np.degrees(angle)

    # if X-Tilt == 0 and Y-Tilt > 0
    mask = np.logical_and(x_tilt == 0, y_tilt > 0)
    angle[mask] = 90.0

    # if X-Tilt == 0 and Y-Tilt < 0
    mask = np.logical_and(x_tilt == 0, y_tilt < 0)
    angle[mask] = -90.0

    # if Y-Tilt == 0
    mask = np.equal(y_tilt, np.zeros(len(y_tilt)))
    angle[mask] = 0.0

    ### Calculate the tilt direction, using the X-Tilt to set the equation
    # default tilt direction equation
    tdir = ne.evaluate('(270 - angle + ccmp) % 360')

    # if X-Tilt >= 0
    tmp = ne.evaluate('(90 - angle + ccmp) % 360')
    mask = np.greater_equal(x_tilt, np.zeros(len(x_tilt)))
    tdir[mask] = tmp[mask]

    return np.round(tdir)
