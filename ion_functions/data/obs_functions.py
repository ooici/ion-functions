#!/usr/bin/env python
"""
@package ion_functions.data.obs_functions
@file ion_functions/data/obs_functions.py
@author Christopher Wingard
@brief Module containing Ocean Bottom Seismometer instrument related functions
"""
import numexpr as ne


def obs_bb_ground_velocity(raw, gain=3.2, sensitivity=1500.):
    """
    Description:

        Calculates the OOI Level 1 (L1) Broadband Ground Velocity core data
        product (GRNDVEL) in units of m/s, using data from the Ocean Bottom
        Broadband Seismometer (OBSBB and OBSBK) instruments.

    Implemented by:

        2014-07-09: Christopher Wingard. Initial Code

    Usage:

        grndvel = obs_bb_ground_velocity(counts, gain, sensitivity)

            where

        grndvel = time-series of ocean bottom seismic signal [m/s] (GRNDVEL_L1)
        raw = raw time-series digitizied in counts [counts] (GRNDVEL_L0)
        gain = Gurlap DM24 fixed gain bit weight [uV/count]
        sensitivity = Gurlap CMG1T sensor sensitivity [V/m/s]

    References:

        OOI (2013). Data Product Specification for Broadband Ground Velocity.
            Document Control Number 1341-00090.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00090_Data_Product_SPEC_GRNDVEL_OOI.pdf)
    """
    # scale the gain and sensitivity ...
    gain = gain * 1.0e-6
    sense = 2. * sensitivity

    # ... and calculate the broadband ground velocity
    grndvel = ne.evaluate("raw * (gain / sense)")
    return grndvel


def obs_bb_ground_acceleration(raw, gain=3.2, sensitivity=0.508):
    """
    Description:

        Calculates the OOI Level 1 (L1) Broadband Ground Acceleration core data
        product (GRNDACC) in units of m/s^2, using data from the Ocean Bottom
        Broadband Seismometer (OBSBB and OBSBK) instruments.

    Implemented by:

        2014-07-09: Christopher Wingard. Initial Code

    Usage:

        grndacc = obs_bb_ground_acceleration(counts, gain, sensitivity)

            where

        grndacc = time-series of ocean bottom seismic signal [m/s^2] (GRNDACC_L1)
        raw = raw time-series digitizied in counts [counts] (GRNDACC_L0)
        gain = Gurlap DM24 fixed gain bit weight [uV/count]
        sensitivity = Gurlap CMG5T sensor sensitivity [V/m/s^2]

    References:

        OOI (2013). Data Product Specification for Broadband Ground Acceleration.
            Document Control Number 1341-00100.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00100_Data_Product_SPEC_GRNDACC_OOI.pdf)
    """
    # scale the gain and sensitivity ...
    gain = gain * 1.0e-6
    sense = 2. * sensitivity

    # ... and calculate the broadband ground acceleration
    grndacc = ne.evaluate("raw * (gain / sense)")
    return grndacc


def obs_sp_ground_velocity(raw, gain=2.84, sensitivity=1200.):
    """
    Description:

        Calculates the OOI Level 1 (L1) Short Period Ground Velocity core data
        product (SGRDVEL) in units of m/s, using data from the Ocean Bottom
        Short Period Seismometer (OBSSP) instruments.

    Implemented by:

        2014-07-09: Christopher Wingard. Initial Code

    Usage:

        sgrdvel = obs_sp_ground_velocity(counts, gain, sensitivity)

            where

        grndacc = time-series of ocean bottom seismic signal [m/s] (SGRDVEL_L1)
        raw = raw time-series digitizied in counts [counts] (SGRDVEL_L0)
        gain = Gurlap DM24 fixed gain bit weight [uV/count]
        sensitivity = Gurlap CMG6T sensor sensitivity [V/m/s]

    References:

        OOI (2013). Data Product Specification for Broadband Ground Acceleration.
            Document Control Number 1341-00100.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00100_Data_Product_SPEC_GRNDACC_OOI.pdf)
    """
    # scale the gain and sensitivity ...
    gain = gain * 1.0e-6
    sense = 2. * sensitivity

    # ... and calculate the short period ground velocity
    sgrdvel = ne.evaluate("raw * (gain / sense)")
    return sgrdvel
