#!/usr/bin/env python
"""
@package ion_functions.data.hyd_functions
@file ion_functions/data/hyd_functions.py
@author Christopher Wingard
@brief Module containing Hydrophone instrument family related functions
"""

import numexpr as ne
import numpy as np

def hyd_acoustic_pwaves(wav, gain):
    """
    Description:

        Calculates the OOI Level 1 (L1) Broadband Acoustic Pressure Waves core
        data product (HYDAPBB), using data from Broadband Hydrophone (HYDBB)
        instruments. The HYDBB instrument senses passive acoustic pressure
        waves from 5 Hz to 100 kHz, at 24-bit resolution.

    Implemented by:

        2014-05-16: Christopher Wingard. Initial Code

    Usage:

        tsv = hyd_acoustic_pwaves(wav, gain)

            where

        tsv = time-series voltage compensated for external gain and wav format
            scaling [Volts] (HYDAPBB_L1)
        wav = raw time-series voltage [Volts] (HYDAPBB_L0)
        gain = external gain setting [dB]

    References:

        OOI (2013). Data Product Specification for Acoustic Pressure Waves.
            Document Control Number 1341-00820.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00820_Data_Product_SPEC_HYDAPBB_OOI.pdf)
    """
    # shape inputs to correct dimensions
    wav = np.atleast_2d(wav)
    nRec = wav.shape[0]

    if np.isscalar(gain) is True:
        gain = np.tile(gain, (nRec, 1))
    else:
        gain = np.reshape(gain, (nRec, 1))

    # Convert the gain from dB to a linear value
    gain = ne.evaluate("10**(gain/20.)")

    # convert the broadband acoustic pressure wave data to Volts
    volts = ne.evaluate("wav * 3.")

    # and correct for the gain
    tsv = ne.evaluate("volts / gain")
    return tsv
