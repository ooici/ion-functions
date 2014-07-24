#!/usr/bin/env python
"""
@package ion_functions.data.hyd_functions
@file ion_functions/data/hyd_functions.py
@author Christopher Wingard
@brief Module containing Hydrophone instrument family related functions
"""
import numexpr as ne
import numpy as np


def hyd_bb_acoustic_pwaves(wav, gain):
    """
    Description:

        Calculates the OOI Level 1 (L1) Broadband Acoustic Pressure Waves core
        data product (HYDAPBB), using data from Broadband Hydrophone (HYDBB)
        instruments. The HYDBB instrument senses passive acoustic pressure
        waves from 5 Hz to 100 kHz, at 24-bit resolution.

    Implemented by:

        2014-05-16: Christopher Wingard. Initial Code

    Usage:

        tsv = hyd_bb_acoustic_pwaves(wav, gain)

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


def hyd_lf_acoustic_pwaves(raw, gain=3.2):
    """
    Description:

        Calculates the OOI Level 1 (L1) Low Frequency Acoustic Pressure Waves
        core data product (HYDAPLF), using data from the Low Frequency
        Hydrophone (HYDLF) instruments.

    Implemented by:

        2014-07-09: Christopher Wingard. Initial Code.

    Usage:

        hydaplf = hyd_lf_acoustic_pwaves(counts, gain)

            where

        hydaplf = time-series of low frequency acoustic pressure waves [V]
            (HYDAPLF_L1)
        raw = raw time-series digitizied in counts [counts] (HYDAPLF_L0)
        gain = Gurlap DM24 fixed gain bit weight [uV/count]

    References:

        OOI (2013). Data Product Specification for Low Frequency Acoustic
            Pressure Waves. Document Control Number 1341-00821.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00821_Data_Product_SPEC_HYDAPLF_OOI.pdf)
    """
    # apply the gain correction to convert the signal from counts to V
    gain = gain * 1.0e-6
    hydaplf = ne.evaluate("raw * gain")
    return hydaplf
