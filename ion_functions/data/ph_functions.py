#!/usr/bin/env python

"""
@package ion_functions.data.ph_functions
@file ion_functions/data/ph_functions.py
@author Christopher Wingard
@brief Module containing pH family instrument related functions
"""

# imports
import numpy as np
import numexpr as ne
import scipy as sp


# wrapper functions to extract L0 parameters from SAMI-II pH instruments (PHSEN)
def ph_434_intensity(light):
    """
    Wrapper function to extract the signal intensity at 434 nm (PH434SI_L0)
    from the pH instrument light measurements. Coded to accept either a
    single record or an array of records.
    """
    light = np.atleast_3d(light).astype(np.float)
    new = np.reshape(light, (-1, 23, 4))
    si434 = new[:, :, 1]
    return si434  # signal intensity, 434 nm (PH578SI_L0)


def ph_578_intensity(light):
    """
    Wrapper function to extract the signal intensity at 578 nm (PH578SI_L0)
    from the pH instrument light measurements. Coded to accept either a
    single record or an array of records.
    """
    light = np.atleast_3d(light).astype(np.float)
    new = np.reshape(light, (-1, 23, 4))
    si578 = new[:, :, 3]
    return si578  # signal intensity, 578 nm (PH578SI_L0)


def ph_thermistor(traw):
    """
    Wrapper function to convert the thermistor data (ABSTHRM_L0) from counts to
    degrees Centigrade from the pH instrument.
    """
    # convert raw thermistor readings from counts to degrees Centigrade
    Rt = ne.evaluate('(traw / (4096.0 - traw)) * 17400.0')
    lRt = np.log(Rt) 
    InvT = ne.evaluate('0.0010183 + 0.000241 * lRt + 0.00000015 * lRt**3')
    therm = ne.evaluate('(1.0 / InvT) - 273.15')

    return therm


def ph_battery(braw):
    """
    Wrapper function to convert the battery voltage from counts to Volts from
    the pH instrument.
    """
    # convert raw battery readings from counts to Volts
    volts = ne.evaluate('braw * 15. / 4096.')
    return volts


#def ph_calc_phwater(ref, light, therm, ea434, eb434, ea578, eb578, psal=35.0):
#    """
#    Wrapper function to vectorize the Sea Water pH calculation defined
#    below in ph_phwater. Coded to accept either a single record or an array of
#    records.
#    """
#    if ref.ndim == 1:
#        # only one record is being presented to the function
#        pH = ph_phwater(ref, light, therm, ea434, eb434, ea578, eb578, psal)
#    else:
#        # multiple records are being passed in, determine how many and create
#        # an empty pH array to hold the final outputs.
#        nRec = ref.shape[0]
#        pH = np.zeros(nRec, dtype=np.float)
#
#        # Salinity can either be a default value of 35 or from a co-located
#        # CTD. Need to replicate it to the same size as nRec if it is just the
#        # default scalar.
#        if np.isscalar(psal) is True:
#            psal = np.tile(psal, (nRec, 1))
#
#        # Iterate through the records
#        for iRec in range(nRec):
#            pH[iRec] = ph_phwater(ref[iRec, :], light[iRec, :], therm[iRec],
#                                  ea434[iRec], eb434[iRec], ea578[iRec],
#                                  eb578[iRec], psal[iRec])
#
#    return pH


def ph_calc_phwater(ref, light, therm, ea434, eb434, ea578, eb578, psal=35.0, ind=1):
    """
    Description:

        OOI Level 1 pH of seawater core data product, which is calculated using
        data from the Sunburst SAMI-II pH instrument (PHSEN). This document is
        intended to be used by OOI programmers to construct appropriate
        processes to create the L1 pH of seawater core data product.

    Implemented by:

        2013-04-19: Christopher Wingard. Initial code.

    Usage:

        ph = ph_phwater(ref, light, therm, ea434, eb434, ea578, eb578, psal=35)

            where

        ph = measured pH of seawater [unitless]
        ref = raw signal and reference measurements during blank cycle [counts]
        light = raw signal and reference measurements during measurement cycle
            [counts]
        therm = thermistor reading at end of measurement cycle [deg_C]
        ea434 = mCP molar absorptivities provided by vendor, specific to a
            reagent bag with a defined shelflife.
        eb434 = mCP molar absorptivities as above
        ea578 = mCP molar absorptivities as above
        eb578 = mCP molar absorptivities as above
        psal = practical salinity estimate used in calculcations, default is
            35.0 [unitless], from a co-located CTD.
        ind = indicator impurity correction, default is 1 [unitless]

    References:

        OOI (2012). Data Product Specification for pH of Seawater. Document
            Control Number 1341-00510. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00510_Data_Product_SPEC_PHWATER_OOI.pdf)
    """
    # reformat all input values to arrays of the correct dimensions, shape, and
    # type, recording the number of input records.
    ref = (np.atleast_2d(ref)).astype(np.float)
    nRec = ref.shape[0]

    light = np.atleast_3d(light).astype(np.float)
    light = np.reshape(light, (nRec, 23, 4))

    therm = np.reshape(therm, (nRec, 1)).astype(np.float)

    ea434 = np.reshape(ea434, (nRec, 1)).astype(np.float)
    eb434 = np.reshape(eb434, (nRec, 1)).astype(np.float)
    ea578 = np.reshape(ea578, (nRec, 1)).astype(np.float)
    eb578 = np.reshape(eb578, (nRec, 1)).astype(np.float)

    if np.isscalar(ind) is True:
        ind = np.tile(ind, nRec).astype(np.int)

    if np.isscalar(psal) is True:
        psal = np.tile(psal, (nRec, 1)).astype(np.float)
    else:
        psal = np.reshape(psal, (nRec, 1)).astype(np.float)

    # Calculate blanks from the 16 sets of reference light measurements
    arr434 = np.array([
        (ref[:, 1] / ref[:, 0]),
        (ref[:, 5] / ref[:, 4]),
        (ref[:, 9] / ref[:, 8]),
        (ref[:, 13] / ref[:, 12]),
    ])
    blank434 = np.reshape(np.mean(arr434, axis=0), (nRec, 1))

    arr578 = np.array([
        (ref[:, 3] / ref[:, 2]),
        (ref[:, 7] / ref[:, 6]),
        (ref[:, 11] / ref[:, 10]),
        (ref[:, 15] / ref[:, 14]),
    ])
    blank578 = np.reshape(np.mean(arr578, axis=0), (nRec, 1))

    # Extract 23 sets of 4 light measurements into arrays corresponding to the
    # raw reference and signal measurements at 434 and 578 nm. Input is an
    # array of length 92 (23 sets * 4 measurements per set). Can reshape and
    # slice to extract the parameters.
    ref434 = light[:, :, 0]   # reference signal, 434 nm
    int434 = light[:, :, 1]   # signal intensity, 434 nm (PH434SI_L0)
    ref578 = light[:, :, 2]   # reference signal, 578 nm
    int578 = light[:, :, 3]   # signal intensity, 578 nm (PH578SI_L0)

    # Absorbance
    A434 = -sp.log10(int434 / ref434)
    A434blank = -sp.log10(blank434)
    abs434 = A434 - A434blank

    A578 = -sp.log10(int578 / ref578)
    A578blank = -sp.log10(blank578)
    abs578 = A578 - A578blank

    R = abs578 / abs434

    # pka from Clayton and Byrne, 1993
    pKa = (1245.69 / (therm + 273.15)) + 3.8275 + (0.0021 * (35. - psal))
    pKa = np.reshape(pKa, (-1, 1))

    # Molar absorptivities
    Ea434 = ea434 - (26. * (therm - 24.788))
    Ea578 = ea578 + (therm - 24.788)
    Eb434 = eb434 + (12. * (therm - 24.788))
    Eb578 = eb578 - (71. * (therm - 24.788))
    e1 = Ea578 / Ea434
    e2 = Eb578 / Ea434
    e3 = Eb434 / Ea434

    V1 = R - e1
    V2 = e2 - R * e3

    # indicator concentration calculations
    HI = (abs434 * Eb578 - abs578 * Eb434) / (Ea434 * Eb578 - Eb434 * Ea578)
    I = (abs578 * Ea434 - abs434 * Ea578) / (Ea434 * Eb578 - Eb434 * Ea578)
    IndConc = HI + I
    pointph = np.real(pKa + sp.log10(V1 / V2))

    # ************************ Initial pH Calcs ************************
    # determine the most linear region of points for pH of seawater
    # calculation, skipping the first 5 points.
    IndConca = IndConc[:, 5:]
    Y = pointph[:, 5:]
    X = np.linspace(1, 18, 18)

    # create arrays for vectorized computations used in sum of squares below.
    # reflows 1D and 2D arrays into 2D and 3D arrays, respectively, shifting
    # each "row" of the arrays by one value, allowing the sum of square
    # calculations to be computed in a vectorized fashion, replacing the for
    # loop that had the computations running on 1:8, 2:9, ... 11:18.
    step = 7  # number of points to use
    count = step + 1
    nPts = np.size(X) - step
    x = np.zeros((nPts, count))
    y = np.zeros((nRec, nPts, count))
    for i in range(nPts):
        x[i, :] = X[i:i+count]
        for j in range(nRec):
            y[j, i, :] = Y[j, i:i+count]

    # compute the range of best fitting points, using array multiplications to
    # determine the best fit via the correlation coefficient.
    sumx = np.sum(x, axis=1)
    sumy = np.sum(y, axis=2)
    sumxy = np.sum(x * y, axis=2)
    sumx2 = np.sum(x**2, axis=1)
    sumy2 = np.sum(y**2, axis=2)
    sumxx = sumx * sumx
    sumyy = sumy * sumy
    ssxy = sumxy - (sumx * sumy) / count
    ssx = sumx2 - (sumxx / count)
    ssy = sumy2 - (sumyy / count)
    r2 = ssxy**2 / (ssx * ssy)

    # Range of seawater points to use
    cutoff1 = np.argmax(r2, axis=1)  # Find the first, best R-squared value
    cutoff2 = cutoff1 + count

    # Indicator and pH range limited to best points
    IndConcS = np.zeros((nRec, count))
    pointphS = np.zeros((nRec, count))
    for i in range(nRec):
        IndConcS[i, :] = IndConca[i, cutoff1[i]:cutoff2[i]]
        pointphS[i, :] = Y[i, cutoff1[i]:cutoff2[i]]

    # ************************* Final pH Calcs *************************
    sumx = np.sum(IndConcS, axis=1)
    sumy = np.sum(pointphS, axis=1)
    sumxy = np.sum(pointphS * IndConcS, axis=1)
    sumx2 = np.sum(IndConcS**2, axis=1)
    sumy2 = np.sum(pointphS**2, axis=1)
    xbar = np.mean(IndConcS, axis=1)
    ybar = np.mean(pointphS, axis=1)
    sumxx = sumx * sumx
    sumyy = sumy * sumy
    ssxy = sumxy - (sumx * sumy) / count
    ssx = sumx2 - (sumxx / count)
    ssy = sumy2 - (sumyy / count)
    slope = ssxy / ssx
    ph = ybar - slope * xbar

    # pH corrections due to indicator impurity
    ind0Flag = (ind == 0) & (ph >= 8.2)
    ind1Flag = (ind == 1) & (ph >= 8.2)
    # If the indicator impurity == 0 and the pH is >= 8.2
    ph[ind0Flag] = ph[ind0Flag] * 0.9723 + 0.2235
    # If the indicator impurity == 1 and the pH is >= 8.2
    ph[ind1Flag] = ph[ind1Flag] * 0.9698 + 0.2484

    return ph
