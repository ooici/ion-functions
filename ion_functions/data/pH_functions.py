#!/usr/bin/env python

"""
@package ion_functions.data.pH_functions
@file ion_functions/data/pH_functions.py
@author Christopher Wingard
@brief Module containing pH family instrument related functions
"""

def pH_phwater(ref, light, tstrt, tend, psal=35.0):
    """
    Description:

        OOI Level 1 pH of seawater core data product, which is calculated using
        data from the Sunburst SAMI-II pH instrument (PHSEN). This document is
        intended to be used by OOI programmers to construct appropriate
        processes to create the L1 pH of seawater core data product. 


    Implemented by:

        2013-04-19: Christopher Wingard. Initial code.

    Usage:

        pH, tfinal = ph_phwater(ref, light, traw, psal=35)

            where

        pH = measured pH of seawater [unitless]
        tfinal = temperature measured at end of cycle [deg_C]
        ref = raw signal and reference measurements during blank cycle [counts] 
        light = raw signal and reference measurements during measurement cycle
            [counts]
        traw = raw thermistor reading at end of measurement cycle [counts]
        psal = practical salinity estimate used in calculcations, default is
            35.0 [unitless]
        
    
    References: 
    
        OOI (2012). Data Product Specification for pH of Seawater. Document
            Control Number 1341-00510. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00510_Data_Product_SPEC_PHWATER_OOI.pdf)
    """
    from scipy import stats
    import numpy as np
    
    # set constants
    cp = 1. # cell path length
    
    # [TODO] these are actually inputs and are instrument/reagent bag specific
    ea434 = 17709.
    ea578 = 107.
    eb434 = 2287.
    eb578 = 38913.

    # calculate blanks from the 16 sets of reference light measurements 
    arr434 = np.array([
        (ref[1] / ref[0]),
        (ref[5] / ref[4]),
        (ref[9] / ref[8]),
        (ref[13] / ref[12]),
    ])
    blank434 = np.mean(arr434)

    arr578 = np.array([
        (ref[3] / ref[2]),
        (ref[7] / ref[6]),
        (ref[11] / ref[10]),
        (ref[15] / ref[14]),
    ])
    blank578 = np.mean(arr578)

    # convert the thermistor readings, taken before and after the end of the
    # measurement cycle, from counts to degrees Centigrade.
    Rt = (tstrt / (4096. - tstrt)) * 17400.
    InvT = 0.0010183 + 0.000241 * np.log(Rt) + 0.00000015 * np.log(Rt)**3
    TempK = 1 / InvT
    tfinal1 = TempK - 273.15
 
    Rt = (tend / (4096. - tend)) * 17400.
    InvT = 0.0010183 + 0.000241 * np.log(Rt) + 0.00000015 * np.log(Rt)**3
    TempK = 1 / InvT
    tfinal2 = TempK - 273.15

    # per comments in DPS, use the final temperature measurement
    # tfinal = np.mean([tfinal1, tfinal2])
    tfinal = tfinal2
    
    # extract 23 sets of 4 light measurements into arrays corresponding to the
    # raw reference and signal measurements at 434 and 578 nm. Input is an
    # array of length 92 (23 sets * 4 measurements per set). Can reshape and
    # slice to extract the parameters.
    new = np.reshape(light,(23,4))
    ref434 = new[:,0]   # reference signal, 434 nm
    int434 = new[:,1]   # signal intensity, 434 nm (PH434SI_L0)
    ref578 = new[:,2]   # reference signal, 578 nm
    int578 = new[:,3]   # signal intensity, 578 nm (PH578SI_L0) 
  
    # Absorbance
    A434 = -np.log10(int434 / ref434)
    A434blank = -np.log10(blank434)
    abs434 = A434 - A434blank
    
    A578 = -np.log10(int578 / ref578)
    A578blank = -np.log10(blank578)
    abs578 = A578 - A578blank

    # pka from Clayton and Byrne, 1993    
    pKa = (1245.69 / (tfinal + 273.15)) + 3.8275 + (0.0021 * (35 - psal))
    R = (A578 - A578blank) / (A434 - A434blank)
    
    # Molar absorptivities
    inta434 = ea434 + 24.5250 * 24.8000
    inta578 = ea578 - 0.5869 * 24.8000
    intb434 = eb434 - 6.2231 * 24.8600
    intb578 = eb578 + 99.6170 * 24.8600
    ea434 = -24.525 * tfinal + inta434
    ea578 = 0.5869 * tfinal + inta578
    eb434 = 6.2231 * tfinal + intb434
    eb578 = -99.6170 * tfinal + intb578
    #ea434 = ea434 - (26 * (tfinal - 24.788))
    #ea578 = ea578 + (tfinal - 24.788)
    #eb434 = eb434 + (12 * (tfinal - 24.788))
    #eb578 = eb578 - (71 * (tfinal - 24.788))
    e1 = ea578 / ea434
    e2 = eb578 / ea434
    e3 = eb434 / ea434
    
    V1 = R - e1
    V2 = e2 - R * e3 
    
    # indicator concentration calculations
    HI = (abs434 * eb578 - abs578 * eb434) / (ea434 * eb578 - eb434 * ea578)
    I = (abs578 * ea434 - abs434 * ea578) / (ea434 * eb578 - eb434 * ea578)
    IndConc = HI + I
    #pointpH = np.real(pKa + np.log10(V1 / V2))
    pointpH = pKa + np.log10(V1 / V2)
    print pointpH
    
    # ************************ Initial pH Calcs ************************
    # determine the most linear region of points for pH of seawater
    # calculation, skipping the first 5 points.
    IndConca = IndConc[5:]
    Y = pointpH[5:]
    X = np.linspace(1, 18, 18)
    
    step = 7 # number of points to use 
    npts = np.size(X) - step
    slp = np.zeros(npts)
    r2 = np.zeros(npts)
    for i in range(npts):
        m, b, r, p, serr = stats.linregress(X[i:i+step], Y[i:i+step])
        slp[i]= m
        r2[i] = r**2
        
    # Range of seawater points to use
    cutoff1 = np.argmax(r2)  # Find the first, best R-squared value
    cutoff2 = cutoff1 + step
    
    # Indicator and pH range limited to best points
    IndConcS = IndConca[cutoff1:cutoff2]
    pointpHS = np.real(Y[cutoff1:cutoff2])

    # ************************* Final pH Calcs *************************
    xbar = np.mean(IndConcS)
    ybar = np.mean(pointpHS)
    m, b, r, p, serr = stats.linregress(IndConcS, pointpHS)
    pH = ybar - m * xbar

    return pH, tfinal