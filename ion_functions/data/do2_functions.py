#!/usr/bin/env python

"""
@package ion_functions.data.do2_functions
@file ion_functions/data/do2_functions.py
@author Stuart Pearce
@brief Module containing Dissolved Oxygen family functions
"""
import numpy as np
import pygsw as gsw


def do2_SVU(calphase, do_temp, csv):
    """
    Description:

        Stern-Volmer-Uchida equation for calculating temperature
        corrected dissolved oxygen concentration. OOI L1 data product. 


    Implemented by:

        2013-04-26: Stuart Pearce. Initial Code.

    Usage:

        DO = do2_SVU(calphase, do_temp, csv)

            where

        DO = dissolved oxygen [micro-mole/L]
        calphase = calibrated phase from an Oxygen sensor [deg]
            (see DOCONCS DPS)
        do_temp = oxygen sensor temperature [deg C],
            (see DOCONCS DPS)
        csv = Stern-Volmer-Uchida Calibration Coefficients array.
            7 element float array, (see DOCONCS DPS)

    Example:
        
        csv = np.array([0.002848, 0.000114, 1.51e-6, 70.42301, -0.10302,
                        -12.9462, 1.265377])
        calphase = 27.799
        do_temp = 19.841
        
        DO = do2_SVU(calphase, do_temp, csv)
        print DO
        > 363.900534505

    References: 
    
        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONOCS_OOI.pdf)
    """
    
    # Calculate DO using Stern-Volmer:
    Ksv = csv[0] + csv[1]*do_temp + csv[2]*(do_temp**2)
    P0 = csv[3] + csv[4]*do_temp
    Pc = csv[5] + csv[6]*calphase
    DO = ((P0/Pc) - 1) / Ksv
    return DO
    
# TODO: the salinity correction is not finished.  Ultimately it needs
# TODO: two encapsulating-in-time CTD samples to be interpolated to O2
# TODO: sample time.  Waiting on Luke and Chris M. to determine how best
# TODO: to do this.  Note: the interpolation should be in another function
# TODO: to accomodate the standalone DOSTA and CTDBP DOSTA.
def do2_salinity_correction(DO, do_t, P, T, SP, lat, lon, pref=0):
    """
    Description:

        Salinity and pressure corrected dissolved oxygen concentration.
        OOI L2 data product DOCONCS.

    Implemented by:

        2013-04-26: Stuart Pearce. Initial Code.

    Usage:

        DOc = do2_salinity_correction(DO,do_t,P,T,SP,lat,lon, pref=0)
        
            where
        
        DOc = corrected dissolved oxygen [micro-mole/kg].
        DO = uncorrected dissolved oxygen [micro-mole/L].
        do_t = Oxygen sensor temperature [deg C].
        P = PRESWAT water pressure [dbar]. (see
            1341-00020_Data_Product_Spec_PRESWAT)
        T = TEMPWAT water temperature [deg C]. (see
            1341-00010_Data_Product_Spec_TEMPWAT)
        SP = PRACSAL practical salinity [unitless]. (see
            1341-00040_Data_Product_Spec_PRACSAL)
        lat, lon = latitude and longitude of the instrument [degrees].
        pref = pressure reference level for potential density [dbar].
            The default is 0 dbar.

    Example:
        DO = 433.88488978325478
        do_t = 1.97
        P = 5.4000000000000004
        T = 1.97
        SP = 33.716000000000001
        lat,lon = -52.82, 87.64
        
        DOc = do2_salinity_correction(DO,do_t,P,T,SP,lat,lon, pref=0)
        print DO
        > 335.967894709

    References: 
    
         OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONOCS_OOI.pdf)
        
    """
    
    # density calculation from GSW toolbox
    SA = gsw.sa_from_sp(SP, P, lon, lat)
    CT = gsw.ct_from_t(SA, T, P)
    pdens = gsw.rho(SA,CT,pref)  # potential referenced to p=0
    
    # Convert from volume to mass units:
    DO = 1000*DO/pdens
    
    # Pressure correction:
    pcomp = 1 + (0.032*P)/1000
    DO = pcomp*DO
    
    # Salinity correction:
    S0 = 0
    ts = np.log((298.15-do_t)/(273.15+do_t))
    B = np.array([-6.24097e-3,
    -6.93498e-3,
    -6.90358e-3,
    -4.29155e-3])
    C0 = -3.11680e-7
    Bts = B[0] + B[1]*ts + B[2]*ts**2 + B[3]*ts**3
    scomp = np.exp((SP-S0)*Bts + C0*(SP**2-S0**2))
    DO = scomp*DO
    return DO