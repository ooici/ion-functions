#!/usr/bin/env python

"""
@package ion_functions.data.do2_functions
@file ion_functions/data/do2_functions.py
@author Stuart Pearce
@brief Module containing Dissolved Oxygen family functions
"""
import numpy as np
import numexpr as ne
import pygsw.vectors as gsw


def do2_SVU(calphase, do_temp, csv):
    """
    Description:

        Stern-Volmer-Uchida equation for calculating temperature
        corrected dissolved oxygen concentration. OOI L1 data product.

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

    Implemented by:
        2013-04-26: Stuart Pearce. Initial Code.

    References:
        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)
    """

    # Calculate DO using Stern-Volmer:
    Ksv = csv[0] + csv[1]*do_temp + csv[2]*(do_temp**2)
    P0 = csv[3] + csv[4]*do_temp
    Pc = csv[5] + csv[6]*calphase
    DO = ne.evaluate('((P0/Pc) - 1) / Ksv')
    return DO


def do2_salinity_correction(DO, do_t, P, T, SP, lat, lon, pref=0):
    """
    Description:

        Salinity and pressure corrected dissolved oxygen concentration.
        OOI L2 data product DOCONCS.

    Usage:

        DOc = do2_salinity_correction(DO,do_t,P,T,SP,lat,lon, pref=0)

            where

        DOc = corrected dissolved oxygen [micro-mole/kg].
        DO = uncorrected dissolved oxygen [micro-mole/L].
        do_t = Oxygen sensor temperature [deg C].
        P = PRESWAT water pressure [dbar]. (see
            1341-00020_Data_Product_Spec_PRESWAT). Interpolated to the
            same timestamp as DO.
        T = TEMPWAT water temperature [deg C]. (see
            1341-00010_Data_Product_Spec_TEMPWAT). Interpolated to the
            same timestamp as DO.
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

    Implemented by:
        2013-04-26: Stuart Pearce. Initial Code.

    References:
         OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)

    """

    # density calculation from GSW toolbox
    SA = gsw.sa_from_sp(SP, P, lon, lat)
    CT = gsw.ct_from_t(SA, T, P)
    pdens = gsw.rho(SA, CT, pref)  # potential referenced to p=0

    # Convert from volume to mass units:
    DO = ne.evaluate('1000*DO/pdens')

    # Pressure correction:
    DO = ne.evaluate('(1 + (0.032*P)/1000) * DO')

    # Salinity correction:
    S0 = 0
    ts = ne.evaluate('log((298.15-do_t)/(273.15+do_t))')
    B0 = -6.24097e-3
    B1 = -6.93498e-3
    B2 = -6.90358e-3
    B3 = -4.29155e-3
    C0 = -3.11680e-7
    Bts = ne.evaluate('B0 + B1*ts + B2*ts**2 + B3*ts**3')
    DO = ne.evaluate('exp((SP-S0)*Bts + C0*(SP**2-S0**2)) * DO')
    return DO


def do2_dofst_volt(voltage_counts, Voffset, Soc, A, B, C, E, P, T, SP, lat, lon):
    """do2_dofst_volt

    Takes voltage counts measured from a DOFST-A (SBE 43) Oxygen sensor
    attached to a CTDPF-A (SBE 16+ V2) CTD, and converts the counts to a
    voltage and then voltage to dissolved oxygen in units of
    micromoles/kg for the OOI level 2 data product DOCONCF L2 (fast
    response oxygen) in combination with salinity, temperature, and
    pressure from the CTD.

    A Wrapper function for "dofst_calc".

    Usage:
        DO = do2_dofst_volt(volt_counts,Voffset,Soc,A,B,C,E,P,T,SP,lat,lon)

            where

        DO = corrected dissolved oxygen [micro-mole/kg].
        volt_counts = Oxygen sensor voltage [V].
        Voffset = Voltage offset [V].
        Soc = Oxygen signal slope
        A = Residual temperature correction factor A
        B = Residual temperature correction factor B
        C = Residual temperature correction factor C
        E = Pressure correction factor
        P = PRESWAT water pressure [dbar]. (see
            1341-00020_Data_Product_Spec_PRESWAT)
        T = TEMPWAT water temperature [deg C]. (see
            1341-00010_Data_Product_Spec_TEMPWAT)
        SP = PRACSAL practical salinity [unitless]. (see
            1341-00040_Data_Product_Spec_PRACSAL)
        lat, lon = latitude and longitude of the instrument [degrees].

    Example:
        v_counts = 16384
        P = 201.2
        T = 30.3
        SP = 31.2
        lat,lon = 39.0, -70.5
        A = -3.1867e-3, B = 1.7749e-4, C = -3.5718e-6
        E = 0.036, Voffset = -0.5186, Soc = 0.4396

        DO = do2_dofst_volt(v_counts,Voffset,Soc,A,B,C,E,P,T,SP,lat,lon)
        print DO
        > 61.89990653

    See Also: dofst_calc
    """
    # convert voltage counts to volts
    volts = voltage_counts / 13107.

    do, do_int = dofst_calc(volts, Voffset, Soc, A, B, C, E, P, T, SP, lat, lon)
    return do


def do2_dofst_frequency(frequency, Foffset, Soc, A, B, C, E, P, T, SP, lat, lon):
    """do2_dofst_frequency

    Takes a frequency measured from a DOFST-K (SBE 43F) Oxygen sensor
    connected to a CTDPF-CKL (SBE 52-MP) profiling CTD, and converts the
    frequency to dissolved oxygen in units of micromoles/kg for the OOI
    level 2 data product DOCONCF L2 (fast response oxygen) in
    combination with salinity, temperature, and pressure from the CTD.

    A Wrapper function for "dofst_calc".

    Usage:
        DO = do2_dofst_frequency(frequency,Foffset,Soc,A,B,C,E,P,T,SP,lat,lon)

            where

        DO = corrected dissolved oxygen [micro-mole/kg].
        frequency = Oxygen sensor frequency [Hz].
        Foffset = Frequency offset [Hz].
        Soc = Oxygen signal slope
        A = Residual temperature correction factor A
        B = Residual temperature correction factor B
        C = Residual temperature correction factor C
        E = Pressure correction factor
        P = PRESWAT water pressure [dbar]. (see
            1341-00020_Data_Product_Spec_PRESWAT)
        T = TEMPWAT water temperature [deg C]. (see
            1341-00010_Data_Product_Spec_TEMPWAT)
        SP = PRACSAL practical salinity [unitless]. (see
            1341-00040_Data_Product_Spec_PRACSAL)
        lat, lon = latitude and longitude of the instrument [degrees].

    Example:
        f = 4354
        P = 60.5200
        T = 15.5257
        SP = 34.1145
        lat,lon = 45.0, -125.0
        A = -4.1168e-3, B = 2.4818e-4, C = -3.8820e-6
        E = 0.036, Foffset = -839.55, Soc = 2.9968e-4

        DO = do2_dofst_frequency(f,Foffset,Soc,A,B,C,E,P,T,SP,lat,lon)
        print DO
        > 256.97434863158

    See Also: dofst_calc
    """
    do, do_int = dofst_calc(frequency, Foffset, Soc, A, B, C, E, P, T, SP, lat, lon)
    return do


# DOFST main sub-function
def dofst_calc(do_raw, offset, Soc, A, B, C, E, P, T, SP, lat, lon):
    """
    Description:

        Salinity and pressure corrected dissolved oxygen concentration.
        OOI L2 data product DOCONCF.

    Usage:

        DO = dostf_calculation(do_raw,offset,Soc,A,B,C,E,P,T,SP,lat,lon)

            where

        DO = corrected dissolved oxygen [micro-mole/kg].
        do_raw = Oxygen sensor voltage or frequency [V] or [Hz].
        offset = Voltage or Frequency offset [V] or [Hz].
        Soc = Oxygen signal slope
        A = Residual temperature correction factor A
        B = Residual temperature correction factor B
        C = Residual temperature correction factor C
        E = Pressure correction factor
        P = PRESWAT water pressure [dbar]. (see
            1341-00020_Data_Product_Spec_PRESWAT)
        T = TEMPWAT water temperature [deg C]. (see
            1341-00010_Data_Product_Spec_TEMPWAT)
        SP = PRACSAL practical salinity [unitless]. (see
            1341-00040_Data_Product_Spec_PRACSAL)
        lat, lon = latitude and longitude of the instrument [degrees].

    Example:
        do_raw = 4354  # frequency in Hz
        P = 60.5200
        T = 15.5257
        SP = 34.1145
        lat,lon = 45.0, -125.0
        A = -4.1168e-3, B = 2.4818e-4, C = -3.8820e-6
        E = 0.036, Foffset = -839.55, Soc = 2.9968e-4

        DO = dofst_calc(do_raw,Foffset,Soc,A,B,C,E,P,T,SP,lat,lon)
        print DO
        > 256.97434863158

    Implemented by:
        2013-08-20: Stuart Pearce. Initial Code.

    References:
         OOI (2013). Data Product Specification for Fast Dissolved
            Oxygen. Document Control Number 1341-00521.
            https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00521_Data_Product_SPEC_DOCONCF_OOI.pdf)
    """
    # Get potential density using the TEOS-10 toolbox
    SA = gsw.sa_from_sp(SP, P, lon, lat)
    pot_rho_t = gsw.pot_rho_t_exact(SA, T, P, 0)

    # Oxygen saturation value after Garcia and Gordon (1992)
    #   empirical polynomial coefficients (not calibration coeffs)
    A0 = 2.00907
    A1 = 3.22014
    A2 = 4.0501
    A3 = 4.94457
    A4 = -0.256847
    A5 = 3.88767
    B0 = -0.00624523
    B1 = -0.00737614
    B2 = -0.010341
    B3 = -0.00817083
    C0 = -0.000000488682
    temp_K = T + 273.15  # temperature in Kelvin
    Ts = np.log((298.15 - T) / (temp_K))
    Oxsol = np.exp(
        A0 + A1*Ts + A2*Ts**2 + A3*Ts**3 + A4*Ts**4 + A5*Ts**5 +
        SP * (B0 + B1*Ts + B2*Ts**2 + B3*Ts**3) +
        C0*SP**2)

    # Intermediate step: Dissolved Oxygen concentration in [mL/L]
    DO_int = Soc * (do_raw + offset) * Oxsol * (1.0 + A*T + B*T**2 + C*T**3) * np.exp((E * P)/temp_K)

    # Correct DO_int for Potential Density and convert to [micromole/Kg]
    DO = DO_int * 44660. / (pot_rho_t)
    return (DO, DO_int)
