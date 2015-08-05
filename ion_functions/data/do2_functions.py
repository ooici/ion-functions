#!/usr/bin/env python

"""
@package ion_functions.data.do2_functions
@file ion_functions/data/do2_functions.py
@author Stuart Pearce, Russell Desiderio
@brief Module containing Dissolved Oxygen family functions
"""
import numpy as np
import numexpr as ne
import pygsw.vectors as gsw

from ion_functions.data.generic_functions import replace_fill_with_nan

"""
DOSTA Processing Configurations:

    DOSTA configured for analog output to CTD voltage channels:
        T_optode_degC = dosta_Topt_volt_to_degC(T_optode_volts)
        DOCONCS-DEG_L0 = dosta_phase_volt_to_degree(DOCONCS-VLT_L0)
        DOCONCS_L1 = do2_SVU(DOCONCS-DEG_L0, T_optode_degC, ...)
        DOXYGEN_L2 = do2_salinity_correction(DOCONCS_L1, ...)

    DOSTA configured for digital output to CTD RS-232:
        DOCONCS_L1 = o2_counts_to_uM(DOCONCS-CNT_L0)
        DOXYGEN_L2 = do2_salinity_correction(DOCONCS_L1, ...)

    DOSTA, autonomous operation, digital output:
        DOCONCS_L1 = do2_SVU(DOCONCS-DEG_L0, T_optode_degC, ...)
        DOXYGEN_L2 = do2_salinity_correction(DOCONCS_L1, ...)


DOSTA DATA PRODUCTS:

    DOCONCS-CNT_L0 [counts]:  oxygen concentration, uncorrected for salinity and pressure.
                              (a) parsed from digital DOSTA output when routed through RS-232 CTD.
    DOCONCS-DEG_L0 [degrees]: CalPhase (calibrated phase). Two sources:
                              (a) parsed from autonomous DOSTA digital output.
                              (b) calculated by dosta_phase_volt_to_degree from DOCONCS-VLT_L0
                                  when the DOSTA is connected to a CTD analog voltage channel.
    DOCONCS-VLT_L0 [volts]:   CalPhase (calibrated phase).
                              (a) parsed from analog DOSTA output when routed through CTD.

    DOCONCS_L1 [micro-mole/liter]: oxygen concentration, uncorrected for salinity and pressure.

                                   As of Aug 2015, SAF has not been updated to reflect the
                                   change of units from micro-mole/kg to micro-mole/liter.

                                   Two sources:
                                   (a) calculated by o2_counts_to_uM from DOCONCS-CNT_L0.
                                   (b) calculated by do2_SVU from DOCONCS-DEG_L0.

    DOXYGEN_L2 [micro-mole/kg]:    oxygen concentration, corrected for salinity and pressure
                                   (and inherently temperature).
                                   (a) calculated by do2_salinity_correction from DOCONCS_L1.


    Temperature in DOSTA data product calculations:

    The temperature input to the function do2_SVU should be the sensor's optode foil temperature
    [degC] measured by the optode's thermistor (Topt in the DPS). This variable is directly parsed
    from digital autonomous DOSTA data streams; if instead the DOSTA is connected to a CTD through
    an analog voltage channel, then Topt[degC] is calculated from Topt[V] by the function
    dosta_Topt_volt_to_degC.

    The temperature input to the function do2_salinity should be TEMPWAT from the co-located CTD.


DOFST DATA PRODUCTS:

    DOCONCF_L0 [counts]: represents either voltage_counts or frequency, depending on DOFST series.
    DOCONCF_L2 [micro-mole/kg]:  oxygen concentration, corrected for temperature, salinity and pressure.                                   (a) calculated by do2_salinity_correction from DOCONCS_L1.
                                 calculated by do2_dofst_volt and do2_dofst_frequency
"""


def dosta_phase_volt_to_degree(phase_volt):
    """
    Description:
        Computes the DOCONCS-DEG_L0 data product from DOCONCS-VLT_L0, the analog output
        of a DOSTA Aanderaa Optode connected to a SBE CTD's 0-5 volt analog data channel.

    Usage:

        phase_degree = dosta_phase_volt_to_degree(phase_volt)

            where

        phase_degree = DOCONCS-DEG_L0 [degrees], calibrated phase measured by an
                       Aanderaa optode.
        phase_volt = DOCONCS-VLT_L0 [V], calibrated phase measured by an Aanderaa
                       optode that has been converted to volts for analog output.

    Implemented by:
        2015-08-04: Russell Desiderio. Initial Code.

    References:
        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)
    """
    # These coefficients to convert analog phase from volts to degrees are universal
    # for all Aanderaa optodes. Obtained from Shawn Sneddon at Xylem-Aanderaa.
    phase_degree = 10.0 + 12.0 * phase_volt
    return phase_degree


def dosta_Topt_volt_to_degC(T_optode_volt):
    """
    Description:
        Computes T_optode [degC], the DOSTA foil temperature as measured by its internal thermistor,
        from the analog output of a DOSTA Aanderaa Optode connected to a SBE CTD's 0-5 volt analog
        data channel.

    Usage:

        T_optode_degC = dosta_Topt_volt_to_degC(T_optode_volt)

            where

        T_optode_degC = optode foil temperature measured by an Aanderaa optode [degC].
        T_optode_volt = optode foil temperature measured by an Aanderaa optode
                        that has been converted to volts for analog output.

    Implemented by:
        2015-08-04: Russell Desiderio. Initial Code.

    Notes:

        T_optode is preferred for calculating oxygen concentation from DOSTA Aanderaa optodes
        because the permeability of the sensor's foil to oxygen is sensitive to temperature
        and this measurement is situated directly at the sensor's foil.

    References:
        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)
    """
    # These coefficients to convert analog T_optode from volts to degC are universal
    # for all Aanderaa optodes. Obtained from Shawn Sneddon at Xylem-Aanderaa.
    T_optode_degC = -5.0 + 8.0 * T_optode_volt
    return T_optode_degC


def o2_counts_to_uM(o2_counts):
    """
    Description:
        Computes the DOCONCS_L1 data product from a DOSTA Aanderaa Optode
        connected to a SBE 16+ V2 CTD via RS-232.

    Usage:

        DO = o2_counts_to_uM(o2_counts)

            where

        DO = DOCONCS_L1, dissolved oxygen concentration uncorrected for salinity
                         and pressure effects [micro-mole/L]
        o2_counts = DOCONCS-CNT_L0, counts from the CTD

    Implemented by:
        2013-04-26: Stuart Pearce. Initial Code.
        2015-04-10: Russell Desiderio. Added documentation and fill value handling.

    Notes:

        The DOCONCS_L1 data product has units of micromole/liter; SAF incorrectly
        lists the units for this L1 product as micromole/kg.

    References:
        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)
    """
    # replace fill values with nan
    o2_counts = replace_fill_with_nan(None, o2_counts)

    DO = (o2_counts / 10000.0) - 10.0
    return DO


def do2_SVU(calphase, temp, csv):
    """
    Description:

        Calculates the DOCONCS_L1 data product from autonomously-operated DOSTA
        (Aanderaa) instruments (on coastal surface moorings and gliders) using
        the Stern-Volmer-Uchida equation for calculating temperature corrected
        dissolved oxygen concentration.

    Usage:

        DO = do2_SVU(calphase, temp, csv)

            where

        DO = dissolved oxygen [micro-mole/L], DOCONCS_L1. see Notes.
        calphase = calibrated phase from an Oxygen sensor [deg], DOCONCS-DEG_L0
            (see DOCONCS DPS)
        temp = oxygen sensor foil temperature T(optode) [deg C], (see DOCONCS DPS)
        csv = Stern-Volmer-Uchida Calibration Coefficients array.
            7 element float array, (see DOCONCS DPS)

    Example:
        csv = np.array([0.002848, 0.000114, 1.51e-6, 70.42301, -0.10302,
                        -12.9462, 1.265377])
        calphase = 27.799
        temp = 19.841

        DO = do2_SVU(calphase, temp, csv)
        print DO
        > 363.900534505

    Implemented by:
        2013-04-26: Stuart Pearce. Initial Code.
        2015-04-10: Russell Desiderio. Revised code to work with CI implementation
                    of calibration coefficients: they are to be implemented as time-
                    vectorized arguments (tiled in the time dimension to match the
                    number of data packets). Fix for "blocker #2972".
        2015-08-04: Russell Desiderio. Added documentation.

    Notes:

        The DOCONCS_L1 data product has units of micromole/liter; SAF incorrectly
        lists the units for this L1 product as micromole/kg.

        The DOCONCS_L1 data product is uncorrected for salinity and pressure.

        The optode sensor's thermistor temperature should be used whenever possible
        because for the OOI DOSTAs (model 4831) it is situated directly at the sensor
        foil and the cal coefficients are derived in part to compensate for the change
        in oxygen permeability through the foil as a function of its temperature.

        The time constant of the model 4831 thermistor is < 2 seconds. Because the foil
        and therefore the calphase response time itself is 8 sec or 24 sec depending on
        the particular optode, there is little or no advantage to be gained by using a
        temperature sensor (eg, from a CTD) with a faster response. It is better to make
        sure that the temperature used most accurately reflects the foil temperature.

        On gliders, there is often a difference in CTD and optode temperature readings of
        1 degree Celsius, which translates to about a 5% difference in calculated oxygen
        concentration for a range of typical water column conditions.

    References:
        OOI (2012). Data Product Specification for Oxygen Concentration
            from "Stable" Instruments. Document Control Number
            1341-00520. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level
            >> 1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)
    """
    # this will work for both old and new CI implementations of cal coeffs.
    csv = np.atleast_2d(csv)

    # Calculate DO using Stern-Volmer:
    Ksv = csv[:, 0] + csv[:, 1]*temp + csv[:, 2]*(temp**2)
    P0 = csv[:, 3] + csv[:, 4]*temp
    Pc = csv[:, 5] + csv[:, 6]*calphase
    DO = ((P0/Pc) - 1) / Ksv
    return DO


def do2_salinity_correction(DO, P, T, SP, lat, lon, pref=0):
    """
    Description:

        Calculates the data product DOXYGEN_L2 (renamed from DOCONCS_L2) from DOSTA
        (Aanderaa) instruments by correcting the the DOCONCS_L1 data product for
        salinity and pressure effects and changing units.

    Usage:

        DOc = do2_salinity_correction(DO,P,T,SP,lat,lon, pref=0)

            where

        DOc = corrected dissolved oxygen [micro-mole/kg], DOXYGEN_L2
        DO = uncorrected dissolved oxygen [micro-mole/L], DOCONCS_L1
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

        DOc = do2_salinity_correction(DO,P,T,SP,lat,lon, pref=0)
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

        "Oxygen solubility in seawater: Better fitting equations", 1992,
        Garcia, H.E. and Gordon, L.I. Limnol. Oceanogr. 37(6) 1307-1312.
        Table 1, 5th column.
    """

    # density calculation from GSW toolbox
    SA = gsw.sa_from_sp(SP, P, lon, lat)
    CT = gsw.ct_from_t(SA, T, P)
    pdens = gsw.rho(SA, CT, pref)  # potential referenced to p=0

    # Convert from volume to mass units:
    DO = ne.evaluate('1000*DO/pdens')

    # Pressure correction:
    DO = ne.evaluate('(1 + (0.032*P)/1000) * DO')

    # Salinity correction (Garcia and Gordon, 1992, combined fit):
    S0 = 0
    ts = ne.evaluate('log((298.15-T)/(273.15+T))')
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

        DO = corrected dissolved oxygen [micro-mole/kg], DOCONCF_L2
        volt_counts = Oxygen sensor voltage counts [counts], DOCONCF_L0
        Voffset = Voltage offset [V].
        Soc = Oxygen signal slope [units are inherently inverse volts]
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

    Implemented by:
        2013-08-20: Stuart Pearce. Initial Code.
        2015-08-05: Russell Desiderio. Added fillvalue conversion to Nan.
    """
    # replace fill values with nan
    voltage_counts = replace_fill_with_nan(None, voltage_counts)

    # convert voltage counts to voltss
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

        DO = corrected dissolved oxygen [micro-mole/kg], DOCONCF_L2.
        frequency = Oxygen sensor frequency [Hz], DOCONCF_L0.
        Foffset = Frequency offset [Hz].
        Soc = Oxygen signal slope [units are inherently inverse Hz = seconds]
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

    Implemented by:
        2013-08-20: Stuart Pearce. Initial Code.
        2015-08-05: Russell Desiderio. Added fillvalue conversion to Nan.
    """
    # replace fill values with nan
    frequency = replace_fill_with_nan(None, frequency)

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
        do_raw = Oxygen sensor voltage [V] or frequency [Hz].
        offset = Voltage [V] or Frequency [Hz] offset.
        Soc = Oxygen signal slope [either inverse volts or seconds]
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

         "Oxygen solubility in seawater: Better fitting equations", 1992,
        Garcia, H.E. and Gordon, L.I. Limnol. Oceanogr. 37(6) 1307-1312.
        Table 1, 1st column.
   """
    # Get potential density using the TEOS-10 toolbox
    SA = gsw.sa_from_sp(SP, P, lon, lat)
    pot_rho_t = gsw.pot_rho_t_exact(SA, T, P, 0)

    # Oxygen saturation value using Garcia and Gordon (1992) fit to Benson and Krause data
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
