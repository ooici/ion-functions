#!/usr/bin/env python
"""
@package ion_functions.data.msp_functions
@file ion_functions/data/msp_functions.py
@author Craig Risien
@brief Module containing MASSP instrument related functions and wrapper functions


MASSP L1 Data Products

Data Product Specification for Dissolved Gas Concentrations (DISSGAS) from the
MASSP Instrument. Document Control Number 1341-00240.
https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >> Controlled
>> 1000 System Level >> 1341-00240_DPS_DISSGAS.pdf)

The OOI Level 1 Dissolved Gas Concentrations (DISSGAS) core data product is
produced by the MASSP instrument class. The data for the computation of this L1
core data product are derived from the Residual Gas Analyzer (RGA) integrated in
the MASSP instrument. The resulting L1 DISSGAS core data product is calculated
from the L0 Mass Spectral Intensities and the sample temperature, also measured
by the MASSP instrument, and is composed of the dissolved concentrations (uM) of
the individual gases: methane, ethane, hydrogen, argon, hydrogen sulfide, oxygen
and carbon dioxide. NOTE: For methane, the Nafion mode data is used, while for
the rest of the gasses, Direct mode data is used.

Function Name		  L1 DP Name		Description

calc_dissgas_smpmethcon   DISSGAS-SMPMETHCON	Dissolved Methane Conc. (uM) in Sample
calc_dissgas_smpethcon    DISSGAS-SMPETHNCON	Dissolved Ethane Conc. (uM) in Sample
calc_dissgas_smph2con     DISSGAS-SMPH2CON	Dissolved Hydrogen Conc. (uM) in Sample
calc_dissgas_smparcon     DISSGAS-SMPARCON	Dissolved Argon Conc. (uM) in Sample
calc_dissgas_smph2scon    DISSGAS-SMPH2SCON	Dissolved Hydrogen Sulfide Conc. (uM) in Sample
calc_dissgas_smpo2con     DISSGAS-SMPO2CON	Dissolved Oxygen Conc. (uM) in Sample
calc_dissgas_smpco2con    DISSGAS-SMPCO2CON	Dissolved Carbon Dioxide Conc. (uM) in Sample
calc_dissgas_bkgmethcon   DISSGAS-BKGMETHCON	Dissolved Methane Conc. (uM) in Background Water
calc_dissgas_bkgethcon    DISSGAS-BKGETHNCON	Dissolved Ethane Conc. (uM) in Background Water
calc_dissgas_bkgh2con     DISSGAS-BKGH2CON	Dissolved H2 Conc. (uM) in Background Water
calc_dissgas_bkgarcon     DISSGAS-BKGARCON	Dissolved AR Conc. (uM) in Background Water
calc_dissgas_bkgh2scon    DISSGAS-BKGH2SCON	Dissolved Hydrogen Sulfide Conc. (uM) in Background Water
calc_dissgas_bkgo2con     DISSGAS-BKGCO2CON 	Dissolved Carbon Dioxide Conc. (uM) in Background Water
calc_dissgas_bkgco2con    DISSGAS-BKGO2CON	Dissolved Oxygen Conc. (uM) in Background Water
calc_dissgas_cal1methcon  DISSGAS-CA1METHCON	Dissolved Methane Conc. (uM) in Calibration Solution 1
calc_dissgas_cal1co2con   DISSGAS-CA1CO2CON	Dissolved Carbon Dioxide Conc. (uM) in Calibration Solution 1
calc_dissgas_cal2methcon  DISSGAS-CA2METHCON	Dissolved Methane Conc. (uM) in Calibration Solution 2
calc_dissgas_cal2co2con   DISSGAS-CA2CO2CON	Dissolved Carbon Dioxide Conc. (uM) in Calibration Solution 2

...................................................................................

The auxiliary data product MASSP Calibration Range (CALRANG) is the collection
of parameters associated with the quality status for each gas concentration. A
value of 0 indicates that both the intensity and temperature used are within the
calibration range. A value of -1 indicates that the intensity used was below the
minimum of the calibration range. A value of 1 indicates that the intensity was
higher than the maximum of the calibration range, but the temperature was within
the calibration range. A value of 2 indicates that the intensity was within the
calibration range, but that the temperature was above the calibration range. A
value of 3 indicates that both the intensity and the temperature were above the
calibration range.

Function Name		  AUX L1 DP Name	Description

calc_calrang_smpmethcon   CALRANG-SMPMETH	Quality status for the Methane conc. in the sample water
calc_calrang_smpethcon    CALRANG-SMPETHN	Quality status for the Ethane conc. in the sample water
calc_calrang_smph2con     CALRANG-SMPH2		Quality status for the Hydrogen conc. in the sample water
calc_calrang_smparcon     CALRANG-SMPAR		Quality status for the Argon conc. in the sample water
calc_calrang_smph2scon    CALRANG-SMPH2S	Quality status for the H2S conc. in the sample water
calc_calrang_smpo2con     CALRANG-SMPO2		Quality status for the oxygen conc. in the sample water
calc_calrang_smpco2con    CALRANG-SMPCO2	Quality status for the CO2 conc. in the sample water
calc_calrang_bkgmethcon   CALRANG-BKGMETH	Quality status for the Methane conc. in the background water
calc_calrang_bkgethcon    CALRANG-BKGETHN	Quality status for the Ethane conc. in the background water
calc_calrang_bkgh2con     CALRANG-BKGH2		Quality status for the Hydrogen conc. in the background water
calc_calrang_bkgarcon     CALRANG-BKGAR		Quality status for the Argon conc. in the background water
calc_calrang_bkgh2scon    CALRANG-BKGH2S	Quality status for the H2S conc. in the background water
calc_calrang_bkgo2con     CALRANG-BKGO2		Quality status for the oxygen conc. in the background water
calc_calrang_bkgco2con    CALRANG-BKGCO2	Quality status for the CO2 conc. in the background water
calc_calrang_cal1methcon  CALRANG-CAL1METH	Quality status for the Methane conc. in the calibration fluid 1 water
calc_calrang_cal1co2con   CALRANG-CAL1CO2	Quality status for the CO2 conc. in the calibration fluid 1 water
calc_calrang_cal2methcon  CALRANG-CAL2METH	Quality status for the Methane conc. in the calibration fluid 2 water
calc_calrang_cal2co2con   CALRANG-CAL2CO2	Quality status for the CO2 conc. in the calibration fluid 2 water

...................................................................................

The auxiliary data product MASSP Time Stamp (TSTAMP) is the collection
of parameters associated with the time stamp for each gas concentration.

Function Name		    AUX L1 DP Name	Description

calc_timestamp_smpmethcon   TSTAMP-SMPMETH	Time stamp for the Methane conc. in the sample water
calc_timestamp_smpethcon    TSTAMP-SMPETHN	Time stamp for the Ethane conc. in the sample water
calc_timestamp_smph2con     TSTAMP-SMPH2	Time stamp for the Hydrogen conc. in the sample water
calc_timestamp_smparcon     TSTAMP-SMPAR	Time stamp for the Argon conc. in the sample water
calc_timestamp_smph2scon    TSTAMP-SMPH2S	Time stamp for the H2S conc. in the sample water
calc_timestamp_smpo2con     TSTAMP-SMPO2	Time stamp for the oxygen conc. in the sample water
calc_timestamp_smpco2con    TSTAMP-SMPCO2	Time stamp for the CO2 conc. in the sample water
calc_timestamp_bkgmethcon   TSTAMP-BKGMETH	Time stamp for the Methane conc. in the background water
calc_timestamp_bkgethcon    TSTAMP-BKGETHN	Time stamp for the Ethane conc. in the background water
calc_timestamp_bkgh2con     TSTAMP-BKGH2	Time stamp for the Hydrogen conc. in the background water
calc_timestamp_bkgarcon     TSTAMP-BKGAR	Time stamp for the Argon conc. in the background water
calc_timestamp_bkgh2scon    TSTAMP-BKGH2S	Time stamp for the H2S conc. in the background water
calc_timestamp_bkgo2con     TSTAMP-BKGO2	Time stamp for the oxygen conc. in the background water
calc_timestamp_bkgco2con    TSTAMP-BKGCO2	Time stamp for the CO2 conc. in the background water
calc_timestamp_cal1methcon  TSTAMP-CAL1METH	Time stamp for the Methane conc. in the calibration fluid 1 water
calc_timestamp_cal1co2con   TSTAMP-CAL1CO2	Time stamp for the CO2 conc. in the calibration fluid 1 water
calc_timestamp_cal2methcon  TSTAMP-CAL2METH	Time stamp for the Methane conc. in the calibration fluid 2 water
calc_timestamp_cal2co2con   TSTAMP-CAL2CO2	Time stamp for the CO2 conc. in the calibration fluid 2 water

...................................................................................

Functions that calculate all addtional L1 auxiliary data products are listed below.

The auxiliary data product MASSP Sample Inlet (MSINLET) is the collection of
parameters associated with measurement of sample properties at the time of gas
equilibration across the gas permeable membrane.

Function Name		            AUX L1 DP Name	Description

calc_msinlet_smpphint               MSINLET-SMPPHINT	Sample pH Intensity
calc_msinlet_smpphint_timestamp     TSTAMP-SMPPHINT	Time stamp of Sample pH Intensity
calc_msinlet_bkgphint               MSINLET-BKGPHINT	Background Water pH Intensity
calc_msinlet_bkgphint_timestamp     TSTAMP-BKGPHINT	Time stamp of Background Water pH Intensity
calc_msinlet_cal1phint              MSINLET-CA1PHINT	Calibration Solution 1 pH Intensity
calc_msinlet_cal1phint_timestamp    TSTAMP-CA1PHINT     Time stamp of Calibration Solution 1 pH Intensity
calc_msinlet_cal2phint              MSINLET-CA2PHINT	Calibration Solution 2 pH Intensity
calc_msinlet_cal1phint_timestamp    TSTAMP-CA2PHINT     Time stamp of Calibration Solution 2 pH Intensity
calc_smpnafeff                      NAFEFF  		Nafion Drier Efficiency
calc_smpnafeff_timestamp            TSTAMP-NAFEFF  	Time stamp of Nafion Drier Efficiency

...................................................................................
...................................................................................
...................................................................................

MASSP L2 Data Products

Data Product Specification for Dissolved Gas Concentrations (TOTLGAS) from the
MASSP Instrument. Document Control Number 1341-00XXX.
https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >> Controlled
>> 1000 System Level >> 1341-00XXX_DPS_TOTLGAS.pdf)

The OOI Level 2 Dissolved Gas Concentration (DISSGAS) core data product is
produced by the MASSP instrument class. The data for the computation of this L2
core data product are derived from the L1 core data product DISSGAS. The
resulting L2 TOTLGAS core data product is calculated from the individual
dissolved gas concentrations, the inlet fluid temperature, and the pH of the
fluid also measured by the MASSP instrument, and is composed of the total
concentrations (uM) of the individual gases: hydrogen sulfide and carbon
dioxide.

Function Name		    L2 DP Name		Description

calc_l2_totlgas_smph2scon   TOTLGAS-SMPH2SCON	Total Hydrogen Sulfide Conc. (uM) in Sample Water
calc_l2_totlgas_smpco2con   TOTLGAS-SMPCO2CON	Total Carbon Dioxide Conc. (uM) in Sample Water
calc_l2_totlgas_bkgh2scon   TOTLGAS-BKGH2SCON	Total Hydrogen Sulfide Conc. (uM) in Background Water
calc_l2_totlgas_bkgco2con   TOTLGAS-BKGCO2CON	Total Carbon Dioxide  Conc. (uM) in Background Water

...................................................................................

The auxiliary data product MASSP Time Stamp (TSTAMP) is the collection
of parameters associated with the time stamp for each gas concentration.

Function Name		            AUX L2 DP Name	Description

calc_timestamp_totlgas_smph2scon    TSTAMP-SMPH2SCON	Time stamp for the total H2S conc. in the sample water
calc_timestamp_totlgas_smpco2con    TSTAMP-SMPCO2CON	Time stamp for the total CO2 conc. in the sample water
calc_timestamp_totlgas_bkgh2scon    TSTAMP-BKGH2SCON	Time stamp for the total H2S conc. in the background water
calc_timestamp_totlgas_bkgco2con    TSTAMP-BKGCO2CON	Time stamp for the total CO2 conc. in the background water

...................................................................................

Functions that calculate all addtional L2 auxiliary data products are listed below:

The auxiliary data product MASSP Equilibrated Water (MSWATER) is the collection
of higher level products describing the pH state of the sampled and background
water at equilibration and measurement by the Residual Gas Analyzer (RGA),
onboard the MASSP instrument. These functions are required to calculate the above
L2 data products namely TOTLGAS-SMPH2SCON, TOTLGAS-SMPCO2CON, TOTLGAS-BKGH2SCON,
and TOTLGAS-BKGCO2CON.

Function Name		    AUX L2 DP Name	Description

calc_l2_mswater_smpphval    MSWATER-SMPPHVAL	Mass Spectrometer Equilibrated Sample Water pH Value
calc_l2_mswater_bkgphval    MSWATER-BKGPHVAL	Mass Spectrometer Equilibrated Background Water pH Value

...................................................................................
...................................................................................
...................................................................................

The core functions used by all of the wrapper functions described above are listed below. Note that
the mass-to-charge ratio is denoted as mz.

Function Name		        Description

SamplePreProcess                This subroutine takes in the SAMPLEINT array and produces intermediary
                                    variables sample-mz2, sample-mz18, sample-mz30, sample-mz32,
                                    sample-mz40, sample-mz44, sample-mz15 and sample-mz18Naf, sample-Tnaf,
                                    sample-Tdir as well as MSINLET-SMPPHINT AUX data products.

BackgroundPreProcess            This subroutine takes in the BKGNDINT array and produces intermediary
                                    variables bckgnd-mz2, bckgnd-mz30, bckgnd-mz32, bckgnd-mz40,
                                    bckgnd-mz44, bckgnd-mz15, bckgnd-Tnaf, bckgnd-Tdir, as well as
                                    MSINLET-BKGPHINT AUX data products.

Cal1PreProcess                  This subroutine takes in the DISSGAS-CALINT01 array and produces intermediary
                                    variables cal1-mz44, cal1-mz15, cal1-Tnaf, cal1-Tdir as well as
                                    MSINLET-CA1PHINT AUX data product.

Cal2PreProcess                  This subroutine takes in the DISSGAS-CALINT02 array and produces intermediary
                                    variables cal2-mz44, cal2-mz15, cal2-Tnaf, cal2-Tdir as well as
                                    MSINLET-CA2PHINT AUX data product.

gas_concentration               This sub-routine takes in a column range from DPS Table 1 (refered
                                    as to c1, c2, c3, c4), a corrected intensity (referred as x, from the
                                    Deconvolution subroutine), an averaged temperature (referred as T, see DPS Table 1),
                                    the pressure P of the sampling site, and calculate the final concentration used
                                    for L1 DISSGAS data products. This subroutine also assigns a value to the
                                    corresponding CALRANG parameter (see DPS Table 1) identifying the quality of the
                                    concentration value (indicate if it is out of calibration range for
                                    concentration and/or temperature). The subroutine also uses a temporary
                                    variable, tempCalRang used to compute the final value of CALRANG.

average_mz                      This subroutine takes in an mz as parameter (M), parameter w from the calibration
                                    table and a subset of n scans.

deconvolution_correction        This sub-routine takes in a main variable (see DPS Table 1), a second variable
                                    (see DPS Table 1), and a calibration lookup table (DPS Table 2) column
                                    range (see Table 1).

GasModeDetermination            Takes in the values of sample_valve1, sample_valve2, sample_valve3, and
                                    sample_valve4 and returns the GASMODE array.

SmpModeDetermination            Takes in the values of external_valve1_status, external_valve2_status,
                                    external_valve3_status, external_valve4_status, and external_valve5_status
                                    and returns the SMPMODE array.

"""

# import main python modules
import numpy as np


#Block of functions that calculate the L2 data products
def calc_l2_totlgas_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint,
                              gas_mode_sampleint, port_timestamp_sampleint_mcu,
                              ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    '''
    Below are the steps for calculating L2 TOTLGAS- SMPH2SCON and L2 TOTLGAS-
    BKGH2SCON core data products from L1 DISSGAS-SMPH2SCON and
    DISSGAS-BKGH2SCON, the L1 auxilliary data MSINLET-TEMP, the pressure at the
    site P and the value S from the calibration table, and the higher auxiliary
    data product MSWATER-SMPPHVAL and MSWATER-BKGPHVAL
    '''

    ph_temp_array = calc_l2_mswater_smpphval(port_timestamp_sampleint,
                                             L0_dissgas_sampleint,
                                             gas_mode_sampleint,
                                             port_timestamp_sampleint_mcu,
                                             ph_meter_sampleint_mcu,
                                             inlet_temp_sampleint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    smph2scon = calc_dissgas_smph2scon(port_timestamp_sampleint,
                                       L0_dissgas_sampleint,
                                       gas_mode_sampleint,
                                       port_timestamp_sampleint_mcu,
                                       ph_meter_sampleint_mcu,
                                       inlet_temp_sampleint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    k1T = 10**(-19.83 - (930.8 / (t+273.15)) + (2.8 * np.log(t + 273.15)) -
              (np.sqrt(PSU) * (-0.2391 + 35.685 / (t+273.15))) - (PSU * (0.0109 - (0.3776
               / (t+273.15)))))

    r = ((11.07 + 0.009 * t + 0.000942 * t**2) * 0.0689475729 * pressure +
        (-6.869 * 10**(-6) + 1.2835 * 10**(-7) * t) * pressure**2) / ((t+273.15) * 83.131)

    k1 = np.exp(r) * k1T

    beta = 1 + (k1 / 10**-ph)

    totlgas_smph2scon = beta * smph2scon

    return totlgas_smph2scon


def calc_l2_totlgas_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint,
                              gas_mode_sampleint, port_timestamp_sampleint_mcu,
                              ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    '''
    Below are the steps for calculating L2 TOTLGAS- SMPCO2CON and L2 TOTLGAS-
    BKGCO2CON core data products from L1 DISSGAS-SMPCO2CON and
    DISSGAS-BKGCO2CON, the L1 auxiliary data MSINLET-TEMP, the pressure at the
    site P and the value S from the calibration table and the higher auxiliary
    data MSWATER-SMPPHVAL and MSWATER-BKGPHVAL
    '''

    ph_temp_array = calc_l2_mswater_smpphval(port_timestamp_sampleint,
                                             L0_dissgas_sampleint,
                                             gas_mode_sampleint,
                                             port_timestamp_sampleint_mcu,
                                             ph_meter_sampleint_mcu,
                                             inlet_temp_sampleint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    smpco2con = calc_dissgas_smpco2con(port_timestamp_sampleint,
                                       L0_dissgas_sampleint,
                                       gas_mode_sampleint,
                                       port_timestamp_sampleint_mcu,
                                       ph_meter_sampleint_mcu,
                                       inlet_temp_sampleint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    K1T = np.exp((2.83655 - 2307.1266 / (t + 273.15) - 1.5529413 * np.log(t + 273.15) -
                  (0.20760841 * 4.0484 / (t + 273.15)) * np.sqrt(PSU) + 0.0846834 *
                  PSU - 0.00654208 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    K2T = np.exp((-9.226508 - 3351.616 / (t + 273.15) - 0.2005743 * np.log(t + 273.15) -
                  (0.106901773 * 23.9722 / (t + 273.15)) * np.sqrt(PSU) + 0.1130822 *
                  PSU - 0.00846934 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    r1 = (pressure * (1.758163 - 0.008763 * t - pressure * ((7.32 * 10**-6) -
         (2.0845 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    r2 = (pressure * (1.09075 + 0.00151 * t + pressure * ((2.69 * 10**-6) -
         (3.506 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    K1 = np.exp(r1) * K1T

    K2 = np.exp(r2) * K2T

    alpha = 1 + (K1 / (10**-ph)) + ((K1 * K2) / (10**-ph)**2)

    totlgas_smpco2con = alpha * smpco2con

    return totlgas_smpco2con


def calc_timestamp_totlgas_smph2scon(port_timestamp_sampleint,
                                     L0_dissgas_sampleint,
                                     gas_mode_sampleint,
                                     port_timestamp_sampleint_mcu,
                                     ph_meter_sampleint_mcu,
                                     inlet_temp_sampleint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved hydrogen
    sulfide in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint,
                                        gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_timestamp_totlgas_smpco2con(port_timestamp_sampleint,
                                     L0_dissgas_sampleint,
                                     gas_mode_sampleint,
                                     port_timestamp_sampleint_mcu,
                                     ph_meter_sampleint_mcu,
                                     inlet_temp_sampleint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved carbon dioxide
    in the sample water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint,
                                        gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_l2_totlgas_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                              gas_mode_bkgndint,
                              port_timestamp_bkgndint_mcu,
                              ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    '''
    Below are the steps for calculating L2 TOTLGAS- SMPH2SCON and L2 TOTLGAS-
    BKGH2SCON core data products from L1 DISSGAS-SMPH2SCON and
    DISSGAS-BKGH2SCON, the L1 auxilliary data MSINLET-TEMP, the pressure at the
    site P and the value S from the calibration table, and the higher auxiliary
    data product MSWATER-SMPPHVAL and MSWATER-BKGPHVAL
    '''

    ph_temp_array = calc_l2_mswater_bkgphval(port_timestamp_bkgndint,
                                             L0_dissgas_bkgndint,
                                             gas_mode_bkgndint,
                                             port_timestamp_bkgndint_mcu,
                                             ph_meter_bkgndint_mcu,
                                             inlet_temp_bkgndint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    bkgh2scon = calc_dissgas_bkgh2scon(port_timestamp_bkgndint,
                                       L0_dissgas_bkgndint,
                                       gas_mode_bkgndint,
                                       port_timestamp_bkgndint_mcu,
                                       ph_meter_bkgndint_mcu,
                                       inlet_temp_bkgndint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    k1T = 10**(-19.83 - (930.8 / (t+273.15)) + (2.8 * np.log(t + 273.15)) -
              (np.sqrt(PSU) * (-0.2391 + 35.685 / (t+273.15))) - (PSU * (0.0109 - (0.3776
               / (t+273.15)))))

    r = ((11.07 + 0.009 * t + 0.000942 * t**2) * 0.0689475729 * pressure +
        (-6.869 * 10**(-6) + 1.2835 * 10**(-7) * t) * pressure**2) / ((t+273.15) * 83.131)

    k1 = np.exp(r) * k1T

    beta = 1 + (k1 / 10**-ph)

    totlgas_bkgh2scon = beta * bkgh2scon

    return totlgas_bkgh2scon


def calc_l2_totlgas_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                              gas_mode_bkgndint,
                              port_timestamp_bkgndint_mcu,
                              ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    '''
    Below are the steps for calculating L2 TOTLGAS- SMPCO2CON and L2 TOTLGAS-
    BKGCO2CON core data products from L1 DISSGAS-SMPCO2CON and
    DISSGAS-BKGCO2CON, the L1 auxiliary data MSINLET-TEMP, the pressure at the
    site P and the value S from the calibration table and the higher auxiliary
    data MSWATER-SMPPHVAL and MSWATER-BKGPHVAL
    '''

    ph_temp_array = calc_l2_mswater_bkgphval(port_timestamp_bkgndint,
                                             L0_dissgas_bkgndint,
                                             gas_mode_bkgndint,
                                             port_timestamp_bkgndint_mcu,
                                             ph_meter_bkgndint_mcu,
                                             inlet_temp_bkgndint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    bkgco2con = calc_dissgas_bkgco2con(port_timestamp_bkgndint,
                                       L0_dissgas_bkgndint,
                                       gas_mode_bkgndint,
                                       port_timestamp_bkgndint_mcu,
                                       ph_meter_bkgndint_mcu,
                                       inlet_temp_bkgndint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    K1T = np.exp((2.83655 - 2307.1266 / (t + 273.15) - 1.5529413 * np.log(t + 273.15) -
                  (0.20760841 * 4.0484 / (t + 273.15)) * np.sqrt(PSU) + 0.0846834 *
                  PSU - 0.00654208 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    K2T = np.exp((-9.226508 - 3351.616 / (t + 273.15) - 0.2005743 * np.log(t + 273.15) -
                  (0.106901773 * 23.9722 / (t + 273.15)) * np.sqrt(PSU) + 0.1130822 *
                  PSU - 0.00846934 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    r1 = (pressure * (1.758163 - 0.008763 * t - pressure * ((7.32 * 10**-6) -
         (2.0845 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    r2 = (pressure * (1.09075 + 0.00151 * t + pressure * ((2.69 * 10**-6) -
         (3.506 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    K1 = np.exp(r1) * K1T

    K2 = np.exp(r2) * K2T

    alpha = 1 + (K1 / (10**-ph)) + ((K1 * K2) / (10**-ph)**2)

    totlgas_bkgco2con = alpha * bkgco2con

    return totlgas_bkgco2con


def calc_timestamp_totlgas_bkgh2scon(port_timestamp_bkgndint,
                                     L0_dissgas_bkgndint,
                                     gas_mode_bkgndint,
                                     port_timestamp_bkgndint_mcu,
                                     ph_meter_bkgndint_mcu,
                                     inlet_temp_bkgndint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved hydrogen sulfide
    in the background water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


def calc_timestamp_totlgas_bkgco2con(port_timestamp_bkgndint,
                                     L0_dissgas_bkgndint,
                                     gas_mode_bkgndint,
                                     port_timestamp_bkgndint_mcu,
                                     ph_meter_bkgndint_mcu,
                                     inlet_temp_bkgndint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved carbon dioxide
    in the background water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


#Block of wrapper functions for calculating the pH intensity auxiliary data products and associated timestamps

def calc_l2_mswater_smpphval(port_timestamp_sampleint, L0_dissgas_sampleint,
                             gas_mode_sampleint, port_timestamp_sampleint_mcu,
                             ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass,
                             massp_rga_steps_per_amu, calibration_table,
                             l2_ph_calibration_table):
    '''
    Below are the steps for processing the auxiliary products MSINLET-TEMP,
    MSINLET-SMPPHINT and MSINLET-BKGPHINT into the higher level auxiliary
    products MSWATER-SMPPHVAL MSWATER-BKGPHVAL.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #msinlet_temp is the 12th element of the preprocess_array array
    msinlet_temp = preprocess_array[11]
    #msinlet_smpphint is the 13th element of the preprocess_array array
    msinlet_smpphint = preprocess_array[12]

    A0 = l2_ph_calibration_table[0]
    A1 = l2_ph_calibration_table[1]
    A2 = l2_ph_calibration_table[2]
    a0 = l2_ph_calibration_table[4]
    a1 = l2_ph_calibration_table[3]
    a2 = l2_ph_calibration_table[5]

    pH = (A0 + (A1 * msinlet_temp) + (A2 * msinlet_temp**2)) * ((a2 * msinlet_smpphint**2) + (a1 * msinlet_smpphint) + a0 + 7)

    if pH < 2 or pH > 12:
        l2_msinlet_smpphint = -9999999.0
    else:
        l2_msinlet_smpphint = pH

    return l2_msinlet_smpphint, msinlet_temp


def calc_l2_mswater_bkgphval(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                             gas_mode_bkgndint, port_timestamp_bkgndint_mcu,
                             ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass,
                             massp_rga_steps_per_amu, calibration_table,
                             l2_ph_calibration_table):
    '''
    Below are the steps for processing the auxiliary products MSINLET-TEMP,
    MSINLET-SMPPHINT and MSINLET-BKGPHINT into the higher level auxiliary
    products MSWATER-SMPPHVAL MSWATER-BKGPHVAL.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #msinlet_temp is the 9th element of the preprocess_array array
    msinlet_temp = preprocess_array[8]
    #msinlet_bkgphint is the 10th element of the preprocess_array array
    msinlet_bkgphint = preprocess_array[9]

    A0 = l2_ph_calibration_table[0]
    A1 = l2_ph_calibration_table[1]
    A2 = l2_ph_calibration_table[2]
    a0 = l2_ph_calibration_table[4]
    a1 = l2_ph_calibration_table[3]
    a2 = l2_ph_calibration_table[5]

    pH = (A0 + (A1 * msinlet_temp) + (A2 * msinlet_temp**2)) * ((a2 * msinlet_bkgphint**2) + (a1 * msinlet_bkgphint) + a0 + 7)

    if pH < 2 or pH > 12:
        l2_msinlet_bkgphint = -9999999.0
    else:
        l2_msinlet_bkgphint = pH

    return l2_msinlet_bkgphint, msinlet_temp


def calc_msinlet_smpphint(port_timestamp_sampleint, L0_dissgas_sampleint,
                          gas_mode_sampleint, port_timestamp_sampleint_mcu,
                          ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass,
                          massp_rga_steps_per_amu, calibration_table):
    '''
    Sample pH intensity is output by a sensor onboard the MASSP instrument. It
    is the pH signal intensity of the Sample Water at the time of
    dissolved gas measurement.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint,
                                        gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #msinlet_smpphint is the 13th element of the preprocess_array array
    msinlet_smpphint = preprocess_array[12]

    return msinlet_smpphint


def calc_msinlet_smpphint_timestamp(port_timestamp_sampleint,
                                    L0_dissgas_sampleint, gas_mode_sampleint,
                                    port_timestamp_sampleint_mcu,
                                    ph_meter_sampleint_mcu,
                                    inlet_temp_sampleint_mcu,
                                    massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu, calibration_table):
    '''
    This is a wrapper function to calculate the MASSP timestamp while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_msinlet_bkgphint(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                          gas_mode_bkgndint, port_timestamp_bkgndint_mcu,
                          ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass,
                          massp_rga_steps_per_amu, calibration_table):
    '''
    Background Water pH intensity is output by a sensor onboard the MASSP
    instrument. It is the pH signal intensity of Background Water at
    the time of dissolved gas measurement
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #msinlet_bkgphint is the 10th element of the preprocess_array array
    msinlet_bkgphint = preprocess_array[9]

    return msinlet_bkgphint


def calc_msinlet_bkgphint_timestamp(port_timestamp_bkgndint,
                                    L0_dissgas_bkgndint, gas_mode_bkgndint,
                                    port_timestamp_bkgndint_mcu,
                                    ph_meter_bkgndint_mcu,
                                    inlet_temp_bkgndint_mcu,
                                    massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu, calibration_table):
    '''
    This is a wrapper function to calculate the MASSP timestamp while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #the nafion mode timestamp is the 11th element of the preprocess_array array
    bkg_nafion_timestamp = preprocess_array[10]

    return bkg_nafion_timestamp


def calc_msinlet_cal1phint(port_timestamp_calint01, L0_dissgas_calint01,
                           gas_mode_calint01, port_timestamp_calint01_mcu,
                           ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass,
                           massp_rga_steps_per_amu, calibration_table):
    '''
    Calibration Solution 1 pH intensity is output by a sensor onboard the MASSP
    instrument. It is the pH signal intensity of Calibration Solution 1 at the
    time of dissolved gas measurement
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01,
                                      L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu,
                                      ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table,
                                      calibration_table)

    #msinlet_cal1phint is the 5th element of the preprocess_array array
    msinlet_cal1phint = preprocess_array[4]

    return msinlet_cal1phint


def calc_msinlet_cal1phint_timestamp(port_timestamp_calint01,
                                     L0_dissgas_calint01,
                                     gas_mode_calint01,
                                     port_timestamp_calint01_mcu,
                                     ph_meter_calint01_mcu,
                                     inlet_temp_calint01_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    '''
    This is a wrapper function to calculate the MASSP timestamp while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01,
                                      L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu,
                                      ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table,
                                      calibration_table)

    #the nafion mode timestamp is the 6th element of the preprocess_array array
    cal1_nafion_timestamp = preprocess_array[5]

    return cal1_nafion_timestamp


def calc_msinlet_cal2phint(port_timestamp_calint02, L0_dissgas_calint02,
                           gas_mode_calint02, port_timestamp_calint02_mcu,
                           ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass,
                           massp_rga_steps_per_amu, calibration_table):
    '''
    Calibration Solution 2 pH intensity is output by a sensor onboard the MASSP
    instrument. It is the pH signal intensity of Calibration Solution 2 at the
    time of dissolved gas measurement
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #msinlet_cal2phint is the 5th element of the preprocess_array array
    msinlet_cal2phint = preprocess_array[4]

    return msinlet_cal2phint


def calc_msinlet_cal2phint_timestamp(port_timestamp_calint02,
                                     L0_dissgas_calint02,
                                     gas_mode_calint02,
                                     port_timestamp_calint02_mcu,
                                     ph_meter_calint02_mcu,
                                     inlet_temp_calint02_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    '''
    This is a wrapper function to calculate the MASSP timestamp while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02,
                                      L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu,
                                      ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table,
                                      calibration_table)

    #the direct mode timestamp is the 7th element of the preprocess_array array
    cal2_direct_timestamp = preprocess_array[6]

    return cal2_direct_timestamp


#Block of wrapper functions for calculating the nafion drier efficiency auxiliary data product and associated timestamp

def calc_smpnafeff(port_timestamp_sampleint, L0_dissgas_sampleint,
                   gas_mode_sampleint, port_timestamp_sampleint_mcu,
                   ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                   massp_rga_initial_mass, massp_rga_final_mass,
                   massp_rga_steps_per_amu, calibration_table):
    '''
    The auxiliary data product Nafion Drier Efficiency (NAFEFF) is an indicator
    of the drying efficiency of the nafion drier. The efficiency is represented
    as the percentage of water signal in nafion mode compared to direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #smpnafeff is the 10th element of the preprocess_array array
    smpnafeff = preprocess_array[9]

    return smpnafeff


def calc_smpnafeff_timestamp(port_timestamp_sampleint, L0_dissgas_sampleint,
                             gas_mode_sampleint, port_timestamp_sampleint_mcu,
                             ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass,
                             massp_rga_steps_per_amu, calibration_table):
    '''
    This is a wrapper function to calculate the MASSP timestamp while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the nafion mode timestamp is the 14th element of the preprocess_array array
    smp_nafion_timestamp = preprocess_array[13]

    return smp_nafion_timestamp


#Block of wrapper functions for calculating the L1 data products

def calc_dissgas_smpmethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    This is a wrapper function to calculate the in situ concentration (uM)
    of dissolved methane in the sample water as measured by the MASSP
    instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 11th element of the preprocess_array array
    average_temperature = preprocess_array[10]

    smpmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpmethcon = smpmethcon[0]

    return smpmethcon


def calc_dissgas_smpethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    This is a wrapper function to calculate the in situ concentration (uM)
    of dissolved ethane in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpethcon = smpethcon[0]

    return smpethcon


def calc_dissgas_smph2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    This is a wrapper function to calculate the in situ concentration (uM)
    of dissolved hydrogen in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2con = smph2con[0]

    return smph2con


def calc_dissgas_smparcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved argon in the
    sample water as measured by the MASSP instrument, while
    in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 8th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[7]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smparcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smparcon = smparcon[0]

    return smparcon


def calc_dissgas_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved hydrogen
    sulfide in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2scon = smph2scon[0]

    return smph2scon


def calc_dissgas_smpo2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved oxygen
    in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    #sample_mz34 is the 7th element of the preprocess_array array
    deconvolution_variable = preprocess_array[6]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpo2con = smpo2con[0]

    return smpo2con


def calc_dissgas_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved carbon dioxide
    in the sample water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 9th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[8]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpco2con = smpco2con[0]

    return smpco2con


def calc_dissgas_bkgmethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved methane in the background
    water as measured by the MASSP instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz15 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 8th element of the preprocess_array array
    average_temperature = preprocess_array[7]

    bkgmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgmethcon = bkgmethcon[0]

    return bkgmethcon


def calc_dissgas_bkgethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved ethane in the background
    water as measured by the MASSP instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgethcon = bkgethcon[0]

    return bkgethcon


def calc_dissgas_bkgh2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved hydrogen
    in the background water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 1st element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2con = bkgh2con[0]

    return bkgh2con


def calc_dissgas_bkgarcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved argon
    in the background water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgarcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgarcon = bkgarcon[0]

    return bkgarcon


def calc_dissgas_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved hydrogen sulfide
    in the background water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2scon = bkgh2scon[0]

    return bkgh2scon


def calc_dissgas_bkgo2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved oxygen
    in the background water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 4th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[3]
    #sample_mz34 is the 5th element of the preprocess_array array
    deconvolution_variable = preprocess_array[4]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgo2con = bkgo2con[0]

    return bkgo2con


def calc_dissgas_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved carbon dioxide
    in the background water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgco2con = bkgco2con[0]

    return bkgco2con


def calc_dissgas_cal1methcon(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                             port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved methane in the Calibration
    Solution 1 water as measured by the MASSP instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    cal1methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                    first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal1methcon = cal1methcon[0]

    return cal1methcon


def calc_dissgas_cal1co2con(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                            port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved carbon dioxide in the Calibration
    Solution 1 water as measured by the MASSP instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    cal1co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal1co2con = cal1co2con[0]

    return cal1co2con


def calc_dissgas_cal2methcon(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                             port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved methane in the Calibration
    Solution 2 water as measured by the MASSP instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    cal2methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                    first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal2methcon = cal2methcon[0]

    return cal2methcon


def calc_dissgas_cal2co2con(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                            port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    The in situ concentration (uM) of dissolved carbon dioxide in the Calibration
    Solution 2 water as measured by the MASSP instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    cal2co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal2co2con = cal2co2con[0]

    return cal2co2con


#Block of wrapper functions for calculating the timestamps of the L1 data products

def calc_timestamp_smpmethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                              port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    This is a wrapper function to calculate the timestamp of the in situ concentration (uM)
    of dissolved methane in the sample water as measured by the MASSP
    instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 14th element of the preprocess_array array
    smp_nafion_timestamp = preprocess_array[13]

    return smp_nafion_timestamp


def calc_timestamp_smpethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                             port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    This is a wrapper function to calculate the timestamp of the in situ concentration (uM)
    of dissolved ethane in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_timestamp_smph2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    This is a wrapper function to calculate the timestamp of the in situ concentration (uM)
    of dissolved hydrogen in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_timestamp_smparcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved argon in the
    sample water as measured by the MASSP instrument, while
    in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_timestamp_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                             port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved hydrogen
    sulfide in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_timestamp_smpo2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved oxygen
    in the sample water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_timestamp_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                             port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved carbon dioxide
    in the sample water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp


def calc_timestamp_bkgmethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                              port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved methane in the background
    water as measured by the MASSP instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 11th element of the preprocess_array array
    bkg_nafion_timestamp = preprocess_array[10]

    return bkg_nafion_timestamp


def calc_timestamp_bkgethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                             port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved ethane in the background
    water as measured by the MASSP instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


def calc_timestamp_bkgh2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved hydrogen
    in the background water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


def calc_timestamp_bkgarcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved argon
    in the background water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


def calc_timestamp_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                             port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved hydrogen sulfide
    in the background water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


def calc_timestamp_bkgo2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved oxygen
    in the background water as measured by the MASSP
    instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


def calc_timestamp_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                             port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved carbon dioxide
    in the background water as measured by the MASSP instrument,
    while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


def calc_timestamp_cal1methcon(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                               port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                               massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved methane in the Calibration
    Solution 1 water as measured by the MASSP instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 6th element of the preprocess_array array
    cal1_nafion_timestamp = preprocess_array[5]

    return cal1_nafion_timestamp


def calc_timestamp_cal1co2con(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                              port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved carbon dioxide in the Calibration
    Solution 1 water as measured by the MASSP instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 7th element of the preprocess_array array
    cal1_direct_timestamp = preprocess_array[6]

    return cal1_direct_timestamp


def calc_timestamp_cal2methcon(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                               port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                               massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved methane in the Calibration
    Solution 2 water as measured by the MASSP instrument, while in Nafion mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 6th element of the preprocess_array array
    cal2_nafion_timestamp = preprocess_array[5]

    return cal2_nafion_timestamp


def calc_timestamp_cal2co2con(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                              port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    '''
    The timestamp of the in situ concentration (uM) of dissolved carbon dioxide in the Calibration
    Solution 2 water as measured by the MASSP instrument, while in Direct mode.
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 7th element of the preprocess_array array
    cal2_direct_timestamp = preprocess_array[6]

    return cal2_direct_timestamp


#Block of wrapper functions for calculating the calibration ranges of the L1 data products

def calc_calrang_smpmethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    The auxiliary data product MASSP Calibration Range (CALRANG) is the
    collection of parameters associated with the quality status for each gas
    concentration. A value of 0 indicates that both the intensity and
    temperature used are within the calibration range. A value of -1 indicates
    that the intensity used was below the minimum of the calibration range. A
    value of 1 indicates that the intensity was higher than the maximum of the
    calibration range, but the temperature was within the calibration range. A
    value of 2 indicates that the intensity was within the calibration range,
    but that the temperature was above the calibration range. A value of 3
    indicates that both the intensity and the temperature were above the
    calibration range.

    Quality status for the Methane concentration in the sample water

    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 11th element of the preprocess_array array
    average_temperature = preprocess_array[10]

    smpmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpmethcon = smpmethcon[1]

    return smpmethcon


def calc_calrang_smpethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    Quality status for the Ethane concentration in the sample water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpethcon = smpethcon[1]

    return smpethcon


def calc_calrang_smph2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    Quality status for the Hydrogen concentration in the sample water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2con = smph2con[1]

    return smph2con


def calc_calrang_smparcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    Quality status for the Argon concentration in the sample water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 8th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[7]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smparcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smparcon = smparcon[1]

    return smparcon


def calc_calrang_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    Quality status for the H2S concentration in the sample water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2scon = smph2scon[1]

    return smph2scon


def calc_calrang_smpo2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    Quality status for the oxygen concentration in the sample water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    #sample_mz34 is the 7th element of the preprocess_array array
    deconvolution_variable = preprocess_array[6]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpo2con = smpo2con[1]

    return smpo2con


def calc_calrang_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    Quality status for the CO2 concentration in the sample water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 9th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[8]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpco2con = smpco2con[1]

    return smpco2con


def calc_calrang_bkgmethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    Quality status for the Methane concentration in the background water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz15 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 8th element of the preprocess_array array
    average_temperature = preprocess_array[7]

    bkgmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgmethcon = bkgmethcon[1]

    return bkgmethcon


def calc_calrang_bkgethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    Quality status for the Ethane concentration in the background water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgethcon = bkgethcon[1]

    return bkgethcon


def calc_calrang_bkgh2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    Quality status for the Hydrogen concentration in the background water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 1st element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2con = bkgh2con[1]

    return bkgh2con


def calc_calrang_bkgarcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    Quality status for the Argon concentration in the background water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgarcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgarcon = bkgarcon[1]

    return bkgarcon


def calc_calrang_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    Quality status for the H2S concentration in the background water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2scon = bkgh2scon[1]

    return bkgh2scon


def calc_calrang_bkgo2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    '''
    Quality status for the oxygen concentration in the background water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 4th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[3]
    #sample_mz34 is the 5th element of the preprocess_array array
    deconvolution_variable = preprocess_array[4]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgo2con = bkgo2con[1]

    return bkgo2con


def calc_calrang_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    '''
    Quality status for the CO2 concentration in the background water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgco2con = bkgco2con[1]

    return bkgco2con


def calc_calrang_cal1methcon(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                             port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    '''
    Quality status for the Methane concentration in the calibration fluid 1 water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    ca1methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca1methcon = ca1methcon[1]

    return ca1methcon


def calc_calrang_cal1co2con(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                            port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    Quality status for the CO2 concentration in the calibration fluid 1 water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    ca1co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca1co2con = ca1co2con[1]

    return ca1co2con


def calc_calrang_cal2methcon(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                             port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    '''
    Quality status for the Methane concentration in the calibration fluid 2 water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    ca2methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca2methcon = ca2methcon[1]

    return ca2methcon


def calc_calrang_cal2co2con(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                            port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    '''
    Quality status for the CO2 concentration in the calibration fluid 2 water
    '''

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    ca2co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca2co2con = ca2co2con[1]

    return ca2co2con


#Block of subfunctions called by the above wrapper functions that calculate the L1 and auxiliary data products

def gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                      first_column, last_column, sensor_depth, average_temperature):
    '''
    This sub-routine takes in a column. range from Table 1 (refered as to c1, c2, c3, c4),
    a corrected intensity (referred as x, from the Deconvolution subroutine), an averaged
    temperature (referred as T, see Table 1), the pressure P of the sampling site, and
    calculate the final concentration used for L1 DISSGAS data products. This subroutine
    also assigns a value to the corresponding CALRANG parameter (see Table 1) identifying
    the quality of the concentration value (indicate if it is out of calibration range
    for concentration and/or temperature). The subroutine also uses a temporary variable,
    tempCalRang used to compute the final value of CALRANG.

    The following vars will be hard coded for each of the wrapper functions
    deconvolution_variable, calibration_table, first_column, last_column

    '''

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #extract the four columns of the cal table that I need for a particular gas.
    calibration_table = calibration_table[:, first_column:last_column]

    #Check to see if one of the 4 calibration temperatures == the averaged inlet temperature
    ind = np.where(calibration_table[0, :] == average_temperature)[0]
    if np.size(ind) == 1:
        ct1 = np.where(calibration_table[0, :] == average_temperature)[0][0]
        tempCalRange = 1
    #Check to see if the averaged inlet temperature is greater than the highest calibration temperatures
    elif average_temperature >= calibration_table[0, 3]:
        ct1 = 3
        tempCalRange = 2
    #Otherwise figure out which two columns in the calibration table are needed.
    else:
        #find first column
        ct1 = np.where(calibration_table[0, :] < average_temperature)[0][-1]
        #find second column
        ct2 = np.where(calibration_table[0, :] > average_temperature)[0][0]
        concT2_flag = 1
        tempCalRange = 0

    corrected_intensity = deconvolution_correction(intermediate_mass_ratio, deconvolution_variable, calibration_table)

    #Check to see if the corrected intensity falls within the calibration values
    #Minimum values, row 4 in the cal table
    if corrected_intensity < calibration_table[3, ct1] or corrected_intensity < calibration_table[3, ct2]:
        calrange = -1
    #Maximum values, row 5 in the cal table
    elif corrected_intensity > calibration_table[4, ct1] or corrected_intensity > calibration_table[4, ct2]:
        calrange = tempCalRange + 1
    else:
        calrange = tempCalRange

    #P0 is row 6 (with row 0 being the first row) in the cal table
    if corrected_intensity < calibration_table[5, ct1]:
        alpha = calibration_table[6, ct1] + (calibration_table[7, ct1] * pressure) + (calibration_table[8, ct1] * pressure**2) + (calibration_table[9, ct1] * pressure**3)
        beta = calibration_table[14, ct1] + (calibration_table[15, ct1] * pressure) + (calibration_table[16, ct1] * pressure**2) + (calibration_table[17, ct1] * pressure**3)
        delta = calibration_table[22, ct1] + (calibration_table[23, ct1] * pressure) + (calibration_table[24, ct1] * pressure**2) + (calibration_table[25, ct1] * pressure**3)
        gamma = calibration_table[28, ct1] + (calibration_table[29, ct1] * pressure) + (calibration_table[30, ct1] * pressure**2) + (calibration_table[31, ct1] * pressure**3)
        zeta = calibration_table[36, ct1] + (calibration_table[37, ct1] * pressure) + (calibration_table[38, ct1] * pressure**2) + (calibration_table[39, ct1] * pressure**3)
    elif corrected_intensity >= calibration_table[5, ct1]:
        alpha = calibration_table[10, ct1] + (calibration_table[11, ct1] * pressure) + (calibration_table[12, ct1] * pressure**2) + (calibration_table[13, ct1] * pressure**3)
        beta = calibration_table[18, ct1] + (calibration_table[19, ct1] * pressure) + (calibration_table[20, ct1] * pressure**2) + (calibration_table[21, ct1] * pressure**3)
        delta = calibration_table[26, ct1] * np.exp(calibration_table[26, ct1] * pressure)
        gamma = calibration_table[32, ct1] + (calibration_table[33, ct1] * pressure) + (calibration_table[34, ct1] * pressure**2) + (calibration_table[35, ct1] * pressure**3)
        zeta = calibration_table[40, ct1] + (calibration_table[41, ct1] * pressure) + (calibration_table[42, ct1] * pressure**2) + (calibration_table[43, ct1] * pressure**3)

    #Calculate concT1
    concT1 = (alpha * corrected_intensity**2) + (beta * corrected_intensity) + (delta * np.exp(zeta * corrected_intensity)) + gamma

    if concT2_flag == 1:
        if corrected_intensity < calibration_table[5, ct2]:
            alpha = calibration_table[6, ct2] + (calibration_table[7, ct2] * pressure) + (calibration_table[8, ct2] * pressure**2) + (calibration_table[9, ct2] * pressure**3)
            beta = calibration_table[14, ct2] + (calibration_table[15, ct2] * pressure) + (calibration_table[16, ct2] * pressure**2) + (calibration_table[17, ct2] * pressure**3)
            delta = calibration_table[22, ct2] + (calibration_table[23, ct2] * pressure) + (calibration_table[24, ct2] * pressure**2) + (calibration_table[25, ct2] * pressure**3)
            gamma = calibration_table[28, ct2] + (calibration_table[29, ct2] * pressure) + (calibration_table[30, ct2] * pressure**2) + (calibration_table[31, ct2] * pressure**3)
            zeta = calibration_table[36, ct2] + (calibration_table[37, ct2] * pressure) + (calibration_table[38, ct2] * pressure**2) + (calibration_table[39, ct2] * pressure**3)
        elif corrected_intensity >= calibration_table[5, ct2]:
            alpha = calibration_table[10, ct2] + (calibration_table[11, ct2] * pressure) + (calibration_table[12, ct2] * pressure**2) + (calibration_table[13, ct2] * pressure**3)
            beta = calibration_table[18, ct2] + (calibration_table[19, ct2] * pressure) + (calibration_table[20, ct2] * pressure**2) + (calibration_table[21, ct2] * pressure**3)
            delta = calibration_table[26, ct2] * np.exp(calibration_table[26, ct2] * pressure)
            gamma = calibration_table[32, ct2] + (calibration_table[33, ct2] * pressure) + (calibration_table[34, ct2] * pressure**2) + (calibration_table[35, ct2] * pressure**3)
            zeta = calibration_table[40, ct2] + (calibration_table[41, ct2] * pressure) + (calibration_table[42, ct2] * pressure**2) + (calibration_table[43, ct2] * pressure**3)
        #Calculate concT2
        concT2 = (alpha * corrected_intensity**2) + (beta * corrected_intensity) + (delta * np.exp(zeta * corrected_intensity)) + gamma
        #Calculate concT
        concT = concT1 + ((concT2 - concT1) * (average_temperature - calibration_table[0, ct1])) / (calibration_table[0, ct2] - calibration_table[0, ct1])
    else:
        #Calculate concT
        concT = concT1

    if calrange == -1:
        final_conc = 0
    else:
        final_conc = calibration_table[44, ct1] * (concT - calibration_table[45, ct1])

    return final_conc, calrange


def average_mz(mz, data_in, mass_table, window):
    '''
    This subroutine takes in a mass-to-charge ratio mz, a subset of n scans
    and the mass_table and returns an intermediate mz (mass-to-charge) ratio.
    '''
    #find mz +/- window in the mass_table. The window value comes from the L1 Cal Table
    mz_ind = np.where((mass_table >= mz - window) & (mass_table <= mz + window))

    #subset the data_in array so that we are just dealing with the mz values
    #within the time period of interest
    temp_array = np.array(data_in[:, mz_ind])
    temp_array = np.squeeze(temp_array)

    #sort the array so that I can find the median of the three highest
    #values for each scan
    temp_array = np.sort(temp_array)

    #grab the median values
    median_array = temp_array[:, -2]
    #find and replace any negative values with zero
    median_ind = np.where(median_array < 0)
    median_array[median_ind] = 0
    #calculate the mean of the median values
    intermediate_mass_ratio = np.nanmean(median_array)

    return intermediate_mass_ratio


def deconvolution_correction(intermediate_mass_ratio, deconvolution_variable, calibration_table):
    '''
    This sub-routine takes in a main variable (intermediate_mass_ratio: see DPS Table 1), a
    second variable (deconvolution_variable: see DPS Table 1), and a calibration lookup
    table (DPS Table 2) and calculates a corrected intensity.
    '''
    #Equ 4 on page 13 of the DPS
    corrected_intensity = intermediate_mass_ratio - (calibration_table[2, 0] * deconvolution_variable) - calibration_table[1, 0]

    return corrected_intensity


def rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu):
    '''
    This subroutine takes in the values of rga_final_mass, rga_initial_mass, and
    rga_steps_per_amu, calculates the value for Tnb (Total number of values) and
    returns a table of the masses.
    '''

    Tnb = np.int(((massp_rga_final_mass - massp_rga_initial_mass) * massp_rga_steps_per_amu) + 1)

    mass_table = np.ones(Tnb)
    mass_table[0] = massp_rga_initial_mass

    for x in range(1, Tnb):
        mass_table[x] = mass_table[x-1] + (1 / np.float(massp_rga_steps_per_amu))

    mass_table = np.around(mass_table, decimals=1)

    return mass_table


def SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                     port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                     mass_table, calibration_table):
    '''
    This subroutine takes in L0 DISSGAS-SAMPLEINT and produces
    intermediary variables sample-mz2, sample-mz18, sample-mz30,
    sample-mz32, sample-mz40, sample-mz44, sample-mz15 and
    sample-mz18Naf, sample-Tnaf, sample-Tdir as well as
    MSINLET-SMPPHINT AUX data products. This subroutine groups
    the scans into two subsets corresponding to nafion and direct
    mode, and then extracts the mass data needed for L1 computations.

    Definitions:

    port_timestamp = time stamps associated with L0_dissgas_sampleint,
    gas_mode and sample_mode.

    L0_dissgas_sampleint = mass spectral data set for sample fluid
    (array, 10-16 ampere).

    gas_mode = The auxiliary data product Gas Measurement Mode
    (GASMODE) indicates the operating mode of the MASSP
    instrument and can have integer values of 0 and 1 for Direct
    and Nafion modes, respectively, or value of -1 if the instrument
    is in another operating mode.

    inlet_temp = Sample Temperature (oC) is output by a sensor onboard the MASSP instrument.
    It is the temperature of the Sample Water at the time of dissolved gas measurement.
    The value is set to -9999 when the instrument is not sampling.

    ph_meter = Sample pH intensity is output by a sensor onboard the MASSP instrument.
    It is the pH signal intensity (no unit) of the Sample Water at the time of dissolved
    gas measurement.

    mz = Mass to charge ratio

    sample-mz2 = Intensity returned by the averaging subroutine for mz 2 within DISSGAS-SAMPLEINT
    sample-mz15 = Intensity returned by the averaging subroutine for mz 15 within DISSGAS-SAMPLEINT
    sample-mz18Naf = Intensity returned by the averaging subroutine for mz 18 in nafion mode within DISSGAS-SAMPLEINT
    sample-mz18 = Intensity returned by the averaging subroutine for mz 18 within DISSGAS-SAMPLEINT
    sample-mz30 = Intensity returned by the averaging subroutine for mz 30 within DISSGAS-SAMPLEINT
    sample-mz32 = Intensity returned by the averaging subroutine for mz 32 within DISSGAS-SAMPLEINT
    sample-mz34 = Intensity returned by the averaging subroutine for mz 34 within DISSGAS-SAMPLEINT
    sample-mz40 = Intensity returned by the averaging subroutine for mz 40 within DISSGAS-SAMPLEINT
    sample-mz44 = Intensity returned by the averaging subroutine for mz 44 within DISSGAS-SAMPLEINT
    sample-Tdir = Averaged temperature in Sample fluid direct mode
    sample-Tnaf = Averaged temperature in Sample fluid nafion mode

    nafeff = The auxiliary data product Nafion Drier Efficiency (NAFEFF)
    is an indicator of the drying efficiency of the nafion drier. The
    efficiency is represented as the percentage of water signal in nafion
    mode compared to direct mode.

    '''

    #replace bad data with nans
    inlet_temp_sampleint_mcu[inlet_temp_sampleint_mcu == -127] = np.nan
    inlet_temp_sampleint_mcu[inlet_temp_sampleint_mcu == 85] = np.nan

    #find gas_mode_sampleint == 0 (direct mode)
    ind_direct = np.where(gas_mode_sampleint == 0)[0]
    Tchange = port_timestamp_sampleint_mcu[ind_direct[0]]
    Tlast = port_timestamp_sampleint_mcu[ind_direct[-1]]

    #find gas_mode_sampleint == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_sampleint == 1)[0]
    TlastScanNafion = port_timestamp_sampleint_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 180
    idx = (np.abs(port_timestamp_sampleint-(TlastScanNafion - 180))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_sampleint >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_sampleint[nafion_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanNafion. This subsetted array then gets fed into the ave. routine.
    nafion_samples = nafion_samples[:-1, ]

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.squeeze(np.array(port_timestamp_sampleint[nafion_samples_ind]))
    #DPS says to exclude the last scan at TlastScanNafion.
    nafion_mode_timestamp = np.around(np.nanmean(nafion_mode_timestamp[:-1, ]))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    sample_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    mass_charge_ratio = 18
    #not sure that this window of 0.5 is OK but the 18mz window is not specified in the cal table
    window = round(calibration_table[-1, 8], 1)
    sample_mz18naf = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period
    nafion_samples_ind = np.squeeze(np.where((port_timestamp_sampleint_mcu >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint_mcu <= TlastScanNafion)))
    sample_Tnaf = np.nanmean(inlet_temp_sampleint_mcu[nafion_samples_ind[:-1]])

    #ID timestamp closest to Tlast - 180
    idx = (np.abs(port_timestamp_sampleint-(Tlast - 180))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_sampleint >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_sampleint[direct_samples_ind, :]))

    #calculate direct mode timestamp
    direct_mode_timestamp = np.array(port_timestamp_sampleint[direct_samples_ind])
    direct_mode_timestamp = np.around(np.nanmean(np.squeeze(direct_mode_timestamp)))

    mass_charge_ratio = 2
    window = round(calibration_table[-1, 8], 1)
    sample_mz2 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 18
    #not sure that this window is true but the window is not specified in the cal table
    window = round(calibration_table[-1, 8], 1)
    sample_mz18 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 30
    window = round(calibration_table[-1, 4], 1)
    sample_mz30 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 32
    window = round(calibration_table[-1, 20], 1)
    sample_mz32 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 34
    window = round(calibration_table[-1, 16], 1)
    sample_mz34 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 40
    window = round(calibration_table[-1, 12], 1)
    sample_mz40 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    sample_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it sample-Tdir
    direct_samples_ind = np.where((port_timestamp_sampleint_mcu >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint_mcu <= Tlast))
    sample_Tdir = np.nanmean(inlet_temp_sampleint_mcu[direct_samples_ind])

    #average ph_meter_value for time period Tlast-1min:Tlast, call it msinlet_smpphint
    direct_samples_ind = np.where((port_timestamp_sampleint_mcu >= Tlast - 60) & (port_timestamp_sampleint_mcu <= Tlast))
    msinlet_smpphint = np.absolute(np.nanmean(ph_meter_sampleint_mcu[direct_samples_ind]))

    #Calculate NAFEFF, which is an indicator of the drying efficiency of the nafion drier
    nafeff = np.int(100 * (sample_mz18naf / sample_mz18))

    return (sample_mz15, sample_mz18naf, sample_mz2, sample_mz18, sample_mz30,
            sample_mz32, sample_mz34, sample_mz40, sample_mz44, nafeff, sample_Tnaf,
            sample_Tdir, msinlet_smpphint, nafion_mode_timestamp, direct_mode_timestamp)


def BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                         port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                         mass_table, calibration_table):
    '''
    This subroutine takes in L0 BKGNDINT and produces intermediary
    variables bckgnd-mz2, bckgnd-mz30, bckgnd -mz32, bckgnd -mz40,
    bckgnd -mz44, bckgnd -mz15, bckgnd-Tnaf, bckgnd-Tdir, as well
    as MSINLET-BKGPHINT AUX data products. This subroutine groups
    the scans into two subsets corresponding to nafion and direct
    mode, and then extracts the mass data needed for L1 computations.

    '''

    #replace bad data with nans
    inlet_temp_bkgndint_mcu[inlet_temp_bkgndint_mcu == -127] = np.nan
    inlet_temp_bkgndint_mcu[inlet_temp_bkgndint_mcu == 85] = np.nan

    #find gas_mode_bkgndint == 0 (direct mode)
    ind_direct = np.where(gas_mode_bkgndint == 0)[0]
    Tchange = port_timestamp_bkgndint_mcu[ind_direct[0]]
    Tlast = port_timestamp_bkgndint_mcu[ind_direct[-1]]

    #ID timestamp closest to Tlast - 180
    idx = (np.abs(port_timestamp_bkgndint-(Tlast - 180))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_bkgndint >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_bkgndint[direct_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanDirect. This subsetted array then gets fed into the ave. routine.
    direct_samples = direct_samples[:-1, ]

    #calculate direct mode timestamp
    direct_mode_timestamp = np.squeeze(np.array(port_timestamp_bkgndint[direct_samples_ind]))
    #DPS says to exclude the last scan at TlastScanDirect.
    direct_mode_timestamp = np.around(np.nanmean(direct_mode_timestamp[:-1, ]))

    mass_charge_ratio = 2
    window = round(calibration_table[-1, 8], 1)
    bckgnd_mz2 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 30
    window = round(calibration_table[-1, 4], 1)
    bckgnd_mz30 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 32
    window = round(calibration_table[-1, 20], 1)
    bckgnd_mz32 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 34
    window = round(calibration_table[-1, 16], 1)
    bckgnd_mz34 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 40
    window = round(calibration_table[-1, 12], 1)
    bckgnd_mz40 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    bckgnd_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it bckgnd-Tdir
    direct_samples_ind = np.squeeze(np.where((port_timestamp_bkgndint_mcu >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint_mcu <= Tlast)))
    bckgnd_Tdir = np.nanmean(inlet_temp_bkgndint_mcu[direct_samples_ind[:-1]])

    #find gas_mode_bkgndint == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_bkgndint == 1)[0]
    TlastScanNafion = port_timestamp_bkgndint_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 180
    idx = (np.abs(port_timestamp_bkgndint-(TlastScanNafion - 180))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_bkgndint >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_bkgndint[nafion_samples_ind, :]))

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.array(port_timestamp_bkgndint[nafion_samples_ind])
    nafion_mode_timestamp = np.around(np.nanmean(np.squeeze(nafion_mode_timestamp)))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    bckgnd_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period here, call it bckgnd-Tnaf
    nafion_samples_ind = np.where((port_timestamp_bkgndint_mcu >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint_mcu <= TlastScanNafion))
    bckgnd_Tnaf = np.nanmean(inlet_temp_bkgndint_mcu[nafion_samples_ind])

    #average ph_meter_value for time period TlastScanNafion-1min:TlastScanNafion, call it msinlet_bkgphint
    nafion_samples_ind = np.where((port_timestamp_bkgndint_mcu >= TlastScanNafion - 60) & (port_timestamp_bkgndint_mcu <= TlastScanNafion))
    msinlet_bkgphint = np.absolute(np.nanmean(ph_meter_bkgndint_mcu[nafion_samples_ind]))

    return (bckgnd_mz2, bckgnd_mz15, bckgnd_mz30, bckgnd_mz32, bckgnd_mz34,
            bckgnd_mz40, bckgnd_mz44, bckgnd_Tnaf, bckgnd_Tdir, msinlet_bkgphint,
            nafion_mode_timestamp, direct_mode_timestamp)


def Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                   port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                   inlet_temp_calint01_mcu, mass_table, calibration_table):
    '''
    This subroutine takes in L0 DISSGAS-CALINT01 and produces
    intermediary variables cal1-mz44, cal1-mz15, cal1-Tnaf,
    cal1-Tdir as well as MSINLET-CA1PHINT AUX data product.
    This subroutine groups the scans into two subsets
    corresponding to nafion and direct mode, and then extracts
    the mass data needed for L1 computations. This subroutine
    is very similar to the BackgroundPreProcess subroutine,
    with just different intermediary variable assignations.

    '''

    #replace bad data with nans
    inlet_temp_calint01_mcu[inlet_temp_calint01_mcu == -127] = np.nan
    inlet_temp_calint01_mcu[inlet_temp_calint01_mcu == 85] = np.nan

    #find gas_mode_calint01 == 0 (direct mode)
    ind_direct = np.where(gas_mode_calint01 == 0)[0]
    Tchange = port_timestamp_calint01_mcu[ind_direct[0]]
    Tlast = port_timestamp_calint01_mcu[ind_direct[-1]]

    #ID timestamp closest to Tlast - 60
    idx = (np.abs(port_timestamp_calint01-(Tlast - 60))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_calint01 >= port_timestamp_calint01[idx]) & (port_timestamp_calint01 <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_calint01[direct_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanDirect. This subsetted array then gets fed into the ave. routine.
    direct_samples = direct_samples[:-1, ]

    #calculate direct mode timestamp
    direct_mode_timestamp = np.squeeze(np.array(port_timestamp_calint01[direct_samples_ind]))
    #DPS says to exclude the last scan at TlastScanDirect.
    direct_mode_timestamp = np.around(np.nanmean(direct_mode_timestamp[:-1, ]))

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    cal1_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it cal1-Tdir
    direct_samples_ind = np.squeeze(np.where((port_timestamp_calint01_mcu >= port_timestamp_calint01[idx]) & (port_timestamp_calint01_mcu <= Tlast)))
    cal1_Tdir = np.nanmean(inlet_temp_calint01_mcu[direct_samples_ind[:-1]])

    #find gas_mode_calint01 == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_calint01 == 1)[0]
    TlastScanNafion = port_timestamp_calint01_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 60
    idx = (np.abs(port_timestamp_calint01-(TlastScanNafion - 60))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_calint01 >= port_timestamp_calint01[idx]) & (port_timestamp_calint01 <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_calint01[nafion_samples_ind, :]))

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.array(port_timestamp_calint01[nafion_samples_ind])
    nafion_mode_timestamp = np.around(np.nanmean(np.squeeze(nafion_mode_timestamp)))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    cal1_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period here, call it cal1-Tnaf
    nafion_samples_ind = np.where((port_timestamp_calint01_mcu >= port_timestamp_calint01[idx]) & (port_timestamp_calint01_mcu <= TlastScanNafion))
    cal1_Tnaf = np.nanmean(inlet_temp_calint01_mcu[nafion_samples_ind])

    #average ph_meter_value for time period TlastScanNafion-1min:TlastScanNafion, call it msinlet_cal1phint
    nafion_samples_ind = np.where((port_timestamp_calint01_mcu >= TlastScanNafion - 60) & (port_timestamp_calint01_mcu <= TlastScanNafion))
    msinlet_cal1phint = np.absolute(np.nanmean(ph_meter_calint01_mcu[nafion_samples_ind]))

    return (cal1_mz15, cal1_mz44, cal1_Tnaf, cal1_Tdir, msinlet_cal1phint,
            nafion_mode_timestamp, direct_mode_timestamp)


def Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                   port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                   inlet_temp_calint02_mcu, mass_table, calibration_table):
    '''
    This subroutine takes in L0 DISSGAS-CALINT02 and produces
    intermediary variables cal2-mz44, cal2-mz15, cal2-Tnaf,
    cal2-Tdir as well as MSINLET-CA2PHINT AUX data product.
    This subroutine groups the scans into two subsets
    corresponding to nafion and direct mode, and then extracts
    the mass data needed for L1 computations. This subroutine
    is very similar to the SamplePreProcess subroutine, with
    just different intermediary variable assignations.

    '''

    #replace bad data with nans
    inlet_temp_calint02_mcu[inlet_temp_calint02_mcu == -127] = np.nan
    inlet_temp_calint02_mcu[inlet_temp_calint02_mcu == 85] = np.nan

    #find gas_mode_calint02 == 0 (direct mode)
    ind_direct = np.where(gas_mode_calint02 == 0)[0]
    Tchange = port_timestamp_calint02_mcu[ind_direct[0]]
    Tlast = port_timestamp_calint02_mcu[ind_direct[-1]]

    #find gas_mode_calint02 == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_calint02 == 1)[0]
    TlastScanNafion = port_timestamp_calint02_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 60
    idx = (np.abs(port_timestamp_calint02-(TlastScanNafion - 60))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_calint02 >= port_timestamp_calint02[idx]) & (port_timestamp_calint02 <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_calint02[nafion_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanNafion. This subsetted array then gets fed into the ave. routine.
    nafion_samples = nafion_samples[:-1, ]

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.squeeze(np.array(port_timestamp_calint02[nafion_samples_ind]))
    #DPS says to exclude the last scan at TlastScanNafion.
    nafion_mode_timestamp = np.around(np.nanmean(nafion_mode_timestamp[:-1, ]))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    cal2_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period
    nafion_samples_ind = np.squeeze(np.where((port_timestamp_calint02_mcu >= port_timestamp_calint02[idx]) & (port_timestamp_calint02_mcu <= TlastScanNafion)))
    cal2_Tnaf = np.nanmean(inlet_temp_calint02_mcu[nafion_samples_ind[:-1]])

    #ID timestamp closest to Tlast - 60
    idx = (np.abs(port_timestamp_calint02-(Tlast - 60))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_calint02 >= port_timestamp_calint02[idx]) & (port_timestamp_calint02 <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_calint02[direct_samples_ind, :]))

    #calculate direct mode timestamp
    direct_mode_timestamp = np.array(port_timestamp_calint02[direct_samples_ind])
    direct_mode_timestamp = np.around(np.nanmean(np.squeeze(direct_mode_timestamp)))

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    cal2_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it cal2-Tdir
    direct_samples_ind = np.where((port_timestamp_calint02_mcu >= port_timestamp_calint02[idx]) & (port_timestamp_calint02_mcu <= Tlast))
    cal2_Tdir = np.nanmean(inlet_temp_calint02_mcu[direct_samples_ind])

    #average ph_meter_value for time period Tlast-1min:Tlast, call it msinlet_cal2phint
    direct_samples_ind = np.where((port_timestamp_calint02_mcu >= Tlast - 60) & (port_timestamp_calint02_mcu <= Tlast))
    msinlet_cal2phint = np.absolute(np.nanmean(ph_meter_calint02_mcu[direct_samples_ind]))
    #associate msinlet_smpphint with cal2_mz44 time stamp

    return (cal2_mz15, cal2_mz44, cal2_Tnaf, cal2_Tdir, msinlet_cal2phint,
            nafion_mode_timestamp, direct_mode_timestamp)


def GasModeDetermination(sample_valve1, sample_valve2, sample_valve3, sample_valve4):
    '''
    This subroutine takes in the values of sample_valve1, sample_valve2,
    sample_valve3, and sample_valve4 and returns the value for the AUX GASMODE
    data product.
    '''

    data_array_size = np.shape(sample_valve1)
    gasmode_array = np.ones(data_array_size[0])
    gasmode_array[gasmode_array == 1] = np.nan

    ind = np.where(sample_valve4 == 1)
    gasmode_array[ind] = -1
    ind = np.where((sample_valve2 == 1) & (sample_valve1 == 0))
    gasmode_array[ind] = 0
    ind = np.where((sample_valve1 == 1) & (sample_valve2 == 0) & (sample_valve3 == 0))
    gasmode_array[ind] = 1

    return gasmode_array


def SmpModeDetermination(external_valve1_status, external_valve2_status,
                         external_valve3_status, external_valve4_status,
                         external_valve5_status):
    '''
    This subroutine takes in the values of external_valve1_status,
    external_valve2_status, external_valve3_status, external_valve4_status, and
    external_valve5_status and returns the value for the AUX SMPMODE data
    product.
    '''

    data_array_size = np.shape(external_valve1_status)
    smpmode_array = np.ones(data_array_size[0])
    smpmode_array[smpmode_array == 1] = np.nan

    ind = np.where((external_valve1_status == 0) & (external_valve2_status == 0) &
                  (external_valve3_status == 0) & (external_valve4_status == 0) &
                  (external_valve5_status == 0))
    smpmode_array[ind] = 2

    ind = np.where((external_valve1_status == 1) & (external_valve2_status == 0) &
                  (external_valve3_status == 0) & (external_valve4_status == 0) &
                  (external_valve5_status == 1))
    smpmode_array[ind] = 1

    ind = np.where((external_valve1_status == 1) & (external_valve2_status == 1) &
                  (external_valve3_status == 0) & (external_valve4_status == 0) &
                  (external_valve5_status == 0))
    smpmode_array[ind] = -1

    ind = np.where((external_valve1_status == 1) & (external_valve2_status == 0) &
                  (external_valve3_status == 1) & (external_valve4_status == 0) &
                  (external_valve5_status == 0))
    smpmode_array[ind] = -2

    return smpmode_array
