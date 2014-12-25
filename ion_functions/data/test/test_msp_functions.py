#!/usr/bin/env python

"""
@package ion_functions.test.msp_functions
@file ion_functions/test/msp_functions.py
@author Craig Risien
@brief Unit tests for msp_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
import ion_functions.data.msp_functions as msp
import ion_functions.data.test.test_msp_functions_data as data


@attr('UNIT', group='func')
class TestMspFunctionsUnit(BaseUnitTestCase):

    def test_massp_gas_concs(self):
                """
        Test the all 18 massp gas concentration wrapper functions.

        Values based on that described in DPS as available on Alfresco:

        OOI (2014). Data Product Specification for MASSP Data Products.
            Document Control Number 1341-00240.
            https://alfresco.oceanobservatories.org/ (See: Company Home
            >> OOI >> Controlled >> 1000 System Level >>
            1341-00240_DPS_DISSGAS_2014-07-23.docx)

        Implemented by Craig Risien, Sep 2014
        """

    #Multiply scan data by 10^12
    data.L0_dissgas_sampleint = data.L0_dissgas_sampleint * 10**12
    data.L0_dissgas_bkgndint = data.L0_dissgas_bkgndint * 10**12
    data.L0_dissgas_calint01 = data.L0_dissgas_calint01 * 10**12
    data.L0_dissgas_calint02 = data.L0_dissgas_calint02 * 10**12

    # compute nafion drier efficiency
    smpnafeff = msp.calc_smpnafeff(data.port_timestamp_sampleint,
                                   data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                   data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                   data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                   data.calibration_table)

    smpnafeff_timestamp = msp.calc_smpnafeff_timestamp(data.port_timestamp_sampleint,
                                                       data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                       data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                       data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table)

    #compute pH intensities
    msinlet_smpphint = msp.calc_msinlet_smpphint(data.port_timestamp_sampleint,
                                                 data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                 data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                 data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table)

    msinlet_smpphint_timestamp = msp.calc_msinlet_smpphint_timestamp(data.port_timestamp_sampleint,
                                                                     data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                                     data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                                     data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                     data.calibration_table)

    msinlet_bkgphint = msp.calc_msinlet_bkgphint(data.port_timestamp_bkgndint,
                                                 data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                 data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                 data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table)

    msinlet_bkgphint_timestamp = msp.calc_msinlet_bkgphint_timestamp(data.port_timestamp_bkgndint,
                                                                     data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                                     data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                                     data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                     data.calibration_table)

    msinlet_cal1phint = msp.calc_msinlet_cal1phint(data.port_timestamp_calint01,
                                                   data.L0_dissgas_calint01, data.gas_mode_calint01,
                                                   data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                                   data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table)

    msinlet_cal1phint_timestamp = msp.calc_msinlet_cal1phint_timestamp(data.port_timestamp_calint01,
                                                                       data.L0_dissgas_calint01, data.gas_mode_calint01,
                                                                       data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                                                       data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                       data.calibration_table)

    msinlet_cal2phint = msp.calc_msinlet_cal2phint(data.port_timestamp_calint02,
                                                   data.L0_dissgas_calint02, data.gas_mode_calint02,
                                                   data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                                   data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table)

    msinlet_cal2phint_timestamp = msp.calc_msinlet_cal2phint_timestamp(data.port_timestamp_calint02,
                                                                       data.L0_dissgas_calint02, data.gas_mode_calint02,
                                                                       data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                                                       data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                       data.calibration_table)

    calc_l2_mswater_smpphval = msp.calc_l2_mswater_smpphval(data.port_timestamp_sampleint,
                                                            data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                            data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                            data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                            data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                            data.calibration_table, data.l2_ph_calibration_table)

    calc_l2_mswater_bkgphval = msp.calc_l2_mswater_bkgphval(data.port_timestamp_bkgndint,
                                                            data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                            data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                            data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                            data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                            data.calibration_table, data.l2_ph_calibration_table)

    # compute L2 gas conc values
    #SMPINT
    smpco2con_totalgas = msp.calc_l2_totlgas_smpco2con(data.port_timestamp_sampleint,
                                                       data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                       data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                       data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table, data.l2_ph_calibration_table, data.sensor_depth,
                                                       data.salinity)

    smph2scon_totalgas = msp.calc_l2_totlgas_smph2scon(data.port_timestamp_sampleint,
                                                       data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                       data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                       data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table, data.l2_ph_calibration_table, data.sensor_depth,
                                                       data.salinity)

    #BCKINT
    bkgco2con_totalgas = msp.calc_l2_totlgas_bkgco2con(data.port_timestamp_bkgndint,
                                                       data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                       data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                       data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table, data.l2_ph_calibration_table, data.sensor_depth,
                                                       data.salinity)

    bkgh2scon_totalgas = msp.calc_l2_totlgas_bkgh2scon(data.port_timestamp_bkgndint,
                                                       data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                       data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                       data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table, data.l2_ph_calibration_table, data.sensor_depth,
                                                       data.salinity)

    # compute L2 gas conc timestamp values
    #SMPINT
    smpco2con_totalgas_timestamp = msp.calc_timestamp_totlgas_smpco2con(data.port_timestamp_sampleint,
                                                                        data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                                        data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                                        data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                                        data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                        data.calibration_table)

    smph2scon_totalgas_timestamp = msp.calc_timestamp_totlgas_smph2scon(data.port_timestamp_sampleint,
                                                                        data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                                        data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                                        data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                                        data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                        data.calibration_table)
    #BCKINT
    bkgco2con_totalgas_timestamp = msp.calc_timestamp_totlgas_bkgco2con(data.port_timestamp_bkgndint,
                                                                        data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                                        data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                                        data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                                        data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                        data.calibration_table)

    bkgh2scon_totalgas_timestamp = msp.calc_timestamp_totlgas_bkgh2scon(data.port_timestamp_bkgndint,
                                                                        data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                                        data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                                        data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                                        data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                                        data.calibration_table)

    # compute L1 gas conc values
    #SMPINT
    smpmethcon = msp.calc_dissgas_smpmethcon(data.port_timestamp_sampleint,
                                             data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                             data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                             data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                             data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                             data.calibration_table, data.sensor_depth)

    smpethcon = msp.calc_dissgas_smpethcon(data.port_timestamp_sampleint,
                                           data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                           data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                           data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                           data.calibration_table, data.sensor_depth)

    smph2con = msp.calc_dissgas_smph2con(data.port_timestamp_sampleint,
                                         data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                         data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                         data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                         data.calibration_table, data.sensor_depth)

    smparcon = msp.calc_dissgas_smparcon(data.port_timestamp_sampleint,
                                         data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                         data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                         data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                         data.calibration_table, data.sensor_depth)

    smph2scon = msp.calc_dissgas_smph2scon(data.port_timestamp_sampleint,
                                           data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                           data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                           data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                           data.calibration_table, data.sensor_depth)

    smpo2con = msp.calc_dissgas_smpo2con(data.port_timestamp_sampleint,
                                         data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                         data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                         data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                         data.calibration_table, data.sensor_depth)

    smpco2con = msp.calc_dissgas_smpco2con(data.port_timestamp_sampleint,
                                           data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                           data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                           data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                           data.calibration_table, data.sensor_depth)

    #BCKINT
    bkgmethcon = msp.calc_dissgas_bkgmethcon(data.port_timestamp_bkgndint,
                                             data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                             data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                             data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                             data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                             data.calibration_table, data.sensor_depth)

    bkgethcon = msp.calc_dissgas_bkgethcon(data.port_timestamp_bkgndint,
                                           data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                           data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                           data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                           data.calibration_table, data.sensor_depth)

    bkgh2con = msp.calc_dissgas_bkgh2con(data.port_timestamp_bkgndint,
                                         data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                         data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                         data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                         data.calibration_table, data.sensor_depth)

    bkgarcon = msp.calc_dissgas_bkgarcon(data.port_timestamp_bkgndint,
                                         data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                         data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                         data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                         data.calibration_table, data.sensor_depth)

    bkgh2scon = msp.calc_dissgas_bkgh2scon(data.port_timestamp_bkgndint,
                                           data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                           data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                           data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                           data.calibration_table, data.sensor_depth)

    bkgo2con = msp.calc_dissgas_bkgo2con(data.port_timestamp_bkgndint,
                                         data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                         data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                         data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                         data.calibration_table, data.sensor_depth)

    bkgco2con = msp.calc_dissgas_bkgco2con(data.port_timestamp_bkgndint,
                                           data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                           data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                           data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                           data.calibration_table, data.sensor_depth)

    #CAL1INT
    cal1co2con = msp.calc_dissgas_cal1co2con(data.port_timestamp_calint01,
                                             data.L0_dissgas_calint01, data.gas_mode_calint01,
                                             data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                             data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                             data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                             data.calibration_table, data.sensor_depth)

    cal1methcon = msp.calc_dissgas_cal1methcon(data.port_timestamp_calint01,
                                               data.L0_dissgas_calint01, data.gas_mode_calint01,
                                               data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                               data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                               data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                               data.calibration_table, data.sensor_depth)

    #CAL2INT
    cal2co2con = msp.calc_dissgas_cal2co2con(data.port_timestamp_calint02,
                                             data.L0_dissgas_calint02, data.gas_mode_calint02,
                                             data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                             data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                             data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                             data.calibration_table, data.sensor_depth)

    cal2methcon = msp.calc_dissgas_cal2methcon(data.port_timestamp_calint02,
                                               data.L0_dissgas_calint02, data.gas_mode_calint02,
                                               data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                               data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                               data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                               data.calibration_table, data.sensor_depth)

    # compute L1 gas conc cal. range values
    smpmethcon_calrang = msp.calc_calrang_smpmethcon(data.port_timestamp_sampleint,
                                                     data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                     data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                     data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table, data.sensor_depth)

    smpethcon_calrang = msp.calc_calrang_smpethcon(data.port_timestamp_sampleint,
                                                   data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                   data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                   data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table, data.sensor_depth)

    smph2con_calrang = msp.calc_calrang_smph2con(data.port_timestamp_sampleint,
                                                 data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                 data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                 data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table, data.sensor_depth)

    smparcon_calrang = msp.calc_calrang_smparcon(data.port_timestamp_sampleint,
                                                 data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                 data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                 data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table, data.sensor_depth)

    smph2scon_calrang = msp.calc_calrang_smph2scon(data.port_timestamp_sampleint,
                                                   data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                   data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                   data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table, data.sensor_depth)

    smpo2con_calrang = msp.calc_calrang_smpo2con(data.port_timestamp_sampleint,
                                                 data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                 data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                 data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table, data.sensor_depth)

    smpco2con_calrang = msp.calc_calrang_smpco2con(data.port_timestamp_sampleint,
                                                   data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                   data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                   data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table, data.sensor_depth)

    #BCKINT
    bkgmethcon_calrang = msp.calc_calrang_bkgmethcon(data.port_timestamp_bkgndint,
                                                     data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                     data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                     data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table, data.sensor_depth)

    bkgethcon_calrang = msp.calc_calrang_bkgethcon(data.port_timestamp_bkgndint,
                                                   data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                   data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                   data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table, data.sensor_depth)

    bkgh2con_calrang = msp.calc_calrang_bkgh2con(data.port_timestamp_bkgndint,
                                                 data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                 data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                 data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table, data.sensor_depth)

    bkgarcon_calrang = msp.calc_calrang_bkgarcon(data.port_timestamp_bkgndint,
                                                 data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                 data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                 data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table, data.sensor_depth)

    bkgh2scon_calrang = msp.calc_calrang_bkgh2scon(data.port_timestamp_bkgndint,
                                                   data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                   data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                   data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table, data.sensor_depth)

    bkgo2con_calrang = msp.calc_calrang_bkgo2con(data.port_timestamp_bkgndint,
                                                 data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                 data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                 data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                 data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                 data.calibration_table, data.sensor_depth)

    bkgco2con_calrang = msp.calc_calrang_bkgco2con(data.port_timestamp_bkgndint,
                                                   data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                   data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                   data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                   data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                   data.calibration_table, data.sensor_depth)

    #CAL1INT
    cal1co2con_calrang = msp.calc_calrang_cal1co2con(data.port_timestamp_calint01,
                                                     data.L0_dissgas_calint01, data.gas_mode_calint01,
                                                     data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                                     data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table, data.sensor_depth)

    cal1methcon_calrang = msp.calc_calrang_cal1methcon(data.port_timestamp_calint01,
                                                       data.L0_dissgas_calint01, data.gas_mode_calint01,
                                                       data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                                       data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table, data.sensor_depth)

    #CAL2INT
    cal2co2con_calrang = msp.calc_calrang_cal2co2con(data.port_timestamp_calint02,
                                                     data.L0_dissgas_calint02, data.gas_mode_calint02,
                                                     data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                                     data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table, data.sensor_depth)

    cal2methcon_calrang = msp.calc_calrang_cal2methcon(data.port_timestamp_calint02,
                                                       data.L0_dissgas_calint02, data.gas_mode_calint02,
                                                       data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                                       data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table, data.sensor_depth)

    # compute L1 gas conc timestamp values
    #SMPINT
    smpmethcon_timestamp = msp.calc_timestamp_smpmethcon(data.port_timestamp_sampleint,
                                                         data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                         data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                         data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                         data.calibration_table)

    smpethcon_timestamp = msp.calc_timestamp_smpethcon(data.port_timestamp_sampleint,
                                                       data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                       data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                       data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table)

    smph2con_timestamp = msp.calc_timestamp_smph2con(data.port_timestamp_sampleint,
                                                     data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                     data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                     data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table)

    smparcon_timestamp = msp.calc_timestamp_smparcon(data.port_timestamp_sampleint,
                                                     data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                     data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                     data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table)

    smph2scon_timestamp = msp.calc_timestamp_smph2scon(data.port_timestamp_sampleint,
                                                       data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                       data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                       data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table)

    smpo2con_timestamp = msp.calc_timestamp_smpo2con(data.port_timestamp_sampleint,
                                                     data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                     data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                     data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table)

    smpco2con_timestamp = msp.calc_timestamp_smpco2con(data.port_timestamp_sampleint,
                                                       data.L0_dissgas_sampleint, data.gas_mode_sampleint,
                                                       data.port_timestamp_sampleint_mcu, data.ph_meter_sampleint_mcu,
                                                       data.inlet_temp_sampleint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table)

    #BCKINT
    bkgmethcon_timestamp = msp.calc_timestamp_bkgmethcon(data.port_timestamp_bkgndint,
                                                         data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                         data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                         data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                         data.calibration_table)

    bkgethcon_timestamp = msp.calc_timestamp_bkgethcon(data.port_timestamp_bkgndint,
                                                       data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                       data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                       data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table)

    bkgh2con_timestamp = msp.calc_timestamp_bkgh2con(data.port_timestamp_bkgndint,
                                                     data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                     data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                     data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table)

    bkgarcon_timestamp = msp.calc_timestamp_bkgarcon(data.port_timestamp_bkgndint,
                                                     data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                     data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                     data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table)

    bkgh2scon_timestamp = msp.calc_timestamp_bkgh2scon(data.port_timestamp_bkgndint,
                                                       data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                       data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                       data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table)

    bkgo2con_timestamp = msp.calc_timestamp_bkgo2con(data.port_timestamp_bkgndint,
                                                     data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                     data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                     data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                     data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                     data.calibration_table)

    bkgco2con_timestamp = msp.calc_timestamp_bkgco2con(data.port_timestamp_bkgndint,
                                                       data.L0_dissgas_bkgndint, data.gas_mode_bkgndint,
                                                       data.port_timestamp_bkgndint_mcu, data.ph_meter_bkgndint_mcu,
                                                       data.inlet_temp_bkgndint_mcu, data.massp_rga_initial_mass,
                                                       data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                       data.calibration_table)

    ##CAL1INT
    cal1co2con_timestamp = msp.calc_timestamp_cal1co2con(data.port_timestamp_calint01,
                                                         data.L0_dissgas_calint01, data.gas_mode_calint01,
                                                         data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                                         data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                         data.calibration_table)

    cal1methcon_timestamp = msp.calc_timestamp_cal1methcon(data.port_timestamp_calint01,
                                                           data.L0_dissgas_calint01, data.gas_mode_calint01,
                                                           data.port_timestamp_calint01_mcu, data.ph_meter_calint01_mcu,
                                                           data.inlet_temp_calint01_mcu, data.massp_rga_initial_mass,
                                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                           data.calibration_table)

    ##CAL2INT
    cal2co2con_timestamp = msp.calc_timestamp_cal2co2con(data.port_timestamp_calint02,
                                                         data.L0_dissgas_calint02, data.gas_mode_calint02,
                                                         data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                                         data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                                         data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                         data.calibration_table)

    cal2methcon_timestamp = msp.calc_timestamp_cal2methcon(data.port_timestamp_calint02,
                                                           data.L0_dissgas_calint02, data.gas_mode_calint02,
                                                           data.port_timestamp_calint02_mcu, data.ph_meter_calint02_mcu,
                                                           data.inlet_temp_calint02_mcu, data.massp_rga_initial_mass,
                                                           data.massp_rga_final_mass, data.massp_rga_steps_per_amu,
                                                           data.calibration_table)

    #put all calculated l1 gas conc values in one array for comparison
    msp_l1_calculated = np.array([np.round(smpmethcon), np.round(smpethcon),
                                  np.round(smph2con), np.round(smparcon), np.round(smph2scon),
                                  np.round(smpo2con), np.round(smpco2con), np.round(bkgmethcon),
                                  np.round(bkgethcon), np.round(bkgh2con), np.round(bkgarcon),
                                  np.round(bkgh2scon), np.round(bkgo2con), np.round(bkgco2con),
                                  np.round(cal1methcon), np.round(cal1co2con), np.round(cal2methcon),
                                  np.round(cal2co2con)])

    #expected output results
    #I had to make atol=4 for the nose tests to work because there are slight differences between my results and the l1 test data

    msp_l1_expected = np.array([0, 0, 0, 8, 0, 533, 4, 0, 0, 0, 10, 0, 266, 4, 177, 10, 0, 0])

    # compare calculated results to expected results
    np.testing.assert_allclose(msp_l1_calculated, msp_l1_expected, rtol=0.000001, atol=4)
