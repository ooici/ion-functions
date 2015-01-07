The OOI Data Functions are implemented in this group of modules. These modules,
and the functions therein, represent the transforms and calculations applied
to parameters uploaded to the system via either dataset or instrument agents,
and are used to determine various OOI data products. A brief description of the
modules and their contents are provided below. Each set of modules is grouped
according to its Instrument Family as defined in SAF. Note that the functions
calculating SFLPRES data products have been misclassified, and will be found
in the Seafloor Properties (SFL) module instead of the Seafloor Pressure (PRS)
module.
        
Alphabetical by instrument family Public ID.

CO2: Partial Pressure CO2

    * co2_functions.py -- Covers calculation of the L1 PCO2WAT data product
      and extraction of the associated L0 data products.
      
      [TODO: Fill out with module names and descriptions]
      pco2_abs434_ratio -- extracts the CO2ABS1_L0 measurement from the data array
      pco2_abs620_ratio -- extracts the CO2ABS2_L0 measurement from the data array
      pco2_blank -- scales raw blank measurements
      pco2_thermistor -- converts raw thermistor measurements to degrees C
      pco2_pco2wat -- calculates the PCO2WAT_L2 data product. Serves as a
            wrapper function and calls pco2_calc_pco2.
      pco2_calc_pco2 -- called by pco2_pco2wat.
      pco2_ppressure -- computes PCO2ATM_L1 or PCO2SSW_L1 given inputs of
            either the XCO2ATM_L0 or XCO2SSW_L0 and the Gas Stream Pressure
            (PRESAIR_L0).
      pco2_pco2wat -- computes CO2FLUX_L2 using PCO2ATM_L1 and PCO2SSW_L1 as
            inputs.

CTD: Conductivity, Temperature, Depth

    * ctd_functions.py -- Covers calculation of the L1 CONDWAT, TEMPWAT and
      PRESWAT data products, and the L2 PRACSAL and DENSITY data products.
      This module includes the following functions for the CTD instrument
      classes (series are given in parentheses):
      
      ctd_sbe16plus_condwat -- calculates CONDWAT_L1 from CTDBP(CDEFNO) and
            CTDPF(AB)
      ctd_sbe16plus_tempwat -- calculates TEMPWAT_L1 from CTDBP(CDEFNO) and
            CTDPF(AB)
      ctd_sbe16plus_preswat -- calculates PRESWAT_L1 from CTDBP(CDEF) and
            CTDPF(AB)
      ctd_sbe16digi_preswat -- calculates PRESWAT_L1 from CTDBP(N and O
            series) only
      ctd_sbe37im_condwat   -- calculates CONDWAT_L1 from CTDMO all (GHQR)
      ctd_sbe37im_tempwat   -- calculates TEMPWAT_L1 from CTDMO all (GHQR)
      ctd_sbe37im_preswat   -- calculates PRESWAT_L1 from CTDMO all (GHQR)
      ctd_sbe52mp_condwat   -- calculates CONDWAT_L1 from CTDPF (CKL)
      ctd_sbe52mp_tempwat   -- calculates TEMPWAT_L1 from CTDPF (CKL)
      ctd_sbe52mp_preswat   -- calculates PRESWAT_L1 from CTDPF (CKL)
      ctd_pracsal -- calculates PRACSAL_L2 from all CTDs
      ctd_density -- calculates DENSITY_L2 from all CTDs
      
      CTDAV (all series) and CTDGV (all series) produce the L1 data
      products directly.

      New series for which DPAs are not written (probably not required):
      CTDBP (P) SBE16Plus-IM V2
      CTDPF (J) FastCAT 49

DO2: Dissolved Oxygen

    * do2_functions.py -- calculates L1 & L2 DOCONCS and L2 DOCONCF data
      products.  This module includes:
    
      do2_SVU  -- calculates DOCONCS_L1 from DOSTAs
      do2_salinity_correction -- calculates DOCONCS_L2 from DOSTAs
      do2_dofst_volt -- calculates DOCONCF_L2 from a DOFST-As (SBE 43)
      do2_dofst_frequency -- calculates DOCONCF_L2 from a DOFST-Ks (SBE 43F)

FDC: Direct Covariance Flux

    * fdc_functions.py -- Covers calculation of the FDCHP data products. More
      extensive documentation will be found at module top of fdc_functions.py.

      This module includes the following functions:

      Functions to compute the L1 FDCHP data products:
        fdc_tmpatur:        TMPATUR
        fdc_windtur_north:  WINDTUR-VLN
        fdc_windtur_up:     WINDTUR-VLU
        fdc_windtur_west:   WINDTUR-VLW

      Functions to compute the L2 FDCHP data products:
        fdc_fluxhot:            FLUXHOT
        fdc_fluxmom_alongwind:  FLUXMOM-U
        fdc_fluxmom_crosswind:  FLUXMOM-V

      Functions to compute the auxiliary time base data products:
        fdc_time_L1:  TIME_L1-AUX
        fdc_time_L2:  TIME_L2-AUX

FLO: Fluorometer

    * flo_functions.py -- Covers calculation of the L1 CDOMFLO, CHLAFLO and
      FLUBSCT data products from the FLORT and FLORD instrument classes.
      This module contains the following functions:
      
      ### Core functions
      flo_bback_total -- calculates FLUBSCT-BBACK_L1 from FLORT and FLORD
            (all series)
      flo_scat_seawater -- calculates the scattering coefficient of
            seawater from Zhang et al 2009, a required metadata parameter
            for FLUBSCT
      flo_beta -- calculates FLUBSCT-BETA_L1 from FLORT and FLORD
            (all series)
      flo_cdom -- calculates CDOMFLO_L1 from FLORT (all series)
      flo_chla -- calculates CHLAFLO_L1 from FLORT and FLORD (all series)

      ### Auxiallary functions
      flo_scale_and_offset -- applies scale and offset calculations used by
            all WET Labs ECO instruments
      flo_zhang_scatter_coeff -- calculates the volume scattering
            and total scattering coefficients of seawater based on
            calculations defined in Zhang et al 2009. Used in
            flo_bback_total.
      flo_refractive_index -- used by flo_zhang_scatter_coeff
      flo_isotherm_compress -- used by flo_zhang_scatter_coeff
      flo_density_seawater -- used by flo_zhang_scatter_coeff

HYD: Hydrophone

    * hyd_functions.py -- calculates the L1 HYDAPLF and HYDAPBB data products.
      This module includes:
      
      hyd_bb_acoustic_pwaves -- calculates HYDAPBB_L1 from the HYDBB
      hyd_lf_acoustic_pwaves -- calculates HYDAPLF_L1 from the HYDLF

MET: Meteorology

    * met_functions.py -- Covers calculation of the METBK data products. More
      extensive documentation will be found at module top of met_functions.py.
      
      This module includes functions to calculate the following data products,
      listed in alphabetical order. Except as noted, the functions are named
      as "met_prdname", so that the function calculating BARPRES_L1 is named
      met_barpres.

        BARPRES_L1
        BUOYFLS_L2:  added DPA to match FDCHP, not in original DPS
        BUOYFLX_L2:  added DPA to match FDCHP, not in original DPS
        CURRENT_DIR (meta)
        CURRENT_SPD (meta)
        FRSHFLX_L2
        HEATFLX_L2
        LATNFLX_L2
        MOMMFLX_L2
        NETLIRR_L2
        NETSIRR_L2 (this may operationally be an L1 product)
        RAINFLX_L2
        RAINRTE_L2
        RELWIND_DIR-AUX (meta)
        RELWIND_SPD-AUX (meta)
        SALSURF_L2
        SENSFLX_L2
        SPECHUM_L2
        SPHUM2M_L2
        STABLTY_L2:  metadata
        TEMPA2M_L2
        TEMPSKN_L2:  metadata
        TIMEFLX-AUX (meta)
        WIND10M_L2
        WINDAVG-VLE_L1  (function name: met_windavg_mag_corr_east)
        WINDAVG-VLN_L1  (function name: met_windavg_mag_corr_north)

MSP: Mass Spectrometer

    * msp_functions.py -- calculates all of the L1 and L2 data products associated
      with the MASSP instrument. In addition, it includes all of the functions
      needed to calculate the necessary auxiliary data products including timestamp
      and calibration range information.

NIT: Nitrate

    * nit_functions.py -- calculates the L2 temperature and salinity corrected
      dissolved nitrate concentration (NITRTSC) data product.

OBS: Ocean Bottom Seismometer

    * obs_functions.py -- calculates the L1 GRNDVEL, GRNDACC and SGRDVEL
      data products.  This module includes:
    
      obs_bb_ground_velocity -- calculates GRNDVEL_L1 from OBSBB
      obs_bb_ground_acceleration -- calculates GRNDACC_L1 from OBSBB
      obs_sp_ground_velocity -- calculates SGRDVEL_L1 from OBSBB

OPT: Optical Properties

    * opt_functions.py -- Covers calculation of the data products from the
            OPTAA, PARAD, and SPKIR instrument classes. 
      
      This module includes the following functions:
      
      OPTAA core functions used to calculate primary data products OPTATTN
      and OPTABSN:
            opt_beam_attenuation -- wrapper function to calculate
              OPTATTN_L2 from functions below.
            opt_optical_absorption -- wrapper function to calculate
              OPTABSN_L2 from functions below.
            opt_internal_temp -- Calculates internal instrument temperature.
            opt_pd_calc -- Converts raw measurements to either beam
              attenuation or optical absorbtion depending on input.
            opt_tempsal_corr -- Applies temperature and salinity corrections.
            opt_scatter_corr -- Applies proportional scatter correction to
                OPTABSN.
      OPTAA auxiliary functions coded but not used to calculate primary data
      products:
            opt_pressure -- Calculates in situ pressure, if auxiliary sensor is
                installed. this product is not used in the calculation of
                OPTATTN nor OPTABSN.
            opt_external_temp -- Calculates external in situ temperature, if
                auxiliary sensor is installed. Normally this product is not
                used in the calculation of OPTATTN nor OPTABSN; rather,
                TEMPWAT_L1 from a co-located CTD would be used.

    * opt_functions_tscor.py -- Provides an array of the published
      wavelength-dependent temperature and salinity correction factors
      for the OPTAA data products. Values are provided by the vendor.
      Required to calculate OPTATTN_L2 and OPTABSN_L2.
      
      OPT functions calculating the PAR OPTPARW_L1 data product:
            opt_par_satlantic -- computes OPTPARW_L1 from data acquired from a
                Satlantic sensor. 
            opt_par_biospherical_mobile -- computes OPTPARW_L1 from data
                acquired from a Biospherical QSP-2100.
            opt_par_biospherical_wfp -- computes OPTPARW_L1 from data acquired
                from a Biospherical QSP-2200
            opt_par_wetlabs -- computes OPTPARW_L1 from data acquired from a
                WET Labs ECO PAR sensor on CSPPs
       
      OPT function calculating the SPKIR SPECTIR_L1 data product:
            opt_ocr507_irradiance -- computes the downwelling spectral
                irradiance SPECTIR_L1 data product.

PHS: pH

    * ph_functions.py -- Covers calculation of the L1 PHWATER data product
      and extraction of the associated L0 data products.
      
      ph_434_intensity -- extracts the 23 measurements of the 434 nm signal
            intensity [PH434SI_L0] from the 92 light measurements.
      ph_578_intensity -- extracts the 23 measurements of the 578 nm signal
            intensity [PH578SI_L0] from the 92 light measurements.
      ph_battery -- converts the raw battery measurements from counts to
            volts.
      ph_thermistor -- converts the raw thermistor measurements from counts
            [ABSTHRM_L0] to degrees Centigrade.
      ph_calc_phwater -- calculates the OOI Level 1 pH of seawater core
            data product [PHWATER_L1].

PRS: Seafloor Pressure

      NOTE: the SFLPRES data products from the PRESF class of instruments are
      classified in SAF as belonging to the PRS family. However, the functions
      calculating the SLFPRES data products are contained in the SFL module.

    * prs_functions -- Covers calculation of the L1 (BOTTILT) and L2 (BOTSFLU)
      data products collected from the BOTPT instruments. No DPA function is
      required for the BOTPRES_L1 data product ("L0=L1").
      
      prs_bottilt_ccmp -- computes the BOTTILT-CCMP_L1 data product
      prs_bottilt_tmag -- computes the BOTTILT-TMAG_L1 data product
      prs_bottilt_tdir -- computes the BOTTILT-TDIR_L1 data product
       
    * prs_functions_ccmp.py -- dictionary of compass calibration data required by
                               prs_bottilt_ccmp.
    * prs_functions_ccmp_lily_compass_cals.py -- re-organized dictionary of compass
                               calibration data required by the faster version of
                               prs_bottilt_ccmp; the code for the faster version is
                               currently commented out within the prs_bottilt_ccmp
                               code.

      BOTSFLU data products: in development.

SFL: Seafloor Properties

      NOTE: the SFLPRES data products from the PRESF class of instruments are
      contained in this module, even though they are classified in SAF as
      belonging to the PRS (Seafloor Pressure) family.

    * sfl_functions.py -- Covers calculation of L1 and L2 data products collected
      from the THSPH, TRHPH, and PRESF instruments (even though PRESF is in the
      PRS family of instruments).

      THSPH functions and data products:
        sfl_thsph_ph:             THSPHPH-PH_L2
        sfl_thsph_ph_acl:         THSPHPH-PH-ACL_L2
        sfl_thsph_ph_noref:       THSPHPH-PH-NOREF_L2
        sfl_thsph_ph_noref_acl:   THSPHPH-PH-NOREF-ACL_L2
        sfl_thsph_sulfide:        THSPHHS_L2
        sfl_thsph_hydrogen:       THSPHHC_L2

        sfl_thsph_temp_th:        THSPHTE-TH_L1
        sfl_thsph_temp_tl:        THSPHTE-TL_L1
        sfl_thsph_temp_tch:       THSPHTE-TCH_L1
        sfl_thsph_temp_tcl:       THSPHTE-TCL_L1
        sfl_thsph_temp_ref:       THSPHTE-REF_L1
        sfl_thsph_temp_int        THSPHTE-INT_L1

      TRHPH functions and data products:
        sfl_trhph_vfltemp:                 TRHPHTE_L1
        sfl_trhph_vfl_thermistor_temp:     TRHPHTE-T_TS-AUX
        sfl_trhph_vflorp:                  TRHPHEH_L1
        sfl_trhph_chloride:                TRHPHCC_L2

      PRESF functions and data products:
        sfl_sflpres_rtime:                 SFLPRES-RTIME_L1
        sfl_sflpres_tide:                  SFLPRES-TIDE_L1
        sfl_sflpres_wave:                  SFLPRES-WAVE_L1
  
    * sfl_functions_surface.py -- Recreates the 3 arrays of temperature
      (tdat), conductivty (sdat) and chloride (cdat) from the Larson et al
      2007 derived calibration surface provided to CI in a Matlab file
      (Larson_2007surface.mat). This file is required by sfl_trhph_chloride.

VEL: Water Velocity

    * adcp_functions.py -- Covers calculation of the L1 VELPROF, VELTURB
    and ECHOINT data products from the tRDI ADCPs used throughout the
    program (ADCPA, ADCPS, ADCPT and VADCP). This module includes:
    
      **** For instruments programmed in beam coordinates
      adcp_beam_eastward -- calculates VELPROF-VLE_L1
      adcp_beam_northward -- calculates VELPROF-VLN_L1
      adcp_beam_vertical -- calculates VELPROF-VLU_L1
      adcp_beam_error -- calculates VELPROF-ERR_L1
        
      **** For instruments programmed in earth coordinates
      adcp_earth_eastward -- calculates VELPROF-VLE_L1
      adcp_earth_northward -- calculates VELPROF-VLN_L1
      adcp_earth_vertical -- calculates VELPROF-VLU_L1
      adcp_earth_error -- calculates VELPROF-ERR_L1
        
      **** For the VADCP programmed in beam coordinates
      vadcp_beam_eastward -- calculates VELTURB-VLE_L1
      vadcp_beam_northward -- calculates VELTURB-VLN_L1
      vadcp_beam_vertical -- calculates VELTURB-VLU_L1
      vadcp_beam_error -- calculates VELTURB-ERR_L1
        
      **** For all ADCPS
      adcp_backscatter -- calculates ECHOINT-B1_L1,
                          calculates ECHOINT-B2_L1,
                          calculates ECHOINT-B3_L1,
                          calculates ECHOINT-B4_L1 for all ADCPs
              
      **** Base functions used by above functions
      adcp_beam2ins -- applies the beam to instrument transform using a 4
            beam solution for instruments programmed in beam coordinates
      vadcp_beam2ins -- applies the beam to instrument transform using a 5
            beam solution for the VADCP programmed in beam coordinates
      adcp_ins2earth -- applies the instrument to Earth transform for all
            instruments originally programmed in beam coordinates.
    
    * vel_functions.py -- Covers calculation of the L1 VEL3D Eastward
    and Northward component data products from the VEL3D-B, VEL3D-CD,
    VEL3D-K, and VELPT instruments.  This module includes
    
      nobska_mag_corr_east -- calculates VELPTTU-VLE_L1 from VEL3D-B
      nobska_mag_corr_north -- calculates VELPTTU-VLN_L1 from VEL3D-B
      nobska_scale_up_vel -- calculates VELPTTU-VLU_L0 from VEL3D-B

      nortek_mag_corr_east -- calculates VELPTTU-VLE_L1 from VEL3D-CDK &
          VELPT
      nortek_mag_corr_north -- calculates VELPTTU-VLN_L1 from VEL3D-CDK &
          VELPT
      nortek_up_vel -- passes through VELPTTU-VLU_L1 from VEL3D-CDK &
          VELPT
      
WAV: Surface Wave Spectra

    * wav_functions.py -- Covers calculation of the WAVSS data products WAVSTAT
      from the TRIAXYS instrument.
      
      This module includes the following functions:

      wav_triaxys_dir_freq:                            WAVSTAT-FDS_L1
      wav_triaxys_nondir_freq:                         WAVSTAT-FND_L1
      wav_triaxys_buoymotion_time:                     WAVSTAT-MOTT_L1
      wav_triaxys_correct_mean_wave_direction:         WAVSTAT-D_L2
      wav_triaxys_correct_directional_wave_direction:  WAVSTAT-DDS_L2
      wav_triaxys_magcor_buoymotion_x:                 WAVSTAT-MOTX_L1
      wav_triaxys_magcor_buoymotion_y:                 WAVSTAT-MOTY_L1




Additional Functions, available in generic_functions.py, provide for transforms
and calculations that apply to multiple instrument families.
