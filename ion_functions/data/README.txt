The OOI Data Functions are implemented in this group of modules. These modules,
and the functions therein, represent the transforms and calculations applied
to parameters uploaded to the system via either dataset or instrument agents,
and are used to determine various OOI data products. A brief description of the
modules and their contents are provided below. Each set of modules is grouped
according to its Instrument Family as defined in SAF.
        
Conductivity, Temperature, Depth (CTD)

    * ctd_functions.py -- Covers calculation of the L1 CONDWAT, TEMPWAT and
      PRESWAT data products, and the L2 PRACSAL and DENSITY data products.
      This module includes the following functions:
      
      ctd_sbe16plus_condwat -- calculates CONDWAT_L1 from CTDBP(CDEFNO) and
            CTDPF(AB)
      ctd_sbe16plus_tempwat -- calculates TEMPWAT_L1 from CTDBP(CDEFNO) and
            CTDPF(AB)
      ctd_sbe16plus_preswat -- calculates PRESWAT_L1 from CTDBP(CDEF) and
            CTDPF(AB)
      ctd_sbe16digi_preswat -- calculates PRESWAT_L1 from CTDBP(N and O
            series) only
      ctd_sbe37im_condwat   -- calculates CONDWAT_L1 from CTDMO all series
      ctd_sbe37im_tempwat   -- calculates TEMPWAT_L1 from CTDMO all series
      ctd_sbe37im_preswat   -- calculates PRESWAT_L1 from CTDMO all series
      ctd_sbe52mp_condwat   -- calculates CONDWAT_L1 from CTDPF all series
      ctd_sbe52mp_tempwat   -- calculates TEMPWAT_L1 from CTDPF all series
      ctd_sbe52mp_preswat   -- calculates PRESWAT_L1 from CTDPF all series
      ctd_pracsal -- calculates PRACSAL_L2 from all CTDs
      ctd_density -- calculates DENSITY_L2 from all CTDs
      
      CTDAV (all series) and CTDGV (all series) produce the L1 data
      products directly.

Dissolved Oxygen (DO)

    * do2_functions.py -- calculates L1 & L2 DOCONCS and L2 DOCONCF data
      products.  This module includes:
    
      do2_SVU  -- calculates DOCONCS_L1 from DOSTAs
      do2_salinity_correction -- calculates DOCONCS_L2 from DOSTAs
      do2_dofst_volt -- calculates DOCONCF_L2 from a DOFST-As (SBE 43)
      do2_dofst_frequency -- calculates DOCONCF_L2 from a DOFST-Ks (SBE 43F)
     
Fluorometer (FLO)

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

Hydrophone (HYD)

    * hyd_functions.py -- calculates the L1 HYDAPLF and HYDAPBB data products.
      This module includes:
      
      hyd_bb_acoustic_pwaves -- calculates HYDAPBB_L1 from the HYDBB
      hyd_lf_acoustic_pwaves -- calculates HYDAPLF_L1 from the HYDLF
      
Meteorology (MET)

    * met_functions.py -- Covers calculation of the L1 WINDAVG Eastward
      and Northward component data products from METBK instruments.
      This module includes:
      
      windavg_mag_corr_east -- calculates WINDAVG-VLE_L1 from METBKs
      windavg_mag_corr_north -- calculates WINDAVG-VLN_L1 from METBKs

Ocean Bottom Seismometer (OBS)

    * obs_functions.py -- calculates the L1 GRNDVEL, GRNDACC and SGRDVEL
      data products.  This module includes:
    
      obs_bb_ground_velocity -- calculates GRNDVEL_L1 from OBSBB
      obs_bb_ground_acceleration -- calculates GRNDACC_L1 from OBSBB
      obs_sp_ground_velocity -- calculates SGRDVEL_L1 from OBSBB

Optical Properties (OPT)

    * opt_functions.py -- Covers calculation of the L2 OPTATTN and OPTABSN
      data products, as well as products from various irradiance sensors.
      
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
      
      OPT functions calculating the PAR OPTPARW_L1 data product:
            opt_par_satlantic -- computes OPTPARW_L1 from data acquired from a
                Satlantic sensor. 
            opt_par_biospherical_mobile -- computes OPTPARW_L1 from data
                acquired from a Biospherical QSP-2100.
            opt_par_biospherical_wfp -- computes OPTPARW_L1 from data acquired
                from a Biospherical QSP-2200
       
            opt_ocr507_irradiance -- computes the downwelling spectral
                irradiance SPECTIR_L1 data product.

Partial Pressure CO2 (CO2)

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
     
pH (pH)

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

Seafloor Pressure (PRS)

    * prs_functions -- Covers calculation of the L1 data products collected
      from the BOTPT (BOTTILT), PRESF and PREST (both SFLPRES) instruments.
      Note, functions for PRESF and PREST are implemented as Preload
      ParameterFunctions, rather than herein. Also, no function is required
      for the BOTPRES_L1 data product, as that value is output directly by
      the instrument.
      
      prs_bottilt_ccmp -- computes the BOTTILT-CCMP_L1 data product
      prs_bottilt_tmag -- computes the BOTTILT-TMAG_L1 data product
      prs_bottilt_tdir -- computes the BOTTILT-TDIR_L1 data product
       
Seafloor Properties (SFL)

    * sfl_functions.py -- Covers calculation of L1 data products collected
      from the TRHPH instrument (TRHPHTE, TRHPHCC, TRHPHEH)
      
      sfl_trhph_vfltemp -- computes the TRHPHTE_L1 data product
      sfl_trhph_chloride -- computes the TRHPHCC_L2 data product
      sfl_trhph_vflorp -- computes the TRHPHEH_L1 data product
      sfl_trhph_vfl_thermistor_temp -- computes a diagnostic data product,
          which, however, is not a core data product.
        
    * sfl_functions_surface.py -- Recreates the 3 arrays of temperature
      (tdat), conductivty (sdat) and chloride (cdat) from the Larson et al
      2007 derived calibration surface provided to CI in a Matlab file
      (Larson_2007surface.mat). 
          
Water Velocity (VEL)

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
      
Additional Functions, available in generic_functions.py, provide for transforms
and calculations that apply to multiple instrument families.
