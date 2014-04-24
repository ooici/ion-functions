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

Dissolved Oxygen (DO)

     * do2_functions.py -- calculates L1 & L2 DOCONCS and L2 DOCONCF data
       products.  This module includes:
     
       do2_SVU  -- calculates DOCONCS_L1 from DOSTAs
       do2_salinity_correction -- calculates DOCONCS_L2 from DOSTAs
       do2_dofst_volt -- calculates DOCONCF_L2 from a DOFST-As (SBE 43)
       do2_dofst_frequency -- calculates DOCONCF_L2 from a DOFST-Ks (SBE 43F)
     
Optical Properties (OPT)

     * opt_functions.py -- Covers calculation of the L2 OPTATTN and OPTABSN
       data products, as well as products from various irradiance sensors.
       
       This module includes the following functions:
       
        OPTAA core functions used to calculate primary data products OPTATTN
        and OPTABSN:
            opt_beam_attenuation -- wrapper function to calculate OPTATTN_L2
                from functions below.
            opt_optical_absorption -- wrapper function to calculate OPTABSN_L2
                from functions below.
            opt_internal_temp -- Calculates internal instrument temperature.
            opt_pd_calc -- Converts raw measurements to either beam attenuation
                or optical absorbtion depending on input.
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
       
        opt_ocr507_irradiance -- computes the downwelling spectral irradiance
                SPECTIR_L1 data product.

Partial Pressure CO2 (CO2)

        * co2_functions.py -- Covers calculation of the L1 PCO2WAT data product
          and extraction of the associated L0 data products.
          
          [TODO: Fill out with module names and descriptions]
          pco2_abs434_ratio
          pco2_abs620_ratio
          pco2_abs434_blank
          pco2_abs620_blank
          pco2_thermistor
          pco2_pco2wat
          pco2_calc_pco2
     
pH (pH)

        * ph_functions.py -- Covers calculation of the L1 PHWATER data product
          and extraction of the associated L0 data products.
          
          [TODO: Fill out with module names and descriptions]
          ph_434_intensity
          ph_578_intensity
          ph_thermistor
          ph_phwater

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

        * adcp_functions.py -- Covers calculation of the L1 VELPROF and ECHOINT
          data products from ADCPs. This module includes:
          
          adcp_beam_eastward
          adcp_beam_northward
          adcp_beam_vertical
          adcp_beam_error
          adcp_earth_eastward
          adcp_earth_northward
          adcp_beam2ins
          adcp_ins2earth
        
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
       
Meteorology (MET)

        * met_functions.py -- Covers calculation of the L1 WINDAVG Eastward
          and Northward component data products from METBK instruments.
          This module includes:
          
          windavg_mag_corr_east -- calculates WINDAVG-VLE_L1 from METBKs
          windavg_mag_corr_north -- calculates WINDAVG-VLN_L1 from METBKs

Additional Functions, available in generic_functions.py, provide for transforms
and calculations that apply to multiple instrument families.
