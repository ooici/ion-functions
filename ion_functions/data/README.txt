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
       
       ctd_sbe16plus_condwat -- calculates CONDWAT_L1 from CTDBP series
       ctd_sbe16plus_tempwat -- calculates TEMPWAT_L1 from CTDBP series
       ctd_sbe16plus_preswat -- calculates PRESWAT_L1 from CTDBP series
       ctd_sbe16digi_preswat -- calculates PRESWAT_L1 from CTDBP series
       ctd_pracsal -- calculates PRACSAL_L2 from all CTDs
       ctd_density -- calculates DENSITY_L2 from all CTDs
       
       CTDMO, CTDGV, CTDAV produce either the L1 data products directly (CTDAV
       and CTDPF), or the functions are represented by numeric expressions in
       Preload (CTDMO). CTDPF is an unknown as of 2013-05-10, since the DPSs
       don't cover this one.
       
Dissolved Oxygen

     * do2_functions.py
     
       [TODO: Fill out with module names and descriptions]
     
Optical Properties (OPT)

     * opt_functions.py -- Covers calculation of the L2 OPTATTN and OPTABSN
       data products. This module includes the following functions:
       
       opt_beam_attenuation -- wrapper function to calculate OPTATTN_L2 from
          functions below.
       opt_optical_absorption -- wrapper function to calculate OPTABSN_L2 from
          functions below.
       opt_pressure -- Calculates measured pressure, if sensor is installed.
       opt_internal_temp -- Calculates internal instrument temperature.
       opt_external_temp -- Calculates external, in situ temperature.
       opt_pd_calc -- Converts raw measurements to either beam attenuation or
          optical absorbtion
       opt_tempsal_corr -- Applies temperature and salinity corrections.
       opt_scatter_corr -- Applies proportional scatter correction to OPTABSN.
       
     * opt_functions_tscor.py -- Provides an array of the published
       wavelength-dependent temperature and salinity correction factors. Values
       are provided by the vendor.

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
       ParameterFunctions, rather than herein.
       
       [TODO: fill out with module names and descriptions]
       prs_bottilt_ccmp
       prs_bottilt_tmag
       prs_bottilt_tdir
       
Seafloor Properties (SFL)

     * sfl_functions.py -- Covers calculation of L1 data products collected
       from the TRHPH instrument (TRHPHTE, TRHPHCC, TRHPHEH)
       
       [TODO: Fill out with module names and descriptions]
       sfl_trhph_vfltemp
       sfl_trhph_chlorconc
       sfl_trhph_chloride
     
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
       adcp_magvar
     
     * vel_functions.py
     
       [TODO: Fill out with module names and descriptions]

Additional Functions, available in generic_functions.py, provide for transforms
and calculations that apply to multiple instrument families.