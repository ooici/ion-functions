The OOI Data Functions are implemented in this group of modules. These modules,
and the functions therein, represent the transforms and calculations applied
to parameters uploaded to the system via either dataset or instrument agents,
and are used to determine various OOI data products. A brief description of the
modules and their current status are provided below. Each set of modules is
grouped according to its Instrument Family as defined in SAF. 
        
Conductivity, Temperature, Depth (CTD)

     * ctd_functions.py -- covers calculation of the L1 CONDWAT, TEMPWAT and
       PRESWAT data products, and the L2 PRACSAL and DENSITY data products.
       
Dissolved Oxygen

     * do2_functions.py
     
Optical Properties (OPT)

     * opt_functions.py

Partial Pressure CO2 (CO2)

     * co2_functions.py
     
pH (pH)

     * ph_functions.py

Seafloor Properties (SFL)

     * sfl_functions.py
          
Water Velocity (VEL)

     * adcp_functions.py 
     
     * vel_functions.py

Additional Functions, available in generic_functions.py, provide for transforms
and calculations that apply to multiple instrument families.
