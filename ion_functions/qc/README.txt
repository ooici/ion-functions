The OOI QA/QC Functions are implemented in this module. These functions
represent the automated QA/QC functions applied to various OOI data products. A
brief description of the functions and their current status are provided below.

Built-In Numpy Functions. These are not implemented as functions, although test
functions, based on the respective DPS documents, are).
     
     * modulus -- will wrap values back into the nominal range (y). Used, in
       particular, for compass adjustments (test function completed, April
       2013). 
     
     * interp -- computes a 1-D linear interpolation of the data set (x,y) to
       find the values of yi at points xi (test function completed, April
       2013).
     
     * polyval -- evaluates a polynomial of order n at point x, where
       polynomial has coefficients p and length(p) = n+1 (test function
       completed April 2013).

QC Flag Functions. These are the generic, automated QA/QC functions that will
be applied to various OOI data products. 
     
     * dataqc_globalrangetest -- generates a QC flag for a data point
       indicating whether it falls within a give range (function and test
       function implemented in February 2013).
     
     * dataqc_localrangetest -- quality control algorithm testing if
       measurements fall into a user-defined valid range. This range is not
       constant but varies with measurement location (function and test
       function implemented in April 2013).
     
     * dataqc_spiketest -- generates a QC flag for individual data values that
       deviate significantly from surrounding data values (function and test
       function implemented in February 2013).
     
     * dataqc_polytrendtest -- tests a time series for whether the data contain
       a significant portion of a polynomial (function and test function
       implemented in February 2013). 
     
     * dataqc_stuckvaluetest -- generates a QC flag for a repeated occurrence
       of one value in a time series (function and test function implemented in
       February 2013).
     
     * dataqc_gradienttest -- quality control algorithm testing if changes
       between successive data points fall within a certain range (function and
       test function implemented in April 2013).
     
     * dataqc_propogateflags -- propagates "bad" qc flags (from an arbitrary
       number of source datasets) to another (derived) dataset (function and
       test function implemented in April 2013).

Data Processing Functions.
     
     * dataqc_solarelevation -- Computes instantaneous no-sky solar radiation
       and altitude from date and time stamp and position data (function and
       test function implemented in April 2013).
     
     * dataqc_condcompress -- Implementation of the Sea-Bird conductivity
       compressibility correction, scaling the input conductivity based on
       ratio of the original pressure and the updated pressure. 
     
Additional Functions, available in ../utils.py, provide Matlab-based test
utilities (e.g. isvector) used in the various QC functions. These are intended
to duplicate the functionality of these functions as called by the original DPS
code. 
     