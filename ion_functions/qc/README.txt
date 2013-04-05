The OOI QA/QC Functions are implemented in this module. These functions
represent the automated QA/QC functions applied to various OOI data products. A
brief description of the functions and their current status are provided below.

Built-In Numpy Functions. These are not implemented as functions, although test
functions, based on the respective DPS documents, are).
     
     * modulus -- will wrap values back into the nominal range (y). Used, in
       particular, for compass adjustments (Test function completed, April
       2013). 
     
     * interp -- computes a 1-D linear interpolation of the data set (X,Y) to
       find the values of YI at points XI (Test function completed, April 2013).
     
     * polyval -- evaluates a polynomial of order N at point x, where polynomial
       has coefficients P and length(P) = N+1 (Test function competed April
       2013).

QC Flag Functions. These are the automated QA/QC functions. 
     
     * dataqc_globalrangetest -- generates a QC flag for a data point indicating
       whether it falls within a give range (Function and test function
       implemented in February 2013).
     
     * dataqc_localrangetest
     
     * dataqc_spiketest -- generates a QC flag for individual data values that
       deviate significantly from surrounding data values (Function and test
       function implemented in February 2013).
     
     * dataqc_polytrendtest -- tests a time series for whether the data contain
       a significant portion of a polynomial (Function and test function
       implemented in February 2013). 
     
     * dataqc_stuckvaluetest -- generates a QC flag for a repeated occurrence of
       one value in a time series (Function and test function implemented in
       February 2013).
     
     * dataqc_gradienttest
     
     * dataqc_combinedflags

Data Processing Functions
     
     * dataqc_solorelev
     
     * dataqc_condcompress
     
Additional Functions. Available in ../utils.py providing Matlab-based test
utilities (e.g. isvector) used in the various QC functions. 
     