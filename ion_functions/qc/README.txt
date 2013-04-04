The OOI QA/QC Functions represented herein are described below along with their
status. These functions represent the automated QA/QC functions applied to OOI
data products.

Built-In Numpy Functions (not implemented as functions, although test functions,
based on DPS documents, are).
     * modulus -- Will wrap values back into the nominal range (y). Used, in
       particular, for compass adjustments. 
     * interp1 -- 1-D linear interpolation.
     * polyval -- evaluates a polynomial of order N at point x, where polynomial
       has coefficients P and length(P) = N+1.

QC Flag Functions
     * dataqc_globalrangetest -- 
     * dataqc_localrangetest
     * dataqc_spiketest
     * dataqc_polytrendtest
     * dataqc_stuckvaluetest
     * dataqc_gradienttest
     * dataqc_combinedflags

Data Processing Functions
     * dataqc_solorelev
     * dataqc_condcompress
     