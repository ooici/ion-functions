The .m files in this directory are the Matlab code found in the Data
Product Specifications (DPS) for the Dissolved Oxygen (DO2) family of
OOI instruments.

dosv - calculates Salinity and Pressure corrected Oxygen Concentration
    through 2 steps.
    1. Stern-Volmer-Uchida equation. Uses phase, temperature, and SVU
    calibration coefficients to calculate Temperature corrected O2
    concentration.
    2. Uses data from a collocated CTD to correct for
    Salinity and Pressure