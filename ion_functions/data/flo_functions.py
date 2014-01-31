#!/usr/bin/env python
"""
@package ion_functions.data.flo_functions
@file ion_functions/data/flo_functions.py
@author Craig Risien
@brief Module containing calculations related to instruments in the fluorescence family.
"""

# Import NumExpr

import numexpr as ne
    
    
def flo_scale_and_offset(counts_output, counts_dark, scale_factor):
    """
    Description:

        This scale and offset function is a simple numeric expression that can
        be applied to the CHLAFLO, CDOMFLO, FLUBSCT data products
        
    Implemented by:

        2014-01-30: Craig Risien. Initial Code
        
    Usage:

        value = flo_scale_and_offset(counts_output, counts_dark, scale_factor)

            where

        value = output value
        counts_output = calibrated sample voltage output [counts]
        counts_dark = measured signal output of fluormeter in clean water with black tape over the detector [counts]
        scale_factor = multiplier [micrograms / liter / volts]

    References:
    
        N/A
    """

    value = ne.evaluate('(counts_output - counts_dark) * scale_factor')
    return value


def flo_chla(counts_output, counts_dark, scale_factor):
    """
    Description:

        The OOI Level1 Fluorometric Chlorophyll-a Concentration core data
        product is a measure of how much light has been re-emitted after
        being absorbed by chlorophyll-a molecules found in all phytoplankton.
        By measuring the intensity and nature of this fluorescence,
        phytoplankton biomass can be estimated. The concentration of
        chlorophyll-a is a proxy for the abundance of phytoplankton in the
        water column, and thus the amount of primary productivity that can be
        empirically achieved. Chlorophyll absorbs photons in the visible
        spectrum (400-700nm) and fluoresces visible blue light.

    Implemented by:

        2014-01-30: Craig Risien. Initial Code
        
    Usage:

        chla_conc = flo_chla(counts_output, counts_dark, scale_factor)

            where

        chla_conc = Fluorometric Chlorophyll-a Concentration (CHLAFLO_L1) [micrograms / liter]
        counts_output = calibrated sample voltage output (CHLAFLO_L0) [counts]
        counts_dark = measured signal output of fluormeter in clean water with black tape over the detector [counts]
        scale_factor = multiplier [micrograms / liter / volts]

    References:
    
        OOI (2012). Data Product Specification for Fluorometric Chlorophyll-a Concentration. Document Control Number
        1341-00530. https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00530_Data_Product_SPEC_CHLAFLO_OOI.pdf)
    """

    chla_conc = flo_scale_and_offset(counts_output, counts_dark, scale_factor)
    return chla_conc


def flo_cdom(counts_output, counts_dark, scale_factor):
    """
    Description:

        The OOI Level 1 Fluorometric CDOM concentration core data product is a measure of how much
        light has been re-emitted from refractory colored organic compounds found in the colored
        dissolved organic matter (CDOM) in seawater. This data product describes as a measure of the
        amount of tannins (polyphenols that bind to proteins and other large molecules) or lignins
        (polymers of phenolic acids) from decaying plant material or byproducts from the decomposition
        of animals. It accounts for the tea-like color of some water masses. CDOM is not particulate, but
        water masses can contain both CDOM and turbidity. CDOM absorbs ultraviolet light and
        fluoresces visible blue light. The fluorescence of CDOM is used in many applications such as
        continuous monitoring of wastewater discharge, natural tracer of specific water bodies, ocean
        color research and the effect of CDOM on satellite imagery, and investigations of CDOM
        concentrations impacting light availability used for primary production.

    Implemented by:

        2014-01-30: Craig Risien. Initial Code
        
    Usage:

        cdom_conc = flo_cdom(counts_output, counts_dark, scale_factor)

            where

        cdom_conc = Fluorometric CDOM Concentration (CDOMFLO_L1) [ppb]
        counts_output = calibrated sample voltage output (CDOMFLO_L0) [counts]
        counts_dark = measured signal output of fluormeter in clean water with black tape over the detector [counts]
        scale_factor = multiplier [ppb / volts]

    References:
    
        OOI (2012). Data Product Specification for Fluorometric Chlorophyll-a Concentration. Document Control Number
        1341-00550. https://alfresco.oceanobservatories.org/
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00550_Data_Product_SPEC_CDOMFLO_OOI.pdf)
    """

    cdom_conc = flo_scale_and_offset(counts_output, counts_dark, scale_factor)
    return cdom_conc


#def flo_flubsct(counts_output, counts_dark, scale_factor):
#    """
#    Description:
#
#        The OOI Level 1 Optical backscatter (red wavelengths) core data product is an estimate of
#        turbidity and suspended solids in seawater that scatter photons of light in the back direction. Red
#        wavelengths of light fall between roughly 630 and 740nm. Turbidity commonly describes water
#        clarity and is a gross assessment of light attenuation factors like suspended solids, but not a
#        direct measurement of them, only their effect (Boss, et al, 2009). Optical backscatter meters
#        measure red light scattered from suspended matter which is a proxy for turbidity and suspended
#        solids. The size, composition and shape of the suspended particles affect the meters response,
#        so pre-deployment field verification is necessary to define a standard that adequately represents
#        the expected type and size of suspended matter found in situ is crucial to achieve the highest
#        quality data.
#
#    Implemented by:
#
#        2014-01-30: Craig Risien. Initial Code
#        
#    Usage:
#
#        backscat = flo_flubsct(counts_output, counts_dark, scale_factor)
#
#            where
#
#        opt_backscat = Optical backscatter (FLUBSCT_L1) [red wavelengths / meter / steradian]
#        counts_output = calibrated sample voltage output (FLUBSCT_L0) [counts]
#        counts_dark = measured signal output of fluormeter in clean water with black tape over the detector [counts]
#        scale_factor = multiplier
#
#    References:
#    
#        OOI (2012). Data Product Specification for Fluorometric Chlorophyll-a Concentration. Document Control Number
#        1341-00540. https://alfresco.oceanobservatories.org/
#        (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00540_Data_Product_SPEC_FLUBSCT_OOI.pdf)
#    """
#
#    backscat = flo_scale_and_offset(counts_output, counts_dark, scale_factor)
#    return backscat