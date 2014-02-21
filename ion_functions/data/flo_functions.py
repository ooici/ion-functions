#!/usr/bin/env python
"""
@package ion_functions.data.flo_functions
@file ion_functions/data/flo_functions.py
<<<<<<< HEAD
@author Christopher Wingard
@brief Module containing Fluorometer Three Wavelength (FLORT) and Fluorometer
    Two Wavelength (FLORD) instrument family related functions
"""
import numpy as np
import numexpr as ne


def flo_bback_total(beta, degC=20.0, psu=32.0, theta=124.0, wlngth=700.0,
                    xfactor=1.08):
    """
    Description:

        OOI Level 1 Optical Backscatter (Red Wavelengths, 700 nm with a 117
        degree scattering angle) data product (FLUBSCT), which is calculated
        using data from the WET Labs, Inc. ECO fluorometer family of
        instruments.

    Implemented by:

        2013-07-16: Christopher Wingard. Initial Code

    Usage:

        bback = flo_bback_total(beta, degC, psu, theta, wlngth, xfactor)

            where

        bback = total optical backscatter (FLUBSCT) [m-1]
        beta = total volume scattering coefficient [m-1 sr-1]
        degC = in situ water temperature from co-located CTD, if available [deg_C].
            Default is 20.
        psu = in situ salinity from a co-located CTD, if available [psu].
            Default is 32.
        wlngth = optical backscatter measurement wavelength [nm]. All OOI FLORT
            and FLORD instruments use 700 nm.
        theta = optical backscatter scattering angle [degrees]. All OOI FLORT
            and FLORD instruments use 117 degrees.
        xfactor = X (Chi) factor for high angular resolution. For 117 degree
            scattering angle X = 1.10.

    References:

        OOI (2012). Data Product Specification for Optical Backscatter (Red
            Wavelengths). Document Control Number 1341-00540.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00540_Data_Product_SPEC_FLUBSCT_OOI.pdf)
    """
    # calculate the volume scattering coefficient (m-1 sr-1, betasw) and the
    # optical backscatter (m-1, bbacksw) of seawater using data from a
    # co-located CTD. Values below are computed using provided code from Zhang
    # et al 2009.
    #betasw, bsw = flo_zhang_scatter_coeffs(degC, psu, theta, wlngth)

    # recompute bsw from equation in the the manual...
    betasw = 1.38 * (wlngth / 500.0)**-4.32 * (1 + 0.3 * psu / 37) * 10**-4 \
        * (1 + ((1 - 0.09)/(1 + 0.09)) * np.cos(np.radians(theta))**2)
    bsw = 0.0029308 * (wlngth / 500.0)**-4.24

    # calculate the volume scattering of particles (total - seawater)
    betap = beta - betasw

    # calculate the particulate backscatter coefficient (m-1)
    bbackp = xfactor * 2.0 * np.pi * betap

    # calculate the total backscatter coefficient (particulates + seawater)
    bback = bbackp + (bsw / 2.0)

    print betasw, bsw, beta, betap, bbackp, bback
    return bback


def flo_zhang_scatter_coeffs(degC, psu, theta=124.0, wlngth=700.0, delta=0.039):
    """
    Description:

        Computes the scattering coefficients of seawater (both the volume
        scattering and the total backscatter of seawater) based on computation
        of Zhang et al 2009 as presented in the DPS for optical backscatter
        (red wavelengths).

        This code is derived from Matlab code developed and made available
        online by:

            Dr. Xiaodong Zhang
            Associate Professor
            Department of Earth Systems Science and Policy
            University of North Dakota 
            http://www.und.edu/instruct/zhang/

    Implemented by:

        2013-07-15: Christopher Wingard. Initial Code

    Usage:

        betasw, bsw = flo_zhang_scatter_coeffs(degC, psu, theta, wlngth, delta)

            where

        betasw = volume scattering coefficient of pure seawater [m-1 sr-1]
        bsw = total scattering coefficient of pure seawater [m-1]
        degC = in situ water temperature from co-located CTD [deg_C]
        psu = in situ salinity from a co-located CTD [psu]
        theta = optical backscatter angle [degrees]. All OOI FLORT
            and FLORD instruments use 117 degrees.
        wlngth = optical backscatter measurement wavelength [nm]. All OOI FLORT
            and FLORD instruments use 700 nm.
        delta = depolarization ratio [unitless]. Default of 0.039 is assumed.

    References:

        OOI (2012). Data Product Specification for Optical Backscatter (Red
            Wavelengths). Document Control Number 1341-00540.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00540_Data_Product_SPEC_FLUBSCT_OOI.pdf)
    """
    # values of the constants
    Na = 6.0221417930e23  # Avogadro's constant
    Kbz = 1.3806503e-23  # Boltzmann constant
    degK = degC + 273.15  # Absolute tempearture
    M0 = 0.018  # Molecular weigth of water in kg/mol
    pi = np.pi

    # convert the scattering angle from degrees to radians
    rad = np.radians(theta)

    # calculate the absolute refractive index of seawater and the partial
    # derivative of seawater refractive index with regards to salinity.
    nsw, dnds = flo_refractive_index(wlngth, degC, psu)

    # isothermal compressibility is from Lepple & Millero (1971,Deep
    # Sea-Research), pages 10-11 The error ~ +/-0.004e-6 bar^-1
    icomp = flo_isotherm_compress(degC, psu)

    # density of seawater from UNESCO 38 (1981). Note, this could/should be
    # changed to use the Gibbs Seawater Toolbox routines built in to
    # ion-functions.
    rho = flo_density_seawater(degC, psu)

    # water activity data of seawater is from Millero and Leung (1976, American
    # Journal of Science, 276, 1035-1077). Table 19 was reproduced using
    # Eq.(14,22,23,88,107) that were fitted to polynominal equation.
    # dlnawds is partial derivative of the natural logarithm of water activity
    # with regards to salinity.
    dlnawds = ne.evaluate('(-5.58651e-4 + 2.40452e-7 * degC - 3.12165e-9 * degC**2 + 2.40808e-11 * degC**3) +'
                          '1.5 * (1.79613e-5 - 9.9422e-8 * degC + 2.08919e-9 * degC**2 - 1.39872e-11 * degC**3) *'
                          'psu**0.5 + 2 * (-2.31065e-6 - 1.37674e-9 * degC - 1.93316e-11 * degC**2) * psu')

    # density derivative of refractive index from PMH model
    dfri = ne.evaluate('(nsw**2 - 1.0) * (1.0 + 2.0/3.0 * (nsw**2 + 2.0)'
                       '* (nsw/3.0 - 1.0/3.0 / nsw)**2)')

    # volume scattering at 90 degree due to the density fluctuation
    beta_df = ne.evaluate('pi**2 / 2.0 * (wlngth*1e-9)**-4 * Kbz * degK * icomp '
                          '* dfri**2 * (6.0 + 6.0 * delta) / (6.0 - 7.0 * delta)')

    # volume scattering at 90 degree due to the concentration fluctuation
    flu_con = ne.evaluate('psu * M0 * dnds**2 / rho / -dlnawds / Na')
    beta_cf = ne.evaluate('2.0 * pi**2 * (wlngth * 1e-9)**-4 * nsw**2 * flu_con'
                          '* (6.0 + 6.0 * delta) / (6.0 - 7.0 * delta)')

    # total volume scattering at 90 degree
    beta90sw = beta_df + beta_cf

    # total scattering coefficient of seawater (m-1)
    bsw = ne.evaluate('8.0 * pi / 3.0 * beta90sw * ((2.0 + delta) / (1.0 + delta))')

    # total volume scattering coefficient of seawater (m-1 sr-1)
    betasw = ne.evaluate('beta90sw * (1.0 + ((1.0 - delta) / (1.0 + delta)) * cos(rad)**2)')

    return betasw, bsw


def flo_refractive_index(wlngth, degC, psu):
    """
    Helper function for flo_zhang_scatter_coeffs
    @param wlngth backscatter measurement wavlength (nm)
    @param degC in situ water temperature (deg_C)
    @param psu in site practical salinity (psu)
    @retval nsw absolute refractive index of seawater
    @retval dnds partial derivative of seawater refractive index with regards to
        seawater.
    """
    # refractive index of air is from Ciddor (1996, Applied Optics).
    n_air = ne.evaluate('1.0 + (5792105.0 / (238.0185 - 1 / (wlngth/1e3)**2)'
                        '+ 167917.0 / (57.362 - 1 / (wlngth/1e3)**2)) / 1e8')

    # refractive index of seawater is from Quan and Fry (1994, Applied Optics)
    n0 = 1.31405
    n1 = 1.779e-4
    n2 = -1.05e-6
    n3 = 1.6e-8
    n4 = -2.02e-6
    n5 = 15.868
    n6 = 0.01155
    n7 = -0.00423
    n8 = -4382.0
    n9 = 1.1455e6
    nsw = ne.evaluate('n0 + (n1 + n2 * degC + n3 * degC**2) * psu + n4 * degC**2'
                      '+ (n5 + n6 * psu + n7 * degC) / wlngth + n8 / wlngth**2'
                      '+ n9 / wlngth**3')

    # pure seawater
    nsw = ne.evaluate('nsw * n_air')
    dnds = ne.evaluate('(n1 + n2 * degC + n3 * degC**2 + n6 / wlngth) * n_air')

    return nsw, dnds


def flo_isotherm_compress(degC, psu):
    """
    Helper function for flo_zhang_scatter_coeffs
    @param degC in situ water temperature
    @param psu in site practical salinity
    @retval iso_comp seawater isothermal compressibility
    """
    # pure water secant bulk Millero (1980, Deep-sea Research)
    kw = ne.evaluate('19652.21 + 148.4206 * degC - 2.327105 * degC**2'
                     '+ 1.360477e-2 * degC**3 - 5.155288e-5 * degC**4')

    # seawater secant bulk
    a0 = ne.evaluate('54.6746 - 0.603459 * degC + 1.09987e-2 * degC**2'
                     '- 6.167e-5 * degC**3')
    b0 = ne.evaluate('7.944e-2 + 1.6483e-2 * degC - 5.3009e-4 * degC**2')
    ks = ne.evaluate('kw + a0 * psu + b0 * psu**1.5')

    # calculate seawater isothermal compressibility from the secant bulk
    iso_comp = ne.evaluate('1 / ks * 1e-5')  # unit is Pa

    return iso_comp


def flo_density_seawater(degC, psu):
    """
    Helper function for flo_zhang_scatter_coeffs
    @param degC in situ water temperature
    @param psu in site practical salinity
    @retval rho_sw density of seawater
    """
    # density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
    a0 = 8.24493e-1
    a1 = -4.0899e-3
    a2 = 7.6438e-5
    a3 = -8.2467e-7
    a4 = 5.3875e-9
    a5 = -5.72466e-3
    a6 = 1.0227e-4
    a7 = -1.6546e-6
    a8 = 4.8314e-4
    b0 = 999.842594
    b1 = 6.793952e-2
    b2 = -9.09529e-3
    b3 = 1.001685e-4
    b4 = -1.120083e-6
    b5 = 6.536332e-9

    # density for pure water
    rho_w = ne.evaluate('b0 + b1 * degC + b2 * degC**2 + b3 * degC**3'
                        '+ b4 * degC**4 + b5 * degC**5')

    # density for pure seawater
    rho_sw = ne.evaluate('rho_w + ((a0 + a1 * degC + a2 * degC**2'
                         '+ a3 * degC**3 + a4 * degC**4) * psu'
                         '+ (a5 + a6 * degC + a7 * degC**2) * psu**1.5 + a8 * psu**2)')

    return rho_sw


=======
@author Craig Risien
@brief Module containing calculations related to instruments in the fluorescence family.
"""

# Import NumExpr

import numexpr as ne
    
    
>>>>>>> 78e66739c1350fc4ac9fa993a4aed99e0562b92d
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

<<<<<<< HEAD
=======

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
>>>>>>> 78e66739c1350fc4ac9fa993a4aed99e0562b92d
