#!/usr/bin/env python
"""
@package ion_functions.data.ctd_functions
@file ion_functions/data/ctd_functions.py
@author Christopher Wingard
@brief Module containing CTD related data-calculations.
"""

# Import Numexpr and the TEOS-10 GSW libraries
import numexpr
from pygsw import vectors as gsw


def ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3):
    """
    Description:

        OOI Level 1 Water Temperature data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 16Plus instruments and applies to
        CTDBP instruments, all series, and CTDPF instruments, series A and B.

    Implemented by:

        2013-04-12: Luke Campbell. Initial Code
        2013-04-12: Christopher Wingard. Minor edits
        2013-05-10: Christopher Wingard. Minor edits to comments.
        2014-01-31: Russell Desiderio. Standardized comment format.

    Usage:

        t = ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3)

            where

        t = sea water temperature (TEMPWAT_L1) [deg_C]
        t0 = raw temperature (TEMPWAT_L0) [counts]
        a0 = temperature calibration coefficients
        a1 = temperature calibration coefficients
        a2 = temperature calibration coefficients
        a3 = temperature calibration coefficients

    References:

        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)
    """

    mv = numexpr.evaluate('(t0 - 524288) / 1.6e7')
    r = numexpr.evaluate('(mv * 2.9e9 + 1.024e8)/(2.048e4 - mv * 2.0e5)')
    t = numexpr.evaluate('1 / (a0 + a1 * log(r) + a2 * log(r)**2 + a3 * log(r)**3) - 273.15')
    return t


def ctd_sbe37im_tempwat(t0):
    """
    Description:

        OOI Level 1 Water Temperature data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 37IM instruments and applies to
        CTDMO instruments, all series.

    Implemented by:

        2014-02-05: Russell Desiderio. Initial Code

    Usage:

        t = ctd_sbe37im_tempwat(t0)

            where

        t = sea water temperature (TEMPWAT_L1) [deg_C]
        t0 = raw temperature (TEMPWAT_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)
    """

    t = t0 / 10000.0 - 10.0
    return t


def ctd_sbe52mp_tempwat(t0):
    """
    Description:

        OOI Level 1 Water Temperature data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 52MP instruments and applies to
        CTDPF instruments, C,K,L series.

    Implemented by:

        2014-02-17: Russell Desiderio. Initial Code

    Usage:

        t = ctd_sbe52mp_tempwat(t0)

            where

        t = sea water temperature (TEMPWAT_L1) [deg_C]
        t0 = raw temperature (TEMPWAT_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Water Temperature. Document
            Control Number 1341-00010. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)
    """

    t = t0 / 10000.0 - 5.0
    return t


def ctd_sbe16plus_preswat(p0, t0, ptempa0, ptempa1, ptempa2,
                          ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                          pa0, pa1, pa2):
    """
    Description:

        OOI Level 1 Pressure (Depth) data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 16Plus instruments outfitted with
        a strain gauge pressure sensor. This is the default for most of the CTDBP
        instruments (the exceptions are series N and O) and for CTDPF instruments,
        series A and B.

    Implemented by:

        2013-04-12: Chris Wingard. Initial Code.
        2013-05-10: Christopher Wingard. Minor edits to comments.
        2014-01-31: Russell Desiderio. Standardized comment format.

    Usage:

        p = ctd_sbe16plus_preswat(p0, t0, ptempa0, ptempa1, ptempa2,
                          ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                          pa0, pa1, pa2)

            where

        p = sea water pressure (PRESWAT_L1) [dbar]
        p0 = raw pressure (PRESWAT_L0) [counts]
        t0 = raw temperature from pressure sensor thermistor [counts]
        ptempa0 = strain gauge pressure calibration coefficients
        ptempa1 = strain gauge pressure calibration coefficients
        ptempa2 = strain gauge pressure calibration coefficients
        ptca0 = strain gauge pressure calibration coefficients
        ptca1 = strain gauge pressure calibration coefficients
        ptca2 = strain gauge pressure calibration coefficients
        ptcb0 = strain gauge pressure calibration coefficients
        ptcb1 = strain gauge pressure calibration coefficients
        ptcb2 = strain gauge pressure calibration coefficients
        pa0 = strain gauge pressure calibration coefficients
        pa1 = strain gauge pressure calibration coefficients
        pa2 = strain gauge pressure calibration coefficients

    References:

        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
    """
    # compute calibration parameters
    tv = t0 / 13107.0
    t = ptempa0 + ptempa1 * tv + ptempa2 * tv**2
    x = p0 - ptca0 - ptca1 * t - ptca2 * t**2
    n = x * ptcb0 / (ptcb0 + ptcb1 * t + ptcb2 * t**2)

    # compute pressure in psi, rescale and compute in dbar and return
    p_psi = pa0 + pa1 * n + pa2 * n**2
    p_dbar = (p_psi * 0.689475729) - 10.1325
    return p_dbar


def ctd_sbe16digi_preswat(p0, t0, C1, C2, C3, D1, D2, T1, T2, T3, T4, T5):
    """
    Description:

        OOI Level 1 Pressure (Depth) data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 16Plus instruments outfitted with
        a digiquartz pressure sensor. This applies to the CTDBP-N,O instruments
        only.

    Implemented by:

        2013-05-10: Christopher Wingard. Initial Code.
        2013-05-10: Christopher Wingard. Minor edits to comments.
        2014-01-31: Russell Desiderio. Standardized comment format.
        2014-01-31: Russell Desiderio. Modified algorithm to use pressure [Hz] (pf)
                    to calculate pressure period instead of pressure [counts] (p0).
                    See SeaBird 16Plus V2 User Manual (reference (2)), page 57, item 5.

    Usage:

        p = ctd_sbe16digi_preswat(p0,t0,C1,C2,C3,D1,D2,T1,T2,T3,T4,T5)

            where

        p = sea water pressure (PRESWAT_L1) [dbar]
        p0 = raw pressure (PRESWAT_L0) [counts]
        t0 = raw temperature from pressure sensor thermistor [counts]
        C1 = digiquartz pressure calibration coefficients
        C2 = digiquartz pressure calibration coefficients
        C3 = digiquartz pressure calibration coefficients
        D1 = digiquartz pressure calibration coefficients
        D2 = digiquartz pressure calibration coefficients
        T1 = digiquartz pressure calibration coefficients
        T2 = digiquartz pressure calibration coefficients
        T3 = digiquartz pressure calibration coefficients
        T4 = digiquartz pressure calibration coefficients
        T5 = digiquartz pressure calibration coefficients

    References:

        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)

        OOI (2011). SeaBird 16Plus V2 User Manual. 1341-00020_PRESWAT Artifact.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >>
            >> REFERENCE >> Data Product Specification Artifacts >> 1341-00020_PRESWAT >>
            PRESWAT_SeaBird_16PlusV2_2009.pdf)
    """
    # Convert raw pressure input to frequency [Hz]
    pf = p0 / 256.0

    # Convert raw temperature input to voltage
    tv = t0 / 13107.0

    # Calculate U (thermistor temp):
    U = (23.7 * (tv + 9.7917)) - 273.15

    # Calculate calibration parameters
    C = C1 + C2 * U + C3 * U**2
    D = D1 + D2 * U
    T0 = T1 + T2 * U + T3 * U**2 + T4 * U**3 + T5 * U**4

    # Calculate T (pressure period, in microseconds):
    T = (1.0 / pf) * 1.0e6

    # compute pressure in psi, rescale and compute in dbar and return
    p_psi = C * (1.0 - T0**2 / T**2) * (1.0 - D * (1.0 - T0**2 / T**2))
    p_dbar = (p_psi * 0.689475729) - 10.1325
    return p_dbar


def ctd_sbe37im_preswat(p0, p_range_psia):
    """
    Description:

        OOI Level 1 Pressure (Depth) data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 37IM instruments and applies to
        CTDMO instruments, all series.

    Implemented by:

        2014-02-05: Russell Desiderio. Initial Code

    Usage:

        p = ctd_sbe37im_preswat(p0, p_range_psia)

            where

        p = sea water pressure (PRESWAT_L1) [dbar]
        p0 = raw pressure (PRESWAT_L0) [counts]
        p_range_psia = pressure range calibration coefficient [psia]

    References:

        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
    """
    # compute pressure range in units of dbar
    p_range_dbar = (p_range_psia - 14.7) * 0.6894757

    # compute pressure in dbar and return
    p_dbar = p0 * p_range_dbar / (0.85 * 65536.0) - 0.05 * p_range_dbar
    return p_dbar


def ctd_glider_preswat(pr_bar):
    """
    Description:

        OOI Level 1 Pressure (Depth) data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from Seabird CTDs installed on gliders and
        applies to CTDGV instruments.

    Implemented by:

        2015-10-28: Russell Desiderio. Initial Code

    Usage:

        pr_dbar = ctd_glider_preswat(pr_bar)

            where

        pr_dbar = sea water pressure (PRESWAT_L1) [dbar]
        pr_bar = sea water pressure value from glider [bar]

    References:

        OOI (2015). Data Product Specification for Coastal Glider Data Products
            (version 1-03). Document Control Number 1341-00020.
            https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-20001_Data_Product_SPEC_CSGLIDR_OOI.docx)
    """

    pr_dbar = pr_bar * 10.0
    return pr_dbar


def ctd_sbe52mp_preswat(p0):
    """
    Description:

        OOI Level 1 Pressure (Depth) data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 52MP instruments and applies to
        CTDPF instruments, C,K,L series.

    Implemented by:

        2014-02-17: Russell Desiderio. Initial Code

    Usage:

        p = ctd_sbe52mp_preswat(p0)

            where

        p = sea water pressure (PRESWAT_L1) [dbar]
        p0 = raw pressure (PRESWAT_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Pressure (Depth). Document
            Control Number 1341-00020. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
    """

    p_dbar = p0 / 100.0 - 10.0
    return p_dbar


def ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor):
    """
    Description:

        OOI Level 1 Conductivity core data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 16Plus instruments and applies to
        CTDBP instruments, all series, and CTDPF instruments, series A and B.

    Implemented by:

        2013-04-12: Christopher Wingard. Initial Code
        2013-05-10: Christopher Wingard. Minor edits to comments.
        2014-01-31: Russell Desiderio. Standardized comment format.

    Usage:

        c = ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor)

            where

        c = sea water conductivity (CONDWAT_L1) [S m-1]
        c0 = sea water conductivity (CONDWAT_L0) [counts]
        t1 = sea water temperature (TEMPWAT_L1) [deg_C]
        p1 = sea water pressure (PRESWAT_L1) [dbar]
        g = conductivity calibration coefficients
        h = conductivity calibration coefficients
        i = conductivity calibration coefficients
        j = conductivity calibration coefficients
        cpcor = conductivity calibration coefficients
        ctcor = conductivity calibration coefficients

    References:

        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
    """
    # convert raw conductivity measurement to frequency
    f = (c0 / 256.0) / 1000.0

    # calculate conductivity [S m-1]
    c = (g + h * f**2 + i * f**3 + j * f**4) / (1 + ctcor * t1 + cpcor * p1)
    return c


def ctd_sbe37im_condwat(c0):
    """
    Description:

        OOI Level 1 Conductivity core data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 37IM instruments and applies to
        CTDMO instruments, all series.

    Implemented by:

        2014-02-05: Russell Desiderio. Initial Code

    Usage:

        c = ctd_sbe37im_condwat(c0)

            where

        c = sea water conductivity (CONDWAT_L1) [S m-1]
        c0 = sea water conductivity (CONDWAT_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
   """

    c = c0 / 100000.0 - 0.5
    return c


def ctd_sbe52mp_condwat(c0):
    """
    Description:

        OOI Level 1 Conductivity core data product, which is calculated using
        data from the Sea-Bird Electronics conductivity, temperature and depth
        (CTD) family of instruments.

        This data product is derived from SBE 52MP instruments and applies to
        CTDPF instruments, C,K,L series.

    Implemented by:

        2014-02-17: Russell Desiderio. Initial Code

    Usage:

        c = ctd_sbe52mp_condwat(c0)

            where

        c = sea water conductivity (CONDWAT_L1) [S m-1]
        c0 = sea water conductivity (CONDWAT_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Conductivity. Document
            Control Number 1341-00030. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
   """

    c_mmho_cm = c0 / 10000.0 - 0.5
    c_S_m = 0.1 * c_mmho_cm
    return c_S_m


def ctd_pracsal(c, t, p):
    """
    Description:

        OOI Level 2 Practical Salinity core data product, which is calculated
        using the Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version
        3.0, with data from the conductivity, temperature and depth (CTD)
        family of instruments.

    Implemented by:

        2013-03-13: Christopher Wingard. Initial code.
        2013-05-10: Christopher Wingard. Minor edits to comments.
        2014-01-31: Russell Desiderio. Standardized comment format.

    Usage:

        SP = ctd_pracsal(c, t, p)

            where

        SP = practical salinity, PSS-78, (PRACSAL_L2) [unitless]
        c = sea water conductivity (CONDWAT_L1) [S m-1]
        t = sea water temperature (TEMPWAT_L1) [deg_C]
        p = sea water pressure (PRESWAT_L1) [dbar]
    References:

        OOI (2012). Data Product Specification for Salinty. Document Control
            Number 1341-00040. https://alfresco.oceanobservatories.org/ (See: 
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf)
    """

    # Convert L1 Conductivity from S/m to mS/cm
    C10 = c * 10.0

    # Calculate the Practical Salinity (PSS-78) [unitless]
    SP = gsw.sp_from_c(C10, t, p)
    return SP


def ctd_density(SP, t, p, lat, lon):
    """
    Description:

        OOI Level 2 Density core data product, which is calculated using the
        Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version 3.0, with
        data from the conductivity, temperature and depth (CTD) family of
        instruments.

    Implemented by:

        2013-03-11: Christopher Mueller. Initial code.
        2013-03-13: Christopher Wingard. Added commenting and moved to
            ctd_functions
        2013-05-10: Christopher Wingard. Minor edits to comments.
        2014-01-31: Russell Desiderio. Standardized comment format.

    Usage:

        rho = ctd_density(SP, t, p, lat, lon)

            where

        rho = sea water density (DENSITY_L2) [kg m-3]
        SP = practical salinity PSS-78 (PRACSAL_L2) [unitless]
        t = sea water temperature (TEMPWAT_L1) [deg_C]
        p = sea water pressure (PRESWAT_L1) [dbar]
        lat = latitude where input data was collected [decimal degree]
        lon = longitude where input data was collected [decimal degree]

    References:

        OOI (2012). Data Product Specification for Density. Document Control
            Number 1341-00050. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf)
    """
    # Calculate the density [kg m-3]
    rho = gsw.ctd_density(SP, t, p, lat, lon)
    return rho
