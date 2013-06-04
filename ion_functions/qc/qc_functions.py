#!/usr/bin/env python

"""
@package ion_functions.qc_functions
@file ion_functions/qc_functions.py
@author Christopher Mueller
@brief Module containing QC functions ported from matlab samples in DPS documents
"""
from ion_functions.qc.qc_extensions import stuckvalues, spikevalues

import time
import numpy as np
import numexpr as ne
from scipy.interpolate import LinearNDInterpolator
from ion_functions import utils

# try to load the OOI logging module, using default Python logging module if
# unavailable
try:
    from ooi.logging import log
except ImportError:
    from logging import log


def dataqc_globalrangetest_minmax(dat, dat_min, dat_max, strict_validation=False):
    '''
    Python wrapper for dataqc_globalrangetest
    Combines the min/max arguments into list for dataqc_globalrangetest
    '''
    return dataqc_globalrangetest(dat, [dat_min, dat_max], strict_validation=strict_validation)


def dataqc_globalrangetest(dat, datlim, strict_validation=False):
    """
    Description:

        Data quality control algorithm testing if measurements fall into a
        user-defined valid range. Returns 1 for presumably good data and 0 for
        data presumed bad.

    Implemented by:

        2010-07-28: DPS authored by Mathias Lankhorst. Example code provided
        for Matlab.
        2013-04-06: Christopher Wingard. Initial python implementation.
        2013-05-30: Christopher Mueller. Performance improvements by adding
            strict_validation flag.

    Usage:

        qcflag = dataqc_globalrangetest(dat, datlim, strict_validation)

            where

        qcflag = Boolean, 0 if value is outside range, else = 1.

        dat = Input dataset, any scalar or vector. Must be numeric and real.
        datlim = Two-element vector with the minimum and maximum values
            considered to be valid.
        strict_validation = Flag (default is False) to assert testing of input
            types (e.g. isreal, isnumeric)

    References:

        OOI (2012). Data Product Specification for Global Range Test. Document
            Control Number 1341-10004. https://alfresco.oceanobservatories.org
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-10004_Data_Product_SPEC_GLBLRNG_OOI.pdf)
    """
    dat = np.atleast_1d(dat)
    datlim = np.atleast_1d(datlim)

    if strict_validation:
        if not utils.isnumeric(dat).all():
            raise ValueError('\'dat\' must be numeric')

        if not utils.isreal(dat).all():
            raise ValueError('\'dat\' must be real')

        if not utils.isnumeric(datlim).all():
            raise ValueError('\'datlim\' must be numeric')

        if not utils.isreal(datlim).all():
            raise ValueError('\'datlim\' must be real')

        if len(datlim) < 2:  # Must have at least 2 elements
            raise ValueError('\'datlim\' must have at least 2 elements')

    return (datlim.min() <= dat) & (dat <= datlim.max()).astype('int8')


def dataqc_localrangetest(dat, z, datlim, datlimz, strict_validation=False):
    """
    Description:

        Data quality control algorithm testing if measurements fall into a
        user-defined valid range. This range is not constant but varies with
        measurement location. Returns 1 for presumably good data and 0 for data
        presumed bad.

    Implemented by:

        2012-07-17: DPS authored by Mathias Lankhorst. Example code provided
        for Matlab.
        2013-04-06: Christopher Wingard. Initial python implementation.

    Usage:

        qcflag = dataqc_localrangetest(dat, z, datlim, datlimz)

            where

        qcflag = Boolean, 0 if value is outside range, else = 1.

        dat = input data set, a numeric real scalar or column vector.
        z = location of measurement dat. must have same # of rows as dat and
            same # of columns as datlimz
        datlim = two column array with the minimum (column 1) and maximum
            (column 2) values considered valid.
        datlimz = array with the locations where datlim is given. must have
            same # of rows as datlim and same # of columns as z.

    References:

        OOI (2012). Data Product Specification for Local Range Test. Document
            Control Number 1341-10005. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-10005_Data_Product_SPEC_LOCLRNG_OOI.pdf)
    """
    if strict_validation:
        # check if dat and datlim are matrices
        if not utils.isvector(dat):
            raise ValueError('\'dat\' must be a matrix')

        if not utils.ismatrix(datlim):
            raise ValueError('\'datlim\' must be a matrix')

        # check if all inputs are numeric and real
        for k, arg in {'dat': dat, 'z': z, 'datlim': datlim,
                       'datlimz': datlimz}.iteritems():
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))

    # test size and shape of the input arrays datlimz and datlim, setting test
    # variables.
    array_size = datlimz.shape
    if len(array_size) == 1:
        numlim = array_size[0]
        ndim = 1
    else:
        numlim = array_size[0]
        ndim = array_size[1]

    array_size = datlim.shape
    tmp1 = array_size[0]
    tmp2 = array_size[1]
    if tmp1 != numlim:
        raise ValueError('\'datlim\' and \'datlimz\' must '
                         'have the same number of rows.')

    if tmp2 != 2:
        raise ValueError('\'datlim\' must be structured as 2-D array '
                         'with exactly 2 columns and 1 through N rows.')

    # test the size and shape of the z input array
    array_size = z.shape
    if len(array_size) == 1:
        num = array_size[0]
        tmp2 = 1
    else:
        num = array_size[0]
        tmp2 = array_size[1]

    if tmp2 != ndim:
        raise ValueError('\'z\' must have the same number of columns '
                         'as \'datlimz\'.')

    if num != dat.size:
        raise ValueError('Len of \'dat\' must match number of '
                         'rows in \'z\'')

    # test datlim, values in column 2 must be greater than those in column 1
    if not all(datlim[:, 1] > datlim[:, 0]):
        raise ValueError('Second column values of \'datlim\' should be '
                         'greater than first column values.')

    # calculate the upper and lower limits for the data set
    if ndim == 1:
        # determine the lower limits using linear interpolation
        lim1 = np.interp(z, datlimz, datlim[:, 0], left=np.nan, right=np.nan)
        # determine the upper limits using linear interpolation
        lim2 = np.interp(z, datlimz, datlim[:, 1], left=np.nan, right=np.nan)
    else:
        # Compute Delaunay Triangulation and use linear interpolation to
        # determine the N-dimensional lower limits
        F = LinearNDInterpolator(z, dat)
        points = np.reshape(np.array([datlimz, datlim[:, 0]]),
                            [numlim, ndim+1])
        lim1 = F(points)
        # Compute Delaunay Triangulation and use linear interpolation to
        # determine the N-dimensional upper limits
        F = LinearNDInterpolator(z, dat)
        points = np.reshape(np.array([datlimz, datlim[:, 1]]),
                            [numlim, ndim+1])
        lim2 = F(points)

    # replace NaNs from above interpolations
    ff = (np.isnan(lim1)) | (np.isnan(lim2))
    lim1[ff] = np.max(datlim[:, 1])
    lim2[ff] = np.min(datlim[:, 0])

    # compute the qcflags
    qcflag = (dat >= lim1) & (dat <= lim2)
    return qcflag.astype('int8')


def dataqc_spiketest(dat, acc, N=5, L=5, strict_validation=False):
    """
    Description:

        Data quality control algorithm testing a time series for spikes.
        Returns 1 for presumably good data and 0 for data presumed bad.

        The time series is divided into windows of len L (an odd integer
        number). Then, window by window, each value is compared to its (L-1)
        neighboring values: a range R of these (L-1) values is computed (max.
        minus min.), and replaced with the measurement accuracy ACC if ACC>R. A
        value is presumed to be good, i.e. no spike, if it deviates from the
        mean of the (L-1) peers by less than a multiple of the range,
        N*max(R,ACC).

        Further than (L-1)/2 values from the start or end points, the peer
        values are symmetrically before and after the test value. Within that
        range of the start and end, the peers are the first/last L values
        (without the test value itself).

        The purpose of ACC is to restrict spike detection to deviations
        exceeding a minimum threshold value (N*ACC) even if the data have
        little variability. Use ACC=0 to disable this behavior.

    Implemented by:

        2012-07-28: DPS authored by Mathias Lankhorst. Example code provided
        for Matlab.
        2013-04-06: Christopher Wingard. Initial python implementation.
        2013-05-30: Christopher Mueller. Performance optimizations.

    Usage:

        qcflag = dataqc_spiketest(dat, acc, N, L)

            where

        qcflag = Boolean, 0 if value is outside range, else = 1.

        dat = input data set, a numeric, real vector.
        acc = Accuracy of any input measurement.
        N = (optional, defaults to 5) Range multiplier, cf. above
        L = (optional, defaults to 5) Window len, cf. above

    References:

        OOI (2012). Data Product Specification for Spike Test. Document
            Control Number 1341-10006. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-10006_Data_Product_SPEC_SPKETST_OOI.pdf)
    """
    dat = np.atleast_1d(dat)

    if strict_validation:
        if not utils.isnumeric(dat).all():
            raise ValueError('\'dat\' must be numeric')

        if not utils.isreal(dat).all():
            raise ValueError('\'dat\' must be real')

        if not utils.isvector(dat):
            raise ValueError('\'dat\' must be a vector')

        for k, arg in {'acc': acc, 'N': N, 'L': L}.iteritems():
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))
    dat = np.asanyarray(dat, dtype=np.float)
    
    out = spikevalues(dat, L, N, acc)
    return out


def dataqc_polytrendtest(dat, t, ord_n=1, nstd=3, strict_validation=False):
    """
    Description:

        Data quality control algorithm testing if measurements contain a
        significant portion of a polynomial. Returns 1 if this is not the case,
        else 0.

        The purpose of this test is to check if a significant fraction of the
        variability in a time series can be explained by a drift, possibly
        interpreted as a sensor drift. This drift is assumed to be a polynomial
        of order ORD. Use ORD=1 to consider a linear drift

        The time series dat is passed to MatLab's POLYFIT routine to obtain a
        polynomial fit PP to dat, and the difference dat-PP is compared to the
        original dat. If the standard deviation of (dat-PP) is less than that
        of dat by a factor of NSTD, the time series is assumed to contain a
        significant trend (output will be 0), else not (output will be 1).

    Implemented by:

        2012-10-29: DPS authored by Mathias Lankhorst. Example code provided
        for Matlab.
        2013-04-06: Christopher Wingard. Initial python implementation.
        2013-05-30: Christopher Mueller. Performance optimizations.

    Usage:

        qcflag = dataqc_polytrendtest(dat, t, ord_n, nstd, strict_validation)

            where

        qcflag = Boolean, 0 a trend is detected, 1 elsewhere.
        dat = Input dataset, a numeric real vector.
        t = time record associated with dat
        ord_n (optional, defaults to 1) = Polynomial order.
        nstd (optional, defaults to 3) = Factor by how much the standard
            deviation must be reduced before qcflag switches from 1 to 0
        strict_validation (optional, defaults to False) = Flag asserting
            testing of inputs.

    References:

        OOI (2012). Data Product Specification for Trend Test. Document
            Control Number 1341-10007. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-10007_Data_Product_SPEC_TRNDTST_OOI.pdf)
    """
    dat = np.atleast_1d(dat)
    t = np.atleast_1d(t)

    if strict_validation:
        for k, arg in {'dat': dat, 't': t, 'ord_n': ord_n, 'nstd': nstd}.iteritems():
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))

        for k, arg in {'dat': dat, 't': t}.iteritems():
            if not utils.isvector(arg):
                raise ValueError('\'{0}\' must be a vector'.format(k))

        for k, arg in {'ord_n': ord_n, 'nstd': nstd}.iteritems():
            if not utils.isscalar(arg):
                raise ValueError('\'{0}\' must be a scalar'.format(k))

    ord_n = int(round(abs(ord_n)))
    nstd = int(abs(nstd))

    ll = len(dat)
    # Not needed because time is incorporated as 't'
    # t = range(ll)

    pp = np.polyfit(t, dat, ord_n)
    datpp = np.polyval(pp, t)

    # test for a trend
    if np.atleast_1d((np.std(dat - datpp) * nstd) < np.std(dat)).all():
        trndtst = 0
    else:
        trndtst = 1

    # insure output size equals input, even though test yields a single value.
    qcflag = np.ones(ll).astype('int8') * trndtst
    return qcflag


def dataqc_stuckvaluetest(x, reso, num=10, strict_validation=False):
    """
    Description:

        Data quality control algorithm testing a time series for "stuck
        values", i.e. repeated occurences of one value. Returns 1 for
        presumably good data and 0 for data presumed bad.

    Implemented by:

        2012-10-29: DPS authored by Mathias Lankhorst. Example code provided
        for Matlab.
        2013-04-06: Christopher Wingard. Initial python implementation.

    Usage:

        qcflag = =dataqc_stuckvaluetest(x, RESO, NUM);

            where

        qcflag = Boolean output: 0 where stuck values are found, 1 elsewhere.
        x = Input time series (vector, numeric).
        reso = Resolution; repeat values less than reso apart will be
            considered "stuck values".
        num = Minimum number of successive values within reso of each other
            that will trigger the "stuck value". num is optional and defaults
            to 10 if omitted or empty.

    References:

        OOI (2012). Data Product Specification for Stuck Value Test. Document
            Control Number 1341-10008. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-10008_Data_Product_SPEC_STUCKVL_OOI.pdf)
    """
    dat = np.atleast_1d(x)

    if strict_validation:
        if not utils.isnumeric(dat).all():
            raise ValueError('\'x\' must be numeric')

        if not utils.isvector(dat):
            raise ValueError('\'x\' must be a vector')

        if not utils.isreal(dat).all():
            raise ValueError('\'x\' must be real')

        for k, arg in {'reso': reso, 'num': num}.iteritems():
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isscalar(arg):
                raise ValueError('\'{0}\' must be a scalar'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))

    num = np.abs(num)
    dat = np.asanyarray(dat, dtype=np.float)
    ll = len(x)
    if ll < num:
        # Warn - 'num' is greater than len(x), returning zeros
        out = np.zeros(dat.size, dtype='int8')
    else:
        out = stuckvalues(dat, reso, num)

    return out


def dataqc_gradienttest(dat, x, ddatdx, mindx, startdat, toldat, strict_validation=False):
    """
    Description

        Data quality control algorithm testing if changes between successive
        data points fall within a certain range.

        Input data dat are given as a function of coordinate x. The algorithm
        will flag dat values as bad if the change deltaDAT/deltaX between
        successive dat values exceeds thresholds given in ddatdx. Once the
        threshold is exceeded, following dat are considered bad until a dat
        value returns to within toldat of the last known good value.

        It is possible to remove data points that are too close together in x
        coordinates (use mindx).

        By default, the first value of dat is considered good. To change this,
        use startdat and toldat to set as the first good data point the first
        one that comes within toldat of startdat.

    Implemented by:

        2012-07-17: DPS authored by Mathias Lankhorst. Example code provided
        for Matlab.
        2013-04-06: Christopher Wingard. Initial python implementation.

    Usage:

        outdat, outx, outqc = dataqc_gradienttest(dat, x, ddatdx, mindx,
                                                  startdat, toldat);

            where

        outdat = same as dat except that NaNs and values not meeting mindx are
            removed.
        outx = same as x except that NaNs and values not meeting mindx are
            removed.
        outqc = output quality control flags for outdat. 0 means bad data, 1
            means good data.

        dat = input dataset, a numeric real vector.
        x = coordinate (e.g. time, distance) along which dat is given. Must be
            of the same size as dat and strictly increasing.
        ddatdx = two-element vector defining the valid range of ddat/dx
            from one point to the next.
        mindx = scalar. minimum dx for which this test will be applied (data
            that are less than mindx apart will be deleted). defaults to zero
            if NaN/empty.
        startdat = start value (scalar) of dat that is presumed good. defaults
            to first non-NaN value of dat if NaN/empty.
        toldat = tolerance value (scalar) for dat; threshold to within which
            dat must return to be counted as good, after exceeding a ddatdx
            threshold detected bad data.

    References:

        OOI (2012). Data Product Specification for Gradient Test. Document
            Control Number 1341-100010.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10010_Data_Product_SPEC_GRDTEST_OOI.pdf)
    """
    # Sanity checks on dat and x
    dat = np.atleast_1d(dat)
    x = np.atleast_1d(x)

    if strict_validation:
        if not utils.isvector(dat) or not utils.isvector(x):
            raise ValueError('\'dat\' and \'x\' must be vectors')

        if len(dat) != len(x):
            raise ValueError('\'dat\' and \'x\' must be of equal len')

        if not all(np.diff(x) > 0):
            raise ValueError('\'x\' must be montonically increasing')

    # remove any NaNs from input vectors
    dat = dat[~np.isnan(dat)]
    x = x[~np.isnan(x)]

    # flatten input vectors
    dat = dat.flatten()
    x = x.flatten()

    # Check & set mindx
    if utils.isempty(mindx):
        mindx = np.nan

    if np.isnan(mindx):
        mindx = 0

    if strict_validation:
        if not utils.isscalar(mindx):
            raise ValueError('\'mindx\' must be scalar, NaN, or empty.')

    # Apply mindx
    dx = np.diff(x) > mindx
    ff = dx.nonzero()[0]
    gg = np.hstack((np.zeros(1), ff+1)).astype('int8')
    dat = dat[gg]
    x = x[gg]

    # Confirm that there are still data points left, else abort:
    outqc = np.zeros(len(dat), dtype='int8')
    ll = len(dat)
    if ll <= 1:
        log.warn('\'dat\' and \'x\' contain too few points for '
                 'meaningful analysis.')
        outqc = outqc.astype('int8')
        outdat = dat
        outx = x
        return outdat, outx, outqc

    # Check & set startdat, including output for data point 1
    if utils.isempty(startdat):
        startdat = np.nan

    if np.isnan(startdat):
        startdat = dat[0]
        outqc[0] = 1
    else:
        if np.abs(startdat - dat[0]) <= toldat:
            startdat = dat[0]
            outqc[0] = 1
        else:
            outqc[0] = 0

    if strict_validation:
        if not utils.isscalar(startdat):
            raise ValueError('\'startdat\' must be scalar, NaN, or empty.')

    # Main loop, checking for data points 2 through ll
    ii = 1
    while ii < ll:
        if outqc[ii-1] == 0:
            if np.abs(dat[ii] - startdat) <= toldat:
                outqc[ii] = 1
                startdat = dat[ii]
            else:
                outqc[ii] = 0
        else:
            tmp = (dat[ii] - dat[ii-1]) / (x[ii] - x[ii-1])
            if (tmp < ddatdx[0]) or (tmp > ddatdx[1]):
                outqc[ii] = 0
            else:
                outqc[ii] = 1
                startdat = dat[ii]
        ii = ii+1

    outqc = outqc.astype('int8')
    outdat = dat
    outx = x

    return outdat, outx, outqc


def dataqc_solarelevation(lon, lat, dt):
    """
    Description

        Computes instantaneous no-sky solar radiation and altitude from date
        and time stamp and position data. It is put together from expressions
        taken from Appendix E in the 1978 edition of Almanac for Computers,
        Nautical Almanac Office, U.S. Naval Observatory. They are reduced
        accuracy expressions valid for the years 1800-2100. Solar declination
        computed from these expressions is accurate to at least 1'. The solar
        constant (1368.0 W/m^2) represents a mean of satellite measurements
        made over the last sunspot cycle (1979-1995) taken from Coffey et al
        (1995), Earth System Monitor, 6, 6-10.

        This code is a python implementation of soradna1.m available in Air-Sea
        Toolbox.

    Implemented by:

        1997-03-08: Version 1.0 (author unknown) of soradna1.m.
        1998-08-28: Version 1.1 (author unknown) of soradna1.m.
        1999-08-05: Version 2.0 (author unknown) of soradna1.m.

        2013-04-07: Christopher Wingard. Initial python implementation. Note,
        this function is derived from old, unmaintained code. More robust
        implementations exist (e.g. PyEphem and PySolar) that will probably
        calculate these values more accurately.

    Usage:

        z, sorad = dataqc_solarelevation(lon, lat, dt)

            where

        z = solar altitude [degrees]
        sorad = no atmosphere solar radiation [W m^-2]

        lon = longitude (east is positive) [decimal degress]
        lat = latitude [decimal degrees]
        dt = date and time stamp in UTC [seconds since 1970-01-01]

    Examples

        dt = 1329177600     # 2012-02-14 00:00:00
        z, sorad = dataqc_solarelevation(120, 30, dt)
        z = 15.1566, sorad = 366.8129

        OOI (2012). Data Product Specification for Solar Elevation. Document
            Control Number 1341-100011.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10011_Data_Product_SPEC_SOLRELV_OOI.pdf)
    """
    # Test lengths and types of inputs. Latitude and longitude must be the same
    # size and can either be a scalar or a vecotr. The date and time stamp
    # can also be either a scalar or a vector. If all three inputs are vectors,
    # they must be of the same length.
    if len(lon) != len(lat):
        raise ValueError('\'lon\' and \'lat\' must be the same size')

    if utils.isvector(lon) and utils.isvector(lat) and utils.isvector(dt):
        # test their lengths
        if not len(lon) == len(lat) == len(dt):
            raise ValueError('If all inputs are vectors, these must all '
                             'be of the same length')

    # set constants (using values from as_consts.m)
    # ------ short-wave flux calculations
    #   the solar constant [W m^-2] represents a mean of satellite measurements
    #   made over the last sunspot cycle (1979-1995), taken from Coffey et al.
    #   (1995), Earth System Monitor, 6, 6-10.
    solar_const = 1368.0

    # Create a time tuple in UTC from the Epoch time input, and then create
    # scalars or numpy arrays of time elements for subsequent calculations.
    ldt = len(dt)
    yy = np.zeros(ldt, dtype=np.int)
    mn = np.zeros(ldt, dtype=np.int)
    dd = np.zeros(ldt, dtype=np.int)
    hh = np.zeros(ldt, dtype=np.int)
    mm = np.zeros(ldt, dtype=np.int)
    ss = np.zeros(ldt, dtype=np.int)
    for i in range(ldt):
        # create time tuple in UTC
        gtime = time.gmtime(dt[i])
        # create scalar elements
        yy[i] = gtime[0]
        mn[i] = gtime[1]
        dd[i] = gtime[2]
        hh[i] = gtime[3]
        mm[i] = gtime[4]
        ss[i] = gtime[5]

    #constants used in function
    deg2rad = np.pi / 180.0
    rad2deg = 1 / deg2rad

    # compute Universal Time in hours
    utime = hh + (mm + ss / 60.0) / 60.0

    # compute Julian ephemeris date in days (Day 1 is 1 Jan 4713 B.C. which
    # equals -4712 Jan 1)
    jed = (367.0 * yy - np.fix(7.0*(yy+np.fix((mn+9)/12.0))/4.0)
           + np.fix(275.0*mn/9.0) + dd + 1721013 + utime / 24.0)

    # compute interval in Julian centuries since 1900
    jc_int = (jed - 2415020.0) / 36525.0

    # compute mean anomaly of the sun
    ma_sun = 358.475833 + 35999.049750 * jc_int - 0.000150 * jc_int**2
    ma_sun = (ma_sun - np.fix(ma_sun/360.0) * 360.0) * deg2rad

    # compute mean longitude of sun
    ml_sun = 279.696678 + 36000.768920 * jc_int + 0.000303 * jc_int**2
    ml_sun = (ml_sun - np.fix(ml_sun/360.0) * 360.0) * deg2rad

    # compute mean anomaly of Jupiter
    ma_jup = 225.444651 + 2880.0 * jc_int + 154.906654 * jc_int
    ma_jup = (ma_jup - np.fix(ma_jup/360.0) * 360.0) * deg2rad

    # compute longitude of the ascending node of the moon's orbit
    an_moon = (259.183275 - 1800 * jc_int - 134.142008 * jc_int
               + 0.002078 * jc_int**2)
    an_moon = (an_moon - np.fix(an_moon/360.0) * 360.0 + 360.0) * deg2rad

    # compute mean anomaly of Venus
    ma_ven = (212.603219 + 58320 * jc_int + 197.803875 * jc_int
              + 0.001286 * jc_int**2)
    ma_ven = (ma_ven - np.fix(ma_ven/360.0) * 360.0) * deg2rad

    # compute sun theta
    theta = (0.397930 * np.sin(ml_sun) + 0.009999 * np.sin(ma_sun-ml_sun)
             + 0.003334 * np.sin(ma_sun+ml_sun) - 0.000208 * jc_int
             * np.sin(ml_sun) + 0.000042 * np.sin(2*ma_sun+ml_sun) - 0.000040
             * np.cos(ml_sun) - 0.000039 * np.sin(an_moon-ml_sun) - 0.000030
             * jc_int * np.sin(ma_sun-ml_sun) - 0.000014
             * np.sin(2*ma_sun-ml_sun) - 0.000010
             * np.cos(ma_sun-ml_sun-ma_jup) - 0.000010 * jc_int
             * np.sin(ma_sun+ml_sun))

    # compute sun rho
    rho = (1.000421 - 0.033503 * np.cos(ma_sun) - 0.000140 * np.cos(2*ma_sun)
           + 0.000084 * jc_int * np.cos(ma_sun) - 0.000033
           * np.sin(ma_sun-ma_jup) + 0.000027 * np.sin(2.*ma_sun-2.*ma_ven))

    # compute declination
    decln = np.arcsin(theta/np.sqrt(rho))

    # compute equation of time (in seconds of time)
    l = 276.697 + 0.98564734 * (jed-2415020.0)
    l = (l - 360.0 * np.fix(l/360.0)) * deg2rad
    eqt = (-97.8 * np.sin(l) - 431.3 * np.cos(l) + 596.6 * np.sin(2*l)
           - 1.9 * np.cos(2*l) + 4.0 * np.sin(3*l) + 19.3 * np.cos(3*l)
           - 12.7 * np.sin(4*l))
    eqt = eqt / 60.0

    # compute local hour angle from global hour angle
    gha = 15.0 * (utime-12) + 15.0 * eqt / 60.0
    lha = gha - lon

    # compute radius vector
    rv = np.sqrt(rho)

    # compute solar altitude
    sz = (np.sin(deg2rad*lat) * np.sin(decln) + np.cos(deg2rad*lat)
          * np.cos(decln) * np.cos(deg2rad*lha))
    z = rad2deg * np.arcsin(sz)

    # compute solar radiation outside atmosphere (defaults to 0 when solar
    # altitude is below the horizon)
    sorad = (solar_const / rv**2) * np.sin(deg2rad * z)
    sorad[z < 0] = 0

    return (z, sorad)


def dataqc_propogateflags(inflags, strict_validation=False):
    """
    Description:

        Propagate "bad" qc flags (from an arbitrary number of source datasets)
        to another (derived) dataset.

        Consider data from an oceanographic CTD (conductivity, temperature, and
        pressure) instrument. From these three time series, you want to compute
        salinity. If any of the three source data (conductivity, temperature,
        pressure) is of bad quality, the salinity will be bad as well. You can
        feed your QC assessment of the former three into this routine, which
        will then give you the combined assessment for the derived (here:
        salinity) property.

    Implemented by:

        2012-07-17: DPS authored by Mathias Lankhorst. Example code provided
        for Matlab.

        2013-04-06: Christopher Wingard. Initial python implementation.

    Usage:

        outflag = dataqc_propogateflags(inflags)

            where

        outflag = a 1-by-N boolean vector that contains 1 where all of the
            inflags are 1, and 0 otherwise.

        inflags = an M-by-N boolean matrix, where each of the M rows contains
            flags of an independent data set such that "0" means bad data and
            "1" means good data.

    References:

        OOI (2012). Data Product Specification for Combined QC Flags. Document
            Control Number 1341-100012.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10012_Data_Product_SPEC_CMBNFLG_OOI.pdf)
    """
    import numpy as np
    from ion_functions import utils

    if strict_validation:
        if not utils.islogical(inflags):
            raise ValueError('\'inflags\' must be \'0\' or \'1\' '
                             'integer flag array')

    array_size = inflags.shape
    nrows = array_size[0]
    if nrows < 2:
        error('\'inflags\' must be at least a two-dimensional array')

    outflag = np.all(inflags, 0)
    return outflag.astype('int8')


def dataqc_condcompress(p_orig, p_new, c_orig, cpcor=-9.57e-8):
    """
    Description:

        Implementation of the Sea-Bird conductivity compressibility correction,
        scaling the input conductivity based on ratio of the original pressure
        and the updated pressure.

    Implemented by:

        2013-04-07: Christopher Wingard. Initial python implementation.

    Usage:

        c_new = dataqc_condcompress(p_orig, p_new, c_orig, cpcor)

            where

        c_new = updated conductivity record [S/m]

        p_orig = original pressure used to calculate original conductivity,
            this typically the L1a PRESWAT [dbar]
        p_new = updated pressure, typically L1b PRESWAT [dbar]
        c_orig = original conductivty record, typically L1a CONDWAT [S/m]
        cpcor = pressure correction coefficient used to calculate original
            conductivity, default is -9.57e-8

    References:

        OOI (2012). Data Product Specification for Conductivity Compressibility
            Correction. Document Control Number 1341-10030.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10030_Data_Product_SPEC_CNDCMPR_OOI.pdf)
    """
    c_new = c_orig * (1 + cpcor * p_orig) / (1 + cpcor * p_new)
    return c_new
