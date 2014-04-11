#!/usr/bin/env python
"""
@package ion_functions.data.wav_functions
@file ion_functions/data/wav_functions.py
@author Russell Desiderio
@brief Module containing WAVSS wave statistics data-calculations.
"""
import numpy as np

from ion_functions.data.generic_functions import magnetic_declination
from ion_functions.utils import fill_value


def wav_triaxys_dir_freq(nfreq_nondir, nfreq_dir, freq0, delta_freq):
    """
    FLAG:

        The variable nfreq_dir put out by the WAVSS instrument and therefore
        also the data product WAVSTAT-FDS_L1 can vary within each datapacket
        based solely on measured ocean conditions (with unchanged instrument
        settings). The numbers of values in the L0 data products WAVSTAT_PDS
        and WAVSTAT_SDS for each datapacket are also determined by the nfreq_dir
        values, as is true with the WAVSTAT-DDS_L2 data product (see function
        def wav_triaxys_correct_directional_wave_direction).

    Description:

        Function to compute the WAVSTAT-FDS_L1 data product (frequency values for
        directional wave spectral bins) for the WAVSS instrument class (TRIAXYS
        Wave Sensor, manufactured by AXYS Technologies).

    Implemented by:

        2014-04-03: Russell Desiderio.  Initial code.

    Usage:

        fds = wav_triaxys_dir_freq(nfreq_nondir, nfreq_dir, freq0, delta_freq)

            where

        fds =  frequency values for directional wave spectral bins (WAVSTAT-FDS_L1) [Hz]
        nfreq_nondir = number of non-directional wave frequency bins from the value specified in the
            WAVSS $TSPNA (not $TSPMA) data sentence.
        nfreq_dir = number of directional wave frequency bins from the value specified in the WAVSS
            $TSPMA data sentence.
        freq0 = initial frequency value from the value specified in the WAVSS $TSPMA data sentence.
        delta_freq = frequency spacing from the value specified in the WAVSS $TSPMA data sentence.

    References:

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

    """
    # condition input variables.
    # all delta_freq and freq0 values will be floats.
    nfreq_nondir = np.array(nfreq_nondir, ndmin=1)
    nfreq_dir = np.array(nfreq_dir, ndmin=1)
    freq0 = np.array(freq0, ndmin=1)
    delta_freq = np.array(delta_freq, ndmin=1)

    # each data packet may call for a different number of directional frequency values nfreq_dir.
    # however, this number will always be <= nfreq_nondir, and all the nfreq_nondir values will be identical.
    npackets = nfreq_nondir.shape[0]
    fds = np.zeros((npackets, nfreq_nondir[0])) + fill_value

    # for the linspace calculation, which is slightly slower than using arange.
    #freq_end = freq0 + (nfreq_dir - 1) * delta_freq

    for ii in range(npackets):

        #fds[ii, 0:nfreq_dir[ii]] = np.linspace(freq0[ii], freq_end[ii], num=nfreq_dir[ii])
        fds[ii, 0:nfreq_dir[ii]] = freq0[ii] + np.arange(nfreq_dir[ii]) * delta_freq[ii]

    ## return a "rank 1 vector" if fds is a 2D row vector
    #if fds.shape[0] == 1:
    #    fds = np.reshape(fds, (fds.shape[1],))

    return fds


def wav_triaxys_nondir_freq(nfreq, freq0, delta_freq):
    """
    Description:

        Function to compute the WAVSTAT-FND_L1 data product (frequency values for
        non-directional wave spectral bins) for the WAVSS instrument class (TRIAXYS
        Wave Sensor, manufactured by AXYS Technologies).

    Implemented by:

        2014-04-03: Russell Desiderio.  Initial code.

    Usage:

        fnd = wav_triaxys_nondir_freq(nfreq, freq0, delta_freq)

            where

        fnd =  frequency values for non-directional wave spectral bins (WAVSTAT-FND_L1) [Hz]
        nfreq = number of frequency bins from the value specified in the WAVSS $TSPNA data sentence.
        freq0 = initial frequency value from the value specified in the WAVSS $TSPNA data sentence.
        delta_freq = frequency spacing from the value specified in the WAVSS $TSPNA data sentence.

    References:

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

    """
    # condition input variables
    nfreq = np.array(nfreq, ndmin=1)
    freq0 = np.array(freq0, ndmin=1)
    delta_freq = np.array(delta_freq, ndmin=1)

    # each set of data inputs will call for the same number of frequency values nfreq.
    # therefore can vectorize (without using forloop as had to be done for the
    # directional frequencies case) by setting up all variables as 2D arrays
    # of size(npackets, nfreq).
    npackets = nfreq.shape[0]
    n_freqs = nfreq[0]

    # orient arrays such that lead index indexes each set of data inputs
    freq0_2d = np.tile(freq0, (n_freqs, 1)).transpose()
    delta_freq_2d = np.tile(delta_freq, (n_freqs, 1)).transpose()
    steps_2d = np.tile(np.arange(n_freqs), (npackets, 1))
    fnd = freq0_2d + steps_2d * delta_freq_2d

    ## return a "rank 1 vector" if fnd is a 2D row vector
    #if fnd.shape[0] == 1:
    #    fnd = np.reshape(fnd, (fnd.shape[1],))

    return fnd


def wav_triaxys_buoymotion_time(ntp_timestamp, ntime, time0, delta_time):
    """
    Description:

        Function to compute the WAVSTAT-MOTT_L1 data product (time values associated with
        buoy displacement measurements WAVSTAT-MOT[X,Y,Z]) for the WAVSS instrument class
        (TRIAXYS Wave Sensor, manufactured by AXYS Technologies).

    Implemented by:

        2014-04-07: Russell Desiderio.  Initial code.

    Usage:

        mott = wav_triaxys_buoymotion_time(ntp_timestamp, ntime, time0, delta_time):

            where

        mott = NTP times corresponding to buoy displacement data measurements (WAVSTAT-MOTT_L1)
            [secs since 1900-01-01].
        ntp_timestamp = NTP time stamp corresponding to the date and time specified in
            the $TSPHA data sentence [secs since 1900-01-01].
        ntime = number of time values from the value specified in the WAVSS $TSPHA data sentence.
        time0 = time elapsed between ntp_timestamp and time of first WAVSTAT-MOT[XYZ] data point,
            from the value specified in the WAVSS $TSPHA data sentence ("Initial Time") [sec].
        delta_time = time intervals between subsequent buoydisplacement measurement times,
            from the value specified in the WAVSS $TSPHA data sentence ("Time Spacing") [sec].

    References:

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

    """
    # condition input variables;
    # make sure time interval is not type integer
    ntime = np.array(ntime, ndmin=1)
    time0 = np.array(time0, ndmin=1)
    delta_time = np.array(delta_time, dtype='float', ndmin=1)

    # this algorithm is almost identical to that contained in def wav_wavss_nondir_freq above.
    # these are the dimensions of all the 2D arrays used in the calculation
    npackets = ntime.shape[0]
    n_time_values = ntime[0]

    # orient the lead index to iterate over the data packet number
    ntp0_2d = np.tile(ntp_timestamp + time0, (n_time_values, 1)).transpose()
    delta_time_2d = np.tile(delta_time, (n_time_values, 1)).transpose()
    steps_2d = np.tile(np.arange(n_time_values), (npackets, 1))

    mott = ntp0_2d + steps_2d * delta_time_2d

    ## return a "rank 1 vector" if fnd is a 2D row vector
    #if mott.shape[0] == 1:
    #    mott = np.reshape(mott, (mott.shape[1],))

    return mott


def wav_triaxys_correct_mean_wave_direction(dir_raw, lat, lon, ntp_ts):
    """
    Description:

        Function to compute the WAVSTAT-D_L2 data product (mean wave direction corrected for magnetic
        declination) for the WAVSS instrument class (TRIAXYS Wave Sensor, manufactured by AXYS Technologies).

    Implemented by:

        2014-04-08: Russell Desiderio.  Initial code.

    Usage:

        dir_cor = wav_triaxys_correct_mean_wave_direction(dir_raw, lat, lon, ntp_ts)

            where

        dir_cor =  mean wave direction corrected for magnetic declination (WAVSTAT-D_L2) [deg, [0 360)].
        dir_raw =  uncorrected mean wave direction (WAVSTAT-D_L0) [deg, [0 360)].
        lat = latitude of the instrument [decimal degrees].  North is positive, South negative.
        lon = longitude of the instrument [decimal degrees].  East is positive, West negative.
        ntp_ts = NTP time stamp from a data particle [secs since 1900-01-01].

    References:

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

    """
    # calculate the magnetic declination using the WWM2010 model
    # the WAVSS is a surface wave sensor, so that height above sealevel = 0,
    # which is the default value used in the magnetic_declination calculation.
    theta = magnetic_declination(lat, lon, ntp_ts)

    # directions are [0,360) degrees; and magnetic declinations can be positive or negative
    dir_cor = np.mod(dir_raw + theta + 360, 360)

    # return corrected direction
    return dir_cor


def wav_triaxys_correct_directional_wave_direction(dir_raw, lat, lon, ntp_ts):
    """
    FLAG:

        The numbers of values in the L0 and L2 data product WAVSTAT_DDS for each datapacket
        are determined by the values of the nfreq_dir variable, which can vary as a function
        of measured ocean conditions at fixed instrument setting. See also the FLAG note for
        function def wav_triaxys_dir_freq.

    Description:

        Function to compute the WAVSTAT-DDS_L2 data product (directional wave
        directions corrected for magnetic declination) for the WAVSS instrument
        class (TRIAXYS Wave Sensor, manufactured by AXYS Technologies).

    Implemented by:

        2014-04-09: Russell Desiderio.  Initial code.

    Usage:

        dir_cor = wav_triaxys_correct_directional_wave_direction(dir_raw, lat, lon, ntp_ts)

            where

        dir_cor =  directional waves' directions corrected for magnetic declination
            (WAVSTAT-DDS_L2) [deg, [0 360)].
        dir_raw =  uncorrected directional waves' directions (WAVSTAT-DDS_L0) [deg, [0 360)].
        lat = latitude of the instrument [decimal degrees].  North is positive, South negative.
        lon = longitude of the instrument [decimal degrees].  East is positive, West negative.
        ntp_ts = NTP time stamp from a data particle [secs since 1900-01-01].

    References:

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

    """
    # assume that the dir_raw data product comes in as a 2D numpy array with fill values
    # appropriately placed to account for the cases in which the number of reported
    # directional wave frequency bins differs from data packet to data packet (and is
    # less than the number of reported non-directional frequency bins).
    dir_raw = np.array(dir_raw, ndmin=2)

    # change fill values to Nans, so that subsequent array operations will leave the
    # Nan entries unchanged.
    dir_raw[dir_raw == fill_value] = np.nan

    # calculate the magnetic declination using the WWM2010 model
    # the WAVSS is a surface wave sensor, so that height above sealevel = 0,
    # which is the default value used in the magnetic_declination calculation.
    theta = magnetic_declination(lat, lon, ntp_ts)

    # theta in general will be a vector, so replicate it into a matrix to match the dir_raw dimensions.
    theta = np.tile(theta, (dir_raw.shape[1], 1)).transpose()

    # directions are [0,360) degrees; and magnetic declinations can be positive or negative
    dir_cor = np.mod(dir_raw + theta + 360, 360)

    # replace Nans with fills
    dir_cor[np.isnan(dir_cor)] = fill_value

    # return corrected directions
    return dir_cor


def wav_triaxys_magcor_buoymotion_x(x, y, lat, lon, ntp_timestamp):
    """
    Description:

        Function to compute the WAVSTAT-MOTX_L1 data product (eastward buoy displacement)
        for the WAVSS instrument class (TRIAXYS Wave Sensor, manufactured by AXYS Technologies)
        from the WAVSTAT-MOTX_L0 and WAVSTAT-MOTY_L0 data products. All that is required is to
        correct for magnetic declination (variation).

    Implemented by:

        2014-04-10: Russell Desiderio.  Initial code. Uses magnetic declination values calculated
                                        using the WMM 2010 model. WAVSS is a surface sensor, so
                                        that the depth variable for calculating declination is 0
                                        (default value for the magnetic_declination function).

    Usage:

        motx = wav_triaxys_magcor_buoymotion_x(x, y, lat, lon, ntp_timestamp)

            where

        motx =  East displacement of the buoy on which the WAVSS is mounted, corrected for
                magnetic declination (WAVSTAT-MOTX_L1) [m]
        x = uncorrected eastward displacement (WAVSTAT-MOTX_L0) [m]
        y = uncorrected northward displacement (WAVSTAT-MOTY_L0) [m]
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        ntp_timestamp = NTP time stamp corresponding to the date and time specified in
            the $TSPHA data sentence [secs since 1900-01-01].

            Note as to the values of ntp_timestamp used in the calculation:

            The maximum sampling period for this instrument is 35 minutes, during which time
            the magnetic declination will not change. Therefore, to correct for magnetic
            declinaton only one timestamp is required for each ensemble of (x,y) values acquired
            during any given sampling period. All that is necessary, then, is the ntp_timestamp
            specified above, which is the same input ntp_timestamp parameter used in the function
            wav_triaxys_buoymotion_time; it is not necessary to use the vector timestamps in the
            WAVSTAT-MOTT_L1 data product.

    References:

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

    """
    # force shapes of inputs to arrays
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    ntp_timestamp = np.atleast_1d(ntp_timestamp)

    # calculate the magnetic declination using the WWM2010 model.
    # the WAVSS surface wave sensor is at sealevel, which is the default z value for mag dec.
    theta = magnetic_declination(lat, lon, ntp_timestamp)

    # correct for declination by rotating coordinates.
    # the function magnetic_correction_einsum was written to correct (u,v) velocities, but
    # it also applies to (E,N) coordinates.
    motx, _ = magnetic_correction_einsum(theta, x, y)

    # return corrected Eastward buoy displacement(s)
    return motx


def wav_triaxys_magcor_buoymotion_y(x, y, lat, lon, ntp_timestamp):
    """
    Description:

        Function to compute the WAVSTAT-MOTY_L1 data product (northward buoy displacement)
        for the WAVSS instrument class (TRIAXYS Wave Sensor, manufactured by AXYS Technologies)
        from the WAVSTAT-MOTX_L0 and WAVSTAT-MOTY_L0 data products. All that is required is to
        correct for magnetic declination (variation).

    Implemented by:

        2014-04-10: Russell Desiderio.  Initial code. Uses magnetic declination values calculated
                                        using the WMM 2010 model. WAVSS is a surface sensor, so
                                        that the depth variable for calculating declination is 0
                                        (default value for the magnetic_declination function).

    Usage:

        moty = wav_triaxys_magcor_buoymotion_y(x, y, lat, lon, dt)

            where

        moty =  North displacement of the buoy on which the WAVSS is mounted, corrected for
                magnetic declination (WAVSTAT-MOTY_L1) [m]
        x = uncorrected eastward displacement (WAVSTAT-MOTX_L0) [m]
        y = uncorrected northward displacement (WAVSTAT-MOTY_L0) [m]
        lat = instrument's deployment latitude [decimal degrees]
        lon = instrument's deployment longitude [decimal degrees]
        ntp_timestamp = NTP time stamp corresponding to the date and time specified in
            the $TSPHA data sentence [secs since 1900-01-01].

            Note as to the values of ntp_timestamp used in the calculation:
            See Note in Usage section of wav_triaxys_magcor_buoymotion_x.


    References:

        OOI (2012). Data Product Specification for Wave Statistics. Document Control
            Number 1341-00450. https://alfresco.oceanobservatories.org/ (See:
            Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00450_Data_Product_WAVE_STATISTICS_OOI.pdf)

    """
    # force shapes of inputs to arrays
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    ntp_timestamp = np.atleast_1d(ntp_timestamp)

    # calculate the magnetic declination using the WWM2010 model.
    # the WAVSS surface wave sensor is at sealevel, which is the default z value for mag dec.
    theta = magnetic_declination(lat, lon, ntp_timestamp)

    # correct for declination by rotating coordinates.
    # the function magnetic_correction_einsum was written to correct (u,v) velocities, but
    # it also applies to (E,N) coordinates.
    _, moty = magnetic_correction_einsum(theta, x, y)

    # return corrected Northward buoy displacement(s)
    return moty


def magnetic_correction_einsum(theta, u, v):
    """
    Description:

        The executable code in this function is identical to that in the function
        magnetic_correction_vctrzd in the module adcp_functions. At some point in
        the future the function magnetic_correction in generic_functions may be
        deprecated and replaced with this vectorized and much faster version.

        This function corrects velocity profiles for the magnetic variation
        (declination) at the measurement location.  The magnetic declination
        is obtained from the 2010 World Magnetic Model (WMM2010) provided by
        NOAA (see wmm_declination).

        This version handles 'vectorized' input variables without using for
        loops. It was specifically written to handle the case of a 1D array of
        theta values, theta=f(i), with corresponding sets of 'u' and 'v' values
        such that u=f(i,j) and v=f(i,j), where there are j 'u' and 'v' values
        for each theta(i).

    Implemented by:

        2014-04-04: Russell Desiderio. Initial code. This function is used to
                    calculate magnetic corrections by the functions contained
                    in this module instead of the function magnetic_correction
                    found in ion_functions.data.generic_functions.

    Usage:

        u_cor, v_cor = magnetic_correction(theta, u, v)

            where

        u_cor = eastward velocity profiles, in earth coordinates, with
            the correction for magnetic variation applied.
        v_cor = northward velocity profiles, in earth coordinates,
            with the correction for magnetic variation applied.

        theta = magnetic variation based on location (latitude, longitude and
            altitude) and date; units of theta are [degrees]
        u = uncorrected eastward velocity profiles in earth coordinates
        v = uncorrected northward velocity profiles in earth coordinates

    References:

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        OOI (2013). Data Product Specification for Turbulent Velocity Profile
            and Echo Intensity. Document Control Number 1341-00760.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00760_Data_Product_SPEC_VELPROF_OOI.pdf)
    """
    # force shapes of inputs to arrays
    theta = np.atleast_1d(theta)
    uv = np.atleast_2d(u)
    v = np.atleast_2d(v)

    theta_rad = np.radians(theta)
    cosT = np.cos(theta_rad)
    sinT = np.sin(theta_rad)

    # set up rotation matrix
    M = np.array([[cosT, sinT],
                  [-sinT, cosT]])

    # roll axes so that the lead index represents data packet #.
    M = np.rollaxis(M, 2)

    # construct the uncorrected velocity matrix.
    # the coordinate system is 2D, so the middle dimension is sized at 2.
    uv = np.zeros((u.shape[0], 2, u.shape[1]))

    # load the coordinates to be rotated into the appropriate slices
    uv[:, 0, :] = u
    uv[:, 1, :] = v

    # the Einstein summation is here configured to do the matrix
    # multiplication uv_cor(i,k) = M(i,j) * uv(j,k) on each slice h.
    uv_cor = np.einsum('hij,hjk->hik', M, uv)

    # the magnetically corrected u values are:
    u_cor = uv_cor[:, 0, :]

    # the magnetically corrected v values are:
    v_cor = uv_cor[:, 1, :]

    # return corrected u and v values
    return (u_cor, v_cor)
