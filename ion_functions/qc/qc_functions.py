#!/usr/bin/env python

"""
@package ion_functions.qc_functions
@file ion_functions/qc_functions.py
@author Christopher Mueller
@brief Module containing QC functions ported from matlab samples in DPS documents
"""

## DO NOT IMPORT AT THIS LEVEL - Perform imports within each function
def dataqc_globalrangetest_minmax(dat, dat_min, dat_max):
    '''
    Python wrapper for dataqc_globalrangetest
    Combines the min/max arguments into list for dataqc_globalrangetest
    '''
    return dataqc_globalrangetest_minmax(dat, [dat_min,dat_max])

def dataqc_globalrangetest(dat, datlim):
    """
    Global Range Quality Control Algorithm as defined in the DPS for SPEC_GLBLRNG - DCN 1341-10004
    https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/466c4915-c777-429a-8946-c90a8f0945b0/1341-10004_Data_Product_SPEC_GLBLRNG_OOI.pdf

    % DATAQC_GLOBALRANGETEST   Data quality control algorithm testing
    %      if measurements fall into a user-defined valid range.
    %      Returns 1 for presumably good data and 0 for data presumed bad.
    %
    % Time-stamp: <2010-07-28 15:16:00 mlankhorst>
    %
    % USAGE:   out=dataqc_globalrangetest(dat,validrange);
    %
    %          out: Boolean, 0 if value is outside range, else 1.
    %          dat: Input dataset, any scalar, vector, or matrix.
    %               Must be numeric and real.
    %          validrange: Two-element vector with the minimum and
    %               maximum values considered to be valid
    %
    % EXAMPLE:
    %
    %     >> x=[17 16 17 18 25 19];
    %     >> qc=dataqc_globalrangetest(x,[10 20])
    %
    %     qc =
    %
    %          1     1     1     1     0     1
    %
    %
    function out=dataqc_globalrangetest(dat,datlim);

        if ~isnumeric(dat)
            error('DAT must be numeric.')
          end
        if ~all(isreal(dat(:)))
            error('DAT must be real.')
          end
        if ~isnumeric(datlim)
            error('VALIDRANGE must be numeric.')
          end
        if ~all(isreal(datlim(:)))
            error('VALIDRANGE must be real.')
          end
        if length(datlim)~=2
            error('VALIDRANGE must be two-element vector.')
          end
        datlim=[min(datlim(:)) max(datlim(:))];
        out=(dat>=datlim(1))&(dat<=datlim(2))

    """
    import numpy as np
    from ion_functions import utils

    dat_arr = np.atleast_1d(dat)
    datlim_arr = np.atleast_1d(datlim)

    if not utils.isnumeric(dat_arr).all():
        raise ValueError('\'dat\' must be numeric')

    if not utils.isreal(dat_arr).all():
        raise ValueError('\'dat\' must be real')

    if not utils.isnumeric(datlim_arr).all():
        raise ValueError('\'datlim\' must be numeric')

    if not utils.isreal(datlim_arr).all():
        raise ValueError('\'datlim\' must be real')

    if len(datlim_arr) < 2:  # Must have at least 2 elements
        raise ValueError('\'datlim\' must have at least 2 elements')

    return (datlim_arr.min() <= dat) & (dat <= datlim_arr.max()).astype('int8')


def dataqc_spiketest(dat, acc, N=5, L=5):
    """
    Spike Test Quality Control Algorithm as defined in the DPS for SPEC_SPKETST - DCN 1341-10006
    https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/eadad62c-ec80-403d-b3d3-c32c79f9e9e4/1341-10006_Data_Product_SPEC_SPKETST_OOI.pdf

    % DATAQC_SPIKETEST   Data quality control algorithm testing a time
    %                    series for spikes. Returns 1 for presumably
    %                    good data and 0 for data presumed bad.
    %
    % Time-stamp: <2010-07-28 14:25:42 mlankhorst>
    %
    % METHODOLOGY: The time series is divided into windows of length L
    %   (an odd integer number). Then, window by window, each value is
    %   compared to its (L-1) neighboring values: a range R of these
    %   (L-1) values is computed (max. minus min.), and replaced with
    %   the measurement accuracy ACC if ACC>R. A value is presumed to
    %   be good, i.e. no spike, if it deviates from the mean of the
    %   (L-1) peers by less than a multiple of the range, N*max(R,ACC).
    %
    %   Further than (L-1)/2 values from the start or end points, the
    %   peer values are symmetrically before and after the test
    %   value. Within that range of the start and end, the peers are
    %   the first/last L values (without the test value itself).
    %
    %   The purpose of ACC is to restrict spike detection to deviations
    %   exceeding a minimum threshold value (N*ACC) even if the data
    %   have little variability. Use ACC=0 to disable this behavior.
    %
    %
    % USAGE:   out=dataqc_spiketest(dat,acc,N,L);
    %    OR:   out=dataqc_spiketest(dat,acc);
    %
    %          out: Boolean. 0 for detected spike, else 1.
    %          dat: Input dataset, a real numeric vector.
    %          acc: Accuracy of any input measurement.
    %          N (optional, defaults to 5): Range multiplier, cf. above
    %          L (optional, defaults to 5): Window length, cf. above
    %
    % EXAMPLE:
    %
    %    >> x=[-4     3    40    -1     1    -6    -6     1];
    %    >> dataqc_spiketest(x,.1)
    %
    %    ans =
    %
    %         1     1     0     1     1     1     1     1
    %
    function out=dataqc_spiketest(varargin);

    error(nargchk(2,4,nargin,'struct'))
    dat=varargin{1};
    acc=varargin{2};
    N=5;
    L=5;
    switch nargin
        case 3,
            if ~isempty(varargin{3})Data Product Specification for Spike Test
                Ver 1-01 1341-10006 Appendix Page A-2
                N=varargin{3};
            end
        case 4,
            if ~isempty(varargin{3})
                N=varargin{3};
            end
            if ~isempty(varargin{4})
                L=varargin{4};
            end
    end
    if ~isnumeric(dat)
        error('DAT must be numeric.')
    end
    if ~isvector(dat)
        error('DAT must be a vector.')
    end
    if ~isreal(dat)
        error('DAT must be real.')
    end
    if ~isnumeric(acc)
        error('ACC must be numeric.')
    end
    if ~isscalar(acc)
        error('ACC must be scalar.')
    end
    if ~isreal(acc)
        error('ACC must be real.')
    end
    if ~isnumeric(N)
        error('N must be numeric.')
    end
    if ~isscalar(N)
        error('N must be scalar.')
    end
    if ~isreal(N)
        error('N must be real.')
    end
    if ~isnumeric(L)
        error('L must be numeric.')
    end
    if ~isscalar(L)
        error('L must be scalar.')
    end
    if ~isreal(L)
        error('L must be real.')
    end
    L=ceil(abs(L));
    if (L/2)==round(L/2)
        L=L+1;
        warning('L was even; setting L:=L+1')
    end
    if L<3
        L=5;
        warning('L was too small; setting L:=5')
    end
    ll=length(dat);

    L2=(L-1)/2;
    i1=1+L2;
    i2=ll-L2;

    if ll>=L

        for ii=i1:i2
            tmpdat=dat(ii+[-L2:-1 1:L2]);
            R=max(tmpdat)-min(tmpdat);
            R=max([R acc]);
            if (N*R)>abs(dat(ii)-mean(tmpdat))
                out(ii)=1;
            end
        end
        for ii=1:L2
            tmpdat=dat([1:ii-1 ii+1:L]);
            R=max(tmpdat)-min(tmpdat);
            R=max([R acc]);
            if (N*R)>abs(dat(ii)-mean(tmpdat))
                out(ii)=1;
            end
        end
        for ii=ll-L2+1:ll
            tmpdat=dat([ll-L+1:ii-1 ii+1:ll]);
            R=max(tmpdat)-min(tmpdat);
            R=max([R acc]);
            if (N*R)>abs(dat(ii)-mean(tmpdat))
                out(ii)=1;
            end
        end
    else
        warning('L was greater than length of DAT, returning zeros.')
    end

    """

    import numpy as np
    from ion_functions import utils

    dat_arr = np.atleast_1d(dat)

    if not utils.isnumeric(dat_arr).all():
        raise ValueError('\'dat\' must be numeric')

    if not utils.isvector(dat_arr):
        raise ValueError('\'dat\' must be a vector')

    if not utils.isreal(dat_arr).all():
        raise ValueError('\'dat\' must be real')

    for k, arg in {'acc': acc, 'N': N, 'L': L}.iteritems():
        if not utils.isnumeric(arg).all():
            raise ValueError('\'{0}\' must be numeric'.format(k))

        if not utils.isscalar(arg):
            raise ValueError('\'{0}\' must be a scalar'.format(k))

        if not utils.isreal(arg).all():
            raise ValueError('\'{0}\' must be real'.format(k))

    L = np.ceil(np.abs(L))
    if L / 2 == np.round(L / 2):
        L += 1
        # Warn - L was even; setting L = L + 1
    if L < 3:
        L = 5
        # Warn - L was too small; setting L = 5

    ll = len(dat_arr)
    out = np.zeros(dat_arr.size, dtype='int8')

    L2 = int((L - 1) / 2)
    i1 = 1 + L2
    i2 = ll - L2

    if ll >= L:

        for ii in xrange(i1 - 1,i2):  # for ii=i1:i2
            tmpdat = np.hstack((dat_arr[ii - L2:ii], dat_arr[ii + 1:ii + 1 + L2]))  # tmpdat=dat(ii+[-L2:-1 1:L2]);
            R = tmpdat.max() - tmpdat.min()  # R=max(tmpdat)-min(tmpdat);
            R = np.max([R, acc])  # R=max([R acc]);
            if (N * R) > np.abs(dat_arr[ii] - tmpdat.mean()):  # if (N*R)>abs(dat(ii)-mean(tmpdat))
                out[ii] = 1  # out(ii)=1;

        for ii in xrange(L2):  # for ii=1:L2
            tmpdat = np.hstack((dat_arr[:ii], dat_arr[ii+1:L]))  # tmpdat=dat([1:ii-1 ii+1:L]);
            R = tmpdat.max() - tmpdat.min()  # R=max(tmpdat)-min(tmpdat);
            R = np.max([R, acc])  # R=max([R acc]);
            if (N * R) > np.abs(dat_arr[ii] - tmpdat.mean()):  # if (N*R)>abs(dat(ii)-mean(tmpdat))
                out[ii] = 1  # out(ii)=1;

        for ii in xrange(ll - L2, ll):  # for ii=ll-L2+1:ll
            tmpdat = np.hstack((dat_arr[:ii], dat_arr[ii:L]))  # tmpdat=dat([ll-L+1:ii-1 ii+1:ll]);
            R = tmpdat.max() - tmpdat.min()  # R=max(tmpdat)-min(tmpdat);
            R = np.max([R, acc])  # R=max([R acc]);
            if (N * R) > np.abs(dat_arr[ii] - tmpdat.mean()):  # if (N*R)>abs(dat(ii)-mean(tmpdat))
                out[ii] = 1  # out(ii)=1;

    else:
        pass
        # Warn - 'L was greater than length of DAT, returning zeros.'

    return out


def dataqc_stuckvaluetest(x, reso, num=10):
    """
    Stuck Value Test Quality Control Algorithm as defined in the DPS for SPEC_STUCKVL - DCN 1341-10008
    https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/a04acb56-7e27-48c6-a40b-9bb9374ee35c/1341-10008_Data_Product_SPEC_STUCKVL_OOI.pdf

    % DATAQC_STUCKVALUETEST   Data quality control algorithm testing a
    %     time series for "stuck values", i.e. repeated occurences of
    %     one value. Returns 1 for presumably good data and 0 for data
    %     presumed bad.
    %
    % Time-stamp: <2011-10-31 11:20:23 mlankhorst>
    %
    % USAGE:   OUT=dataqc_stuckvaluetest(X,RESO,NUM);
    %
    %       OUT:  Boolean output: 0 where stuck values are found,
    %             1 elsewhere.
    %       X:    Input time series (vector, numeric).
    %       RESO: Resolution; repeat values less than RESO apart will
    %             be considered "stuck values".
    %       NUM:  Minimum number of successive values within RESO of
    %             each other that will trigger the "stuck value". NUM
    %             is optional and defaults to 10 if omitted or empty.
    %
    % EXAMPLE:
    %
    % >> x=[4.83  1.40  3.33  3.33  3.33  3.33  4.09  2.97  2.85  3.67];
    %
    % >> dataqc_stuckvaluetest(x,.001,4)
    %
    % ans =
    %
    %       1     1     0     0     0     0     1     1     1     1
    %
    function out=dataqc_stuckvaluetest(varargin);

    error(nargchk(2,3,nargin,'struct'))
    x=varargin{1};
    reso=varargin{2};
    num=10;
    switch nargin
        case 3,
            if ~isempty(varargin{3})
                num=varargin{3};
            end
    end
    if ~isnumeric(x)
        error('X must be numeric.')
    end
    if ~isvector(x)
        error('X must be a vector.')
    end
    if ~isnumeric(reso)
        error('RESO must be numeric.')
    end
    if ~isscalar(reso)
        error('RESO must be a scalar.')
    end
    if ~isreal(reso)
        error('RESO must be real.')
    end
    reso=abs(reso);
    if ~isnumeric(num)
        error('NUM must be numeric.')
    end
    if ~isscalar(num)
        error('NUM must be a scalar.')
    end
    if ~isreal(num)
        error('NUM must be real.')
    end
    num=abs(num);
    ll=length(x);
    out=zeros(size(x));
    out=logical(out);
    if ll<num
        warning('NUM is greater than length(X). Returning zeros.')
    else
        out=ones(size(x));
        iimax=ll-num+1;
        for ii=1:iimax
            ind=[ii:ii+num-1];
            tmp=abs(x(ii)-x(ind));
            if all(tmp<reso)
                out(ind)=0;
            end
        end
    end
    out=logical(out);
    """

    import numpy as np
    from ion_functions import utils

    dat_arr = np.atleast_1d(x)

    if not utils.isnumeric(dat_arr).all():
        raise ValueError('\'x\' must be numeric')

    if not utils.isvector(dat_arr):
        raise ValueError('\'x\' must be a vector')

    if not utils.isreal(dat_arr).all():
        raise ValueError('\'x\' must be real')

    for k, arg in {'reso': reso, 'num': num}.iteritems():
        if not utils.isnumeric(arg).all():
            raise ValueError('\'{0}\' must be numeric'.format(k))

        if not utils.isscalar(arg):
            raise ValueError('\'{0}\' must be a scalar'.format(k))

        if not utils.isreal(arg).all():
            raise ValueError('\'{0}\' must be real'.format(k))

    num = np.abs(num)
    ll = len(x)
    out = np.zeros(dat_arr.size, dtype='int8')

    if ll < num:
        # Warn - 'num' is greater than length(x), returning zeros
        pass
    else:
        out.fill(1)
        iimax = ll - num+1
        for ii in xrange(iimax):
            slice_ = slice(ii, ii + num)
            tmp = np.abs(dat_arr[ii] - dat_arr[slice_])
            if (tmp < reso).all():
                out[slice_] = 0

    return out


def dataqc_polytrendtest(dat, t, ord_n=1, nstd=3):
    """
    Stuck Value Test Quality Control Algorithm as defined in the DPS for SPEC_TRNDTST - DCN 1341-10007
    https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/c33037ab-9dd5-4615-8218-0957f60a47f3/1341-10007_Data_Product_SPEC_TRNDTST_OOI.pdf

    % DATAQC_POLYTRENDTEST Data quality control algorithm testing
    % if measurements contain a significant portion of a polynomial.
    % Returns 1 if this is not the case, else 0.
    %
    % Time-stamp: <2010-10-29 13:56:46 mlankhorst>
    %
    % RATIONALE: The purpose of this test is to check if a significant
    % fraction of the variability in a time series can be explained
    % by a drift, possibly interpreted as a sensor drift. This drift
    % is assumed to be a polynomial of order ORD. Use ORD=1 to
    % consider a linear drift
    %
    % METHODOLOGY: The time series DAT is passed to MatLab's POLYFIT
    % routine to obtain a polynomial fit PP to DAT, and the
    % difference DAT-PP is compared to the original DAT. If the
    % standard deviation of (DAT-PP) is less than that of DAT by a
    % factor of NSTD, the time series is assumed to contain a
    % significant trend (output will be 0), else not (output will be
    % 1).
    %
    % USAGE: OUT=dataqc_polytrendtest(DAT,ORD,NSTD);
    %
    % OUT: Boolean scalar, 0 if trend is detected, 1 if not.
    %
    % DAT: Input dataset, a numeric real vector.
    % ORD (optional, defaults to 1): Polynomial order.
    % NSTD (optional, defaults to 3): Factor by how much the
    % standard deviation must be reduced before OUT
    % switches from 1 to 0
    %
    function out=dataqc_polytrendtest(varargin);
    error(nargchk(1,3,nargin,'struct'))
    dat=varargin{1};
    if ~isnumeric(dat)
        error('DAT must be numeric.')
    end
    if ~isvector(dat)
        error('DAT must be vector.')
    end
    if ~isreal(dat)
        error('DAT must be real.')
    end
    ord=1;
    nstd=3;
    if nargin==2
        if ~isempty(varargin{2})
            ord=varargin{2};
        end
    end
    if nargin==3
        if ~isempty(varargin{2})
            ord=varargin{2};
        end
        if ~isempty(varargin{3})
            nstd=varargin{3};
        end
    end
    if ~isnumeric(ord)
        error('ORD must be numeric.')
    end
    if ~isscalar(ord)
        error('ORD must be scalar.')
    end
    if ~isreal(ord)
        error('ORD must be real.')
    end
    if ~isnumeric(nstd)
        error('NSTD must be numeric.')
    end
    if ~isscalar(nstd)
        error('NSTD must be scalar.')
    end
    if ~isreal(nstd)
        error('NSTD must be real.')
    end
    ord=round(abs(ord));
    nstd=abs(nstd);
    ll=length(dat);
    x=[1:ll];
    pp=polyfit(x,dat,ord);
    datpp=polyval(pp,x);
    if (nstd*std(dat-datpp))<std(dat)
        out=0;
    else
        out=1;
    end
    """

    import numpy as np
    from ion_functions import utils

    dat_arr = np.atleast_1d(dat)

    if not utils.isnumeric(dat_arr).all():
        raise ValueError('\'dat\' must be numeric')

    if not utils.isvector(dat_arr):
        raise ValueError('\'dat\' must be a vector')

    if not utils.isreal(dat_arr).all():
        raise ValueError('\'dat\' must be real')

    for k, arg in {'ord_n': ord_n, 'nstd': nstd}.iteritems():
        if not utils.isnumeric(arg).all():
            raise ValueError('\'{0}\' must be numeric'.format(k))

        if not utils.isscalar(arg):
            raise ValueError('\'{0}\' must be a scalar'.format(k))

        if not utils.isreal(arg).all():
            raise ValueError('\'{0}\' must be real'.format(k))

    ord_n = int(round(abs(ord_n)))
    nstd = int(abs(nstd))
    # Not needed because time is incorporated as 't'
    # ll = len(dat)
    # t = range(ll)
    pp = np.polyfit(t, dat, ord_n)
    datpp = np.polyval(pp, t)

    if np.atleast_1d((np.std(dat - datpp) * nstd) < np.std(dat)).all():
        return 0

    return 1



