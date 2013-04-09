% DATAQC_GRADIENTTEST Data quality control algorithm testing if
% changes between successive data points fall within a certain
% range.
%
% Time-stamp: <2012-04-26 13:47:12 mlankhorst>
%
% Input data DAT are given as a function of coordinate X. The
% algorithm will flag DAT values as bad if the change
% deltaDAT/deltaX between successive DAT values exceeds thresholds
% given in DDATDX. Once the threshold is exceeded, following DAT
% are considered bad until a DAT value returns to within TOLDAT of
% the last known good value.
%
% It is possible to remove data points that are too close together
% in X coordinates (use MINDX).
%
% By default, the first value of DAT is considered good. To change
% this, use STARTDAT and TOLDAT to set as the first good data point
% the first one that comes within TOLDAT of STARTDAT.
%
% USAGE: [OUTDAT,OUTX,OUTQC]= ...
% dataqc_gradienttest(DAT,X,DDATDX,MINDX,STARTDAT,TOLDAT);
%
% DAT: Input dataset, a numeric real vector.
% X: Coordinate (e.g. time, distance) along which DAT is
% given. Must be of the same size as DAT and strictly
% increasing.
% DDATDX: Two-element vector defining the valid range of
% deltaDAT/deltaX from one point to the next.
% MINDX: Scalar. Minimum deltaX for which this test will
% be applied (data that are less than MINDX apart will be
% deleted). Defaults to zero if NaN/empty.
% STARTDAT: Start value (scalar) of DAT that is presumed
% good. Defaults to first non-NaN value of DAT if NaN/empty.
% TOLDAT: Tolerance value (scalar) for DAT; threshold to within
% which DAT must return to be counted as good, after
% exceeding a DDATDX threshold detected bad data.
%
% OUTDAT: Same as DAT except that NaNs and values not meeting
% MINDX are removed.
% OUTX: Same as X except that NaNs and values not meeting
% MINDX are removed.
% OUTQC: Output quality control flags for OUTDAT. 0 means bad
% data, 1 means good data.
%
%
% EXAMPLES:
%
% Ordinary use, default MINDX and STARTDAT:
%
% [outdat,outx,outqc]= ...
% dataqc_gradienttest([3 5 98 99 4],[1:5],[-50 50],[],[],5)
% outdat = 3 5 98 99 4
% outx = 1 2 3 4 5
% outqc = 1 1 0 0 1
%
%
% Alternate STARTDAT to swap good/bad segments:
%
% [outdat,outx,outqc]= ...
% dataqc_gradienttest([3 5 98 99 4],[1:5],[-50 50],[],100,5)
% outdat = 3 5 98 99 4
% outx = 1 2 3 4 5
% outqc = 0 0 1 1 0
%
%
% Alternate MINDX to remove certain X and DAT:
%
% [outdat,outx,outqc]= ...
% dataqc_gradienttest([3 5 98 99 4],[1 2 3 3.1 4], ...
% [-50 50],0.2,[],5)
% outdat = 3 5 98 4
% outx = 1 2 3 4
% outqc = 1 1 0 1
%
function [outdat,outx,outqc] = dataqc_gradienttest(dat,x,ddatdx,mindx,startdat,toldat);

% Sanity checks on DAT and X:
if ((~isvector(dat))|(~isvector(x)))
    error('DAT and X must be vectors.')
end
if (length(dat))~=(length(x))
    error('DAT and X must be of equal length.')
end
if ~all((diff(x))>0)
    error('X must be strictly monotonically increasing.')
end

ff=find((~isnan(dat))&(~isnan(x)));
dat=dat(ff);
x=x(ff);

dat=dat(:)';
x=x(:)';

% Check & set MINDX
if isempty(mindx)
    mindx=nan;
end
if isnan(mindx)
    mindx=0;
end
if ~isscalar(mindx)
    error('MINDX must be scalar, NaN, or empty.')
end

% Apply MINDX
dx=diff(x);
ff=find(dx>mindx);
gg=[1 ff+1];
dat=dat(gg);
x=x(gg);

% Confirm that there are still data points left, else abort:
outqc=zeros(size(dat));
ll=length(dat);
if ll<=1
    warning(['DAT and X contain too few points for meaningful' ...
            ' analysis.'])
    outdat=dat;
    outx=x;
    return;
end

% Check & set STARTDAT, including output for data point 1:
if isempty(startdat)
    startdat=nan;
end
if isnan(startdat)
    startdat=dat(1);
    outqc(1)=1;
else
    if abs(startdat-dat(1))<=toldat
        startdat=dat(1);
        outqc(1)=1;
    else
        outqc(1)=0;
    end
end
if ~isscalar(startdat)
    error('STARTDAT must be scalar, NaN, or empty.')
end

% Main loop, checking for data points 2 through ll:
ii=2;
while (ii<=ll)
    if outqc(ii-1)==0
        if abs(dat(ii)-startdat)<=toldat
            outqc(ii)=1;
            startdat=dat(ii);
        else
            outqc(ii)=0;
        end
    else
        tmp=(dat(ii)-dat(ii-1))/(x(ii)-x(ii-1));
        if (tmp<ddatdx(1))|(tmp>ddatdx(2))
            outqc(ii)=0;
        else
            outqc(ii)=1;
            startdat=dat(ii);
        end
    end
    ii=ii+1;
end

outqc=logical(outqc);
outdat=dat;
outx=x;