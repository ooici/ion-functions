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
        if ~isempty(varargin{3})
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