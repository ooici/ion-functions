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