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