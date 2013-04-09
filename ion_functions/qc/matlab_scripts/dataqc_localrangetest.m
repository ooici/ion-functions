function out = dataqc_localrangetest(dat,z,datlim,datlimz)

% DATAQC_LOCALRANGETEST: Data quality control algorithm testing
%   if measurements fall into a user-defined valid range.
%   This range is not constant but varies with measurement location.
%   Returns 1 for presumably good data and 0 for data presumed bad.
%
% Time-stamp: <2012-07-06 15:33:41 mlankhorst>
%
% USE CASE SCENARIO: A time series of measurements DAT is given at a
%   time-varying altitude Z, e.g. air temperature from a rising
%   weather balloon. The purpose of this routine is to identify
%   data points that fall within an expected range, which varies
%   with altitude (e.g. at the ground, accept temperatures -20..35
%   degrees C, but at 10000m altitude accept -80..-40 degrees C).
%
%   Z can have more than one dimension, in which case the locations
%   are interpreted as being in a multi-dimensional space. Use this
%   feature e.g. if your temperature ranges are further categorized
%   by season/month, in which case one dimension might be the
%   altitude, and the other the time of year.
%
% USAGE:    OUT=dataqc_localrangetest(DAT,Z,DATLIM,DATLIMZ);
%
%           OUT: Boolean, 0 if value is outside range, else 1.
%
%           DAT: Input dataset, a numeric real scalar or column
%                vector.
%           Z: Location of measurement DAT. Must have same number
%              of rows as DAT and same number of columns as DATLIMZ.
%           DATLIM: Two-column matrix with the minimum (column 1)
%                   and maximum (column 2) values considered valid.
%           DATLIMZ: Matrix with the locations where DATLIM is
%                    given. Must have same number of rows as DATLIM and
%                    same number of columns as Z.
%

checkinput(dat,'DAT');
checkinput(z,'Z');
checkinput(datlim,'DATLIM');
checkinput(datlimz,'DATLIMZ');

[numlim,ndim]=size(datlimz);

[tmp1,tmp2]=size(datlim);
if tmp1~=numlim
    error('DATLIM and DATLIMZ must have same number of rows.')
end
if tmp2~=2
    error('DATLIM must have exactly 2 columns.')
end

[num,tmp2]=size(z);
if tmp2~=ndim
    error('Z must have same number of columns as DATLIMZ.')
end

if ~isvector(dat)
    Data Product Specification for Local Range Test
    Ver 1-00 1341-10005 Page 7 of 7
    error('DAT must be vector.');
end
dat=dat(:);
if num~=length(dat)
    error('Length of DAT must match number of rows in Z.')
end

if ~all(datlim(:,2)>datlim(:,1))
    warning(['Second column values of DATLIM should be greater than' ...
        ' first column values.'])
end

if ndim==1
    lim1=interp1(datlimz,datlim(:,1),z);
    lim2=interp1(datlimz,datlim(:,2),z);
else
    F=TriScatteredInterp(datlimz,datlim(:,1));
    lim1=F(z);
    F=TriScatteredInterp(datlimz,datlim(:,2));
    lim2=F(z);
end

ff=find((isnan(lim1))|(isnan(lim2)));
lim1(ff)=max(datlim(:,2));
lim2(ff)=min(datlim(:,1));

out=(dat>=lim1)&(dat<=lim2);
end %function

function checkinput(x,txt)

if ~isnumeric(x)
    error('%s must be numeric.',txt)
end
if ~ismatrix(x)
    error('%s must be a matrix.',txt)
end
if ~all(isreal(x(:)))
    error('%s must be real.',txt)
end
end %function
