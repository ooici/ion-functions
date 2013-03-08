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