% DATAQC_PROPAGATEFLAGS     Propagate "bad" qc flags (from an
%   arbitrary number of source datasets) to another (derived)
%   dataset.
%
% Time-stamp: <2010-12-16 18:27:54 mlankhorst>
%
% USAGE:    OUTFLAG=DATAQC_PROPAGATEFLAGS(INFLAGS);
%
%   INFLAGS: An M-by-N boolean matrix, where each of the M rows
%            contains flags of an independent data set such that
%            "0" means bad data and "1" means good data.
%   OUTFLAG: A 1-by-N boolean vector that contains 1 where all of
%            the INFLAGS are 1, and 0 otherwise.
%
% EXEMPLAR USE CASE: Consider data from an oceanographic
%   CTD (conductivity-temperature-pressure) instrument. From these
%   three time series, you want to compute salinity. If any of the
%   three source data (conductivity, temperature, pressure) is of
%   bad quality, the salinity will be bad as well. You can feed
%   your QC assessment of the former three into this routine, which
%   will then give you the combined assessment for the derived
%   (here: salinity) property.
%
% EXAMPLE SCREENSHOT:
%
%
%   INFLAGS =
%
%       0 0 1 1
%       1 0 1 0
%
%   >> dataqc_propagateflags(INFLAGS)
%
%   ans =
%
%       0 0 1 0
%
function out=dataqc_propagateflags(in);

if ~islogical(in)
    error('INFLAGS must be of type ''logical''')
end

si=size(in);
if (length(si))~=2
    error('INFLAGS must be two-dimensional array')
end

out=all(in,1);