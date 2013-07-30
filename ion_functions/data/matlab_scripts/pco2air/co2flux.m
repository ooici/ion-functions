% CO2FLUX Flux of CO2 across the air-sea interface
%
% Time-stamp: <2012-03-28 11:53:24 mlankhorst>
%
% Compute flux of carbon dioxide between air and sea from partial
% pressures of CO2 in air and sea, sea water temperature and
% salinity, and wind speed via Bulk formula.
%
% USAGE:    f=co2flux(pco2w,pco2a,u10,t,s);
%
%   f:      CO2 flux, positive from sea to air [mol m-2 s-1]
%   pco2w:  Partial pressure of CO2 in sea water [microatm]
%   pco2a:  Partial pressure of CO2 in air [microatm]
%   u10:    Instantaneous wind speed at 10m above sea level
%           [m s-1]
%   t:      Sea surface temperature [deg C] (note: difference
%           between temperature scales of 1968 and 1990 is
%           negligible for this algorithm)
%   s:      Sea surface salinity [g/kg] (note: difference between
%           absolute/practical salinity is negligible for this
%           algorithm)
%
% References:
%
%   R.F. Weiss (1974): "Carbon Dioxide in Water and Seawater: The
%   Solubility of a Non-Ideal Gas". Marine Chemistry, vol. 2,
%   pp. 203-215.
%
%   R. Wanninkhof (1992): "Relationship Between Wind Speed and Gas
%   Exchange Over the Ocean". Journal of Geophysical Research,
%   vol. 97, no. C5, pp. 7373-7382.
%
%   C. Sweeney, E. Gloor, A. R. Jacobson, R. M. Key, G. McKinley,
%   J. L. Sarmiento, R. Wanninkhof (2007): "Constraining global
%   air-sea gas exchange for CO2 with recent bomb 14C
%   measurements". Global Biogeochemical Cycles, vol. 21,
%   no. GB2015.
%
function f=co2flux(pco2w,pco2a,u10,t,s);

% convert micro-atm to atm:
pco2a=pco2a./1e6;
pco2w=pco2w./1e6

% Compute Schmidt number (after Wanninkhof, 1992, Table A1):
Sc=2073.1-(125.62.*t)+(3.6276.*(t.^2))-(0.043219.*(t.^3));

% Compute gas transfer velocity
% (after Sweeney et al., 2007, Fig. 3 and Table 1):
k=0.27.*(u10.^2).*sqrt(660./Sc);

% convert cm h-1 to m s-1
k=k./(100*3600);

% Compute absolute temperature:
T=t+273.15;

% Compute solubility (after Weiss, 1974, Eqn. 12 and Table I).
% Note that there are two versions, one for units per volume and
% one per mass. Here, the volume version is used.
% mol atm-1 m-3
K0=1000.*exp(-58.0931+(90.5069.*(100./T))+(22.2940.*log(T./100))+ ...
    s.*(0.027766-(0.025888.*(T./100))+ ...
    (0.0050578.*((T./100).^2))));

% mol atm-1 kg-1
% K0=exp(-60.2409+(93.4517.*(100./T))+(23.3585.*log(T./100))+ ...
%   s.*(0.023517-(0.023656.*(T./100))+ ...
%       (0.0047036.*((T./100).^2))));

% Compute flux (after Wanninkhof, 1992, eqn. A2):
f=k.*K0.*(pco2w-pco2a);