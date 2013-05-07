function DO = dosv(Pt,T,S,P,PDENS,csv)
% function DO = dosv(Pt,T,S,P,PDENS,csv)
%
% Calculate DO from optode bphase and temperature
% using the modified Stern-Volmer equation.
%
% INPUTS:
%
% Pt : Optode bphase (degrees)
% T : Temperature (deg C)
% S : Salinity
% P : Pressure (dbar)
% PDENS : Potential density (kg/m^3)
% csv : Stern-Volmer coefficients (7-element vector)
%
% Pt,T,S,P and PDENS may be vectors or matrices.
%
% OUTPUT:
%
% DO : Dissolved oxygen (umol/kg)
%
% rev 04/16/2012 R DRUCKER
% UNIVERSITY OF WASHINGTON

% Check input variables:
if ~isequal(size(Pt),size(T),size(S),size(P),size(PDENS))
    error('Input variable sizes do not match')
end
if length(csv)~=7
    error('Stern-Volmer coefficients incorrect size')
end
[rows,cols] = size(Pt);
Pt = Pt(:);
T = T(:);
S = S(:);
P = P(:);
PDENS = PDENS(:);

% Calculate DO using Stern-Volmer:
Ksv = csv(1) + csv(2)*T + csv(3)*(T.*T);
P0 = csv(4) + csv(5)*T;
Pc = csv(6) + csv(7)*Pt;
DO = ((P0./Pc) - 1) ./ Ksv;

% Convert from volume to mass units:
DO = 1000*DO./PDENS;

% Pressure correction:
pcomp = 1 + (0.032*P)/1000;
DO = pcomp.*DO;

% Salinity correction:
S0 = 0;
ts = log((298.15-T)./(273.15+T));
B = [-6.24097e-3
    -6.93498e-3
    -6.90358e-3
    -4.29155e-3];
C0 = -3.11680e-7;
scomp = exp((S-S0).*(vandermonde(ts,3)*B)+C0*(S.^2-S0.^2));
DO = scomp.*DO;

% Reshape output:
DO = reshape(DO,rows,cols);
function V = vandermonde(x,n)
    x = x(:);
    m = length(x);
    V = NaN(m,n);
    V(:,1) = ones(m,1);
    for k = 2:n+1
        V(:,k) = V(:,k-1).*x;
    end