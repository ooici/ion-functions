function [V_R, C, S, Cl] = Chloride(V_R1,V_R2,V_R3,T,Tdat,Sdat,Cdat);
%This function uses resistivity and temperature inputs to calculate Chloride
%concentration based on the surface data developed in Larson et al, 2007.

%load Larson_2007surface.mat %this loads Tdat, Cdat, Sdat
if V_R2 < 0.75
    V_R = V_R3/5;
elseif 0.75 <= V_R2 < 3.90
    V_R = V_R2;
else
    V_R = V_R1*5;
end %if

C=1/V_R; %conductivity from resistivity

%extract a curve of constant temperature out of the data surface
Scurve = linspace(min(min(Sdat)),max(max(Sdat)),100);
Tcurve = Scurve*0+T;
Ccurve = interp2(Tdat,Sdat,Cdat,Tcurve,Scurve);

if isfinite(Ccurve)
    %now interpolate onto the Scurve/Ccurve
    S = interp1(Ccurve,Scurve,C);
    Cl = (S*1000);
else
    S = NaN;
    Cl = NaN;
end %if