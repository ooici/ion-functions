function [T_s, T_c, T_u, T_lc, T] = resT(V_s, V_c, a, b, c, d, e)

% This function computes Temperature from raw Temp-Resistivity Probe
% measurements using a thermistor and a thermocouple. Initial inputs are
% the raw voltages for the thermistor (V_s) and thermocouple (V_c) and
% laboratory calibration curve coefficients a,b,c,d,e.

% The following coefficients are derived from a 4th degree polynomial fit
% of a calibration curve, where ax^4+bx^3+cx^2+dx+e. These will be
% entered as input metadata. Use the following for testing purposes.
% a=1.98e-9;
% b=-2.45e-6;
% c=9.28e-4;
% d=-0.0888;
% e=.731;

T_s = 27.50133 - 17.2658*V_s + 15.83424/V_s %raw thermistor temperature
T_c=244970*V_c/1000 %raw thermocouple temperature
T_u= T_s + T_c %uncorrected total temperature

T_lc= a*(T_u)^4 + b*(T_u)^3 + c*(T_u)^2 + d*(T_u) + e; %correction
%based on laboratory calibration

T= T_u + T_lc %final, corrected temperature at sensor tip