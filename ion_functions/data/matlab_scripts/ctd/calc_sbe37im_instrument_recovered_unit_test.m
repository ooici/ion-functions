%.. calc_sbe37im_instrument_recovered_unit_test.m
%
%*********************************************************************
%*********************************************************************
%.. THE UNIT TEST INPUTS CAN BE COPY\PASTED FROM THE CODE.
%.. THE UNIT TEST TEST VALUES CAN BE COPY\PASTED FROM THE CODE OUTPUT.
%*********************************************************************
%*********************************************************************
%
%.. the TEMPWAT_L1, CONDWAT_L1, and PRESWAT_L1 DPSs did not 
%.. contain unit tests for INSTRUMENT_RECOVERED data from
%.. the SBE 37-IM ctd (CTDMO series G,H,Q,R).
%
%.. these unit tests are constructed using information from sbe37IM
%.. sensor cal sheets, an application note from the confluence IDD
%.. page for ctdmo_ghqr (the link to this pdf, "Explanation on raw output
%.. format for 37-IMv3.1.pdf", is in the Record Structure section
%.. underneath the 'Science Data' heading), and checked against SBE
%.. data processed values contained in a spreadsheet created by Aidan
%.. Alai at WHOI ("Instrument_Recovered_Data DPA verification.xlsx"),
%.. which also contains the rawdata values. 

%.. Written by:
%.. Russell Desiderio, Oregon State University
%.. desi@coas.oregonstate.edu
%.. 2016-Jun-15


%.. rawdata; rows are from Aidan's 'Corrected DPA Data' sheet.
%.. row #:     274,   21346,     789,     590
temp_L0 = [ 366964,  499888,  465784,  500403];  % temperature
pres_L0 = [ 533152,  571309,  632465,  828170];  % pressure
ptmp_L0 = [   1608,    1452,    1471,    1453];  % pressure thermistor
cond_L0 = [1564991, 1457279, 1484332, 1462659];  % conductivity 

%.. SBE-calculated values; rows are from the 'SBE Processed Data' sheet.
%.. row #:       274,    21346,      789,      590
temp_SBE = [ 10.9818,   3.8488,   5.4520,   3.8255];  % [C]
pres_SBE = [  16.159,  134.950,  325.258,  933.882];  % [db]
cond_SBE = [3.891373, 3.240767, 3.399795, 3.272396];  % [S/m]

%.. for clarity use variable names as on SBE calsheets
%>>>>>>>>>>>>>>>>>>>>>>>> TEMPERATURE >>>>>>>>>>>>>>>>>>>>>>>>
%
%<<< the coded temperature calculation is the same in all three of:
%<<<     (1) aidan's spreadsheet,
%<<<     (2) as denoted in the appnote pdf, and 
%<<<     (3) in the SBE 37IM temperature calibration sheet.
%
%.. cal coeffs
a0 = -1.179278E-04;
a1 = 3.097942E-04;
a2 = -4.688854E-06;
a3 = 2.081274E-07;

n = temp_L0;
ln_n = log(n);
temp_mat = 1.0 ./ (a0 + ln_n .* (a1 + ln_n .* (a2 + ln_n * a3)) ) - 273.15;

temperature = [temp_mat; temp_SBE; temp_mat-temp_SBE];

%>>>>>>>>>>>>>>>>>>>>>>>> PRESSURE >>>>>>>>>>>>>>>>>>>>>>>>
%
%<<< the coded pressure calculation, which differs from the SBE processed
%<<< values at the mm scale (which is at the precision of the SBE values),
%<<< is the same as:
%<<<     (1) denoted in the appnote pdf, and 
%<<<     (2) in the SBE 37IM pressure calibration sheet.
%<<<
%<<< Aidan's calculation divides the pressure thermistor rawdata by
%<<< 13107 (as is done for the SBE16+ DPA). However, the agreement to
%<<< the SBE calculated values is much worse (up to 0.7 db at 1000 db).
%<<< Also, the pressure thermistor temperature value so calculated is
%<<< unphysical (-69.5 C = ptempa0), whereas this value as coded below is
%<<< within a couple of degrees of the temperature reading of the CTD's
%<<< primary temperature sensor.
%
%.. cal coeffs
pa0 = 1.202594e-1;
pa1 = 4.514834e-3;
pa2 = -1.091899e-11;
ptempa0 = -6.953022e1;
ptempa1 = 5.115592e-2;
ptempa2 = -3.918145e-7;
ptca0 = 5.247204e5;
ptca1 = 9.617295e-1;
ptca2 = 6.296724e-3;
ptcb0 = 2.498163e1;
ptcb1 = -2.75e-4;
ptcb2 = 0;

y = ptmp_L0;  % do not divide by 13107
t = ptempa0 + y .* (ptempa1 + y * ptempa2);
x = pres_L0 - ptca0 - t .* (ptca1 + t * ptca2);
n = x * ptcb0 ./ ( ptcb0 + t .* (ptcb1 + t * ptcb2) );
psia = pa0 + n .* (pa1 + n * pa2);
pres_mat = psia * 0.689475729 - 10.1325;

pressure = [pres_mat; pres_SBE; pres_mat-pres_SBE];

%>>>>>>>>>>>>>>>>>>>>>>>> CONDUCTIVITY >>>>>>>>>>>>>>>>>>>>>>>>
%.. the conductivity calculation also depends on (as calculated above):
%..     in situ temperature [degC]
%..     in situ pressure    [dbar]
%
%<<< the coded conductivity calculation is the same as denoted in the
%<<< appnote pdf, which states that the raw conductivity values need to
%<<< be divided by 256 before the formula in the SBE 37IM conductivity
%<<< calibration sheet is used.
%<<<
%<<< Aidan skipped the factor containing wbotc in his calculation of
%<<< conductivity. This causes disagreement in the 5th decimal place;
%<<< including it pushes the disagreement out to the 7th decimal place,
%<<< which is at the precision of the SBE values.
%
%.. cal coeffs
g = -9.899853E-01;
h = 1.314100E-01;
i = -4.181710E-04;
j = 4.723872E-05;
cpcor = -9.570000E-08;
ctcor = 3.250000E-06;
wbotc = 4.842900E-07;

f = (cond_L0/256.0/1000.0) .* sqrt(1.0 + wbotc * temp_mat); 
numerator = g + f .* f .* (h + f .* (i + f * j));
denominator = 1.0 + ctcor * temp_mat + cpcor * pres_mat;
cond_mat = numerator ./ denominator;

conductivity = [cond_mat; cond_SBE; cond_mat-cond_SBE];

%>>>>>>>>>>>>>>>>>>>>>>>> OUTPUT COMPARISON >>>>>>>>>>>>>>>>>>
format_string = '%14.8f%14.8f%14.8f\n';

fprintf('\n');
disp('TEMPERATURE [C]:   matlab, Seabird, matlab - SBE');
fprintf(format_string, temperature);
fprintf('\n');

fprintf('\n');
disp('PRESSURE [dbar]:   matlab, Seabird, matlab - SBE');
fprintf(format_string, pressure);
fprintf('\n');

fprintf('\n');
disp('CONDUCTIVITY [S/m]:   matlab, Seabird, matlab - SBE');
fprintf(format_string, conductivity);
fprintf('\n');


%.. OUTPUT:
% 
% 
% >>calc_sbe37im_instrument_recovered_unit_test
% 
% TEMPERATURE [C]:   matlab, Seabird, matlab - SBE
%    10.98176998   10.98180000   -0.00003002
%     3.84878502    3.84880000   -0.00001498
%     5.45204077    5.45200000    0.00004077
%     3.82553871    3.82550000    0.00003871
% 
% 
% PRESSURE [dbar]:   matlab, Seabird, matlab - SBE
%    16.16196123   16.15900000    0.00296123
%   134.95247848  134.95000000    0.00247848
%   325.26072294  325.25800000    0.00272294
%   933.88494132  933.88200000    0.00294132
% 
% 
% CONDUCTIVITY [S/m]:   matlab, Seabird, matlab - SBE
%     3.89137261    3.89137300   -0.00000039
%     3.24076705    3.24076700    0.00000005
%     3.39979529    3.39979500    0.00000029
%     3.27239616    3.27239600    0.00000016




