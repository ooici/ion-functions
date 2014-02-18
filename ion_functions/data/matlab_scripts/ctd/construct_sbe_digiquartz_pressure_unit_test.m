%.. construct_sbe_digiquartz_pressure_unit_test.m
%
%*********************************************************
%*********************************************************
%.. THE UNIT TEST INPUTS AND OUTPUT ARE LOCATED IN THE   *
%.. DOCUMENTATION AT THE VERY END OF THIS M-FILE.        *
%*********************************************************
%*********************************************************
%
%.. the PRESWAT_L1 DPS did not include a unit test for 
%.. the Paroscientific digiquartz pressure sensors used
%.. with SBE 16Plus CTDs.
%
%.. this unit test is constructed using DPS artifacts found on
%.. alfresco (See: Company Home >> OOI >> REFERENCE >> Data
%.. Product Specification Artifacts >> 1341-00020_PRESWAT >>
%..     (1) P124969A-corrected.pdf (Paroscientific cal sheet for
%..         pressure sensor model 2200A-219, SN 124969, pressure
%..         range 0 to 200 psia (127 dbar)).
%..     (2) PRESWAT_SeaBird_16PlusV2_2009.pdf (User manual).
%
%.. Written by:
%.. Russell Desiderio, Oregon State University
%.. desi@coas.oregonstate.edu
%.. 2014-Feb-03

%.. cal coeffs from cal sheet:
%.. 19-Oct-2012 coeffs, not 29-Nov-2012 update.
C1 = 991.3651;
C2 =   1.0136e-05;
C3 =  -1.18210e-04;
D1 =   0.031072;
D2 =   0.0;
T1 =  27.67412;
T2 =  -1.08033e-04;
T3 =   1.03670e-06;
T4 =   1.68749e-09;
T5 =   0.0;

%.. first, for temperature U = 21.0 deg_C,
%.. calculate derived cal coeffs C, D, and T0  
%.. and check against cal sheet values for this U.
U  = 21.0;
C  = C1 + U*(C2 + U*C3);
D  = D1 + U*D2;
T0 = T1 + U*(T2 + U*(T3 + U*(T4 + U*T5) ) );

disp('*******************************************');
disp(' Back-calculation to make unit test for ');
disp(' SBE digiquartz pressure sensor.');
disp(' ');
disp('Check calculated derived calcoeffs [C D T0]:'); 
fprintf('%10.5f  %9.7f %10.6f\n', [C D T0]);
disp('Cal sheet values [C D T0]:');
disp(' 991.3132   0.031072   27.67232');
disp(' ');
disp(' ');
%.. the calculated values match those on the cal sheet.
%.. retain full precision in calculations to come.

disp(' Back Calculations for unit test:');
%.. For a complete unit test, need to start with SBE
%.. 16Plus V2 pressure thermistor reading in counts
%.. and pressure reading in counts.
%
%.. back-calculate SBE counts tc from pressure
%.. thermistor for U = 21.0 deg_C.
%
%.. the forward calculation from the DPS:
%.. U  = ( 23.7 * (tv+9.7917) ) - 273.15;
%.. tv = tc/13107;
U  = 21.0;
tc = 13107 * ( (U + 273.15)/23.7 - 9.7917);
disp('pressure thermistor counts tc for U=21.0 degC:');
fprintf('tc = %10.1f\n', tc);
%.. answer is tc = 34336.3 counts.
%.. in the SBE data stream, this is coded for in 4 hex
%.. bytes, so max value for tc is 16^4-1 = 65535. the
%.. calculated value is consistent with this.

%.. back-calculate SBE pressure reading in counts,
%.. assuming a pressure of 50 dbar.
%
%.. forward calculation from cal sheet agrees with DPS:
%
%.. factor = 1.0 - (T0/T)^2
%.. Pabs = C * factor * (1.0 - D*factor)
%.. Pdb  = Pabs * 0.689475729 - 10.1325
%
%.. find T, the pressure period, from factor.
Pdb = 50.0;
Pabs = (Pdb + 10.1325)/0.689475729; 
%.. need to solve a quadratic equation, the positive root < 1 
%.. is the one we want. however, instead of the usual equation:
%   factor = ( 1 - sqrt(1 - 4.0*Pabs*D/C) )/2.0/D; i will
%.. use the more robust algorithm.
temporary = -0.5 * ...
    (-1 - sign(-1)*sqrt(1 - 4.0*Pabs*D/C));
factor = temporary/D;
T = T0/(sqrt(1.0 - factor));
disp('pressure period for U=21.0 degC, Pdb=50 dbar:');
fprintf('T = %10.6f usec\n', T);
%
%.. the back-calculation to SBE pressure counts in
%.. the DPS appears to be incorrect: the SBE 16Plus V2 
%.. manual divides the raw counts by 256 to convert 
%.. to frequency, while the DPS omits this step.
%
%.. pf = pc/256.0; % pf in [Hz]; omitted in DPS
%.. T = 10^6/pf;   % T is the period in usec

pc = 256.0 * 10^6/T;
disp('pressure counts for U=21.0 degC, Pdb=50 dbar:');
fprintf('pc = %10.1f\n', pc);
%.. answer is pc = 8,833,628.8 counts.
%.. in the SBE data stream, this is coded for in 6 hex
%.. bytes, so max value for pc is 16^6-1 = 16,777,215.
%.. the calculated value is consistent with this.

%.. counts in the actual SBE data stream will be integers.
%
%.. so, for the unit test, specify as inputs:
%.. tc =   34336 counts.
%.. pc = 8833629 counts.
%
%.. forward calculate the corresponding pressure:
%.. follow the DPS, except also include the divisor
%.. 256 to convert pressure counts to pressure Hz.

%.. start completely clean
clear all
%.. cal coeffs
C1 = 991.3651;
C2 =   1.0136e-05;
C3 =  -1.18210e-04;
D1 =   0.031072;
D2 =   0.0;
T1 =  27.67412;
T2 =  -1.08033e-04;
T3 =   1.03670e-06;
T4 =   1.68749e-09;
T5 =   0.0;
%.. inputs
tc =   34336; % counts.
pc = 8833629; % counts.
%.. calculate thermistor temperature U from tc
tv = tc/13107.0; % [volts]
U  = 23.7 * (tv+9.7917) - 273.15;
%.. calculate intermediate cal coeffs C, D, T0
C  = C1 + U*(C2 + U*C3);
D  = D1 + U*D2;
T0 = T1 + U*(T2 + U*(T3 + U*(T4 + U*T5) ) );
%.. calculate the pressure frequency in Hz
pf = pc/256.0;  % [Hz]
T  = 10^6/pf;   % pressure period [usec]
factor = 1.0 - (T0/T)^2;
Pabs = C * factor * (1.0 - D*factor);
Pdb  = Pabs * 0.689475729 - 10.1325;
disp(' ');

disp(' Forward calculation of pressure.')
disp(' (Should be 50 dbar):');
disp([' calculated pressure [dbar] = ' ...
    num2str(Pdb,'%10.6f')]);
disp('*******************************************');
disp(' ');
disp(' ');
disp(' UNIT TEST SBE inputs:');
fprintf(' thermistor counts tc =%10u.\n', tc);
fprintf(' pressure   counts pc =%10u.\n', pc);
disp(' ');
disp(' UNIT TEST pressure output:');
fprintf(' calculated pressure =%10.6f dbar.\n', Pdb);
disp(' ');
disp(' ');
disp('*******************************************');
%***********************************************************
% RESULT (copy-pasted from screen output above):           *
%                                                          *
%  UNIT TEST SBE inputs:                                   *
%  thermistor counts tc =     34336.                       *
%  pressure   counts pc =   8833629.                       *
%                                                          *
%  UNIT TEST pressure output:                              *
%  calculated pressure = 49.999967 dbar.                   *
%                                                          *
% VALID for the following calibration coefficients:        *
%                                                          *
% C1 = 991.3651;                                           *
% C2 =   1.0136e-05;                                       *
% C3 =  -1.18210e-04;                                      *
% D1 =   0.031072;                                         *
% D2 =   0.0;                                              *
% T1 =  27.67412;                                          *
% T2 =  -1.08033e-04;                                      *
% T3 =   1.03670e-06;                                      *
% T4 =   1.68749e-09;                                      *
% T5 =   0.0;                                              *
%***********************************************************
