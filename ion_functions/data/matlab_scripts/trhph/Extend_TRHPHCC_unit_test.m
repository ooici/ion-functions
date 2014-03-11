%.. Extend_TRHPHCC_unit_test.m
%.. 28-Feb-2014 
%.. R. Desiderio desi@coas.oregonstate.edu
%
%.. the DPS unit tests for the TRHPHCC data product don't 
%..     test the resistivity branching, because the temperature
%..     values for these resistivities are outside the range
%..     of the Temp-[chl]-Cond surface, and therefore return NaNs.
%.. this code checks the unit test in the DPS (Chris Wingard has
%..     already done this and coded the DPA in python) and then
%..     provides non-trivial results for R and T values inside the
%..     surface domain.
%
%.. the matlab code from the Appendix A in the DPS is used (with
%.. cosmetic changes for clarity)

load Larson_2007surface.mat
%.. variables loaded:
%.. Tdat, a 2D array of temperature [degC]
%..     a 200 element row vector replicated to 200 rows;
%.. Sdat, a 2D array of chloride concentration [mol/kg]
%..     a 200 element column vector replicated to 200 columns;
%.. Cdat, a 2D array of experimentally determined conductivity values.

%.. the data limits are:
%
%                      min       max
% Tdat              103.9750  382.2280
% Sdat                0.0540    0.6340
% Cdat                0.4204    7.3838
%
% v_r = 1/Cdat        0.1354    2.3787

%.. original unit test data and results.
%.. note that all of the nan values could result from T < 103.795.
%.. therefore, v_r selection was and is not tested by these data.
%
%                           v_r1   v_r2   v_r3   temp  [chl]
original_unit_test_data = [0.906  4.095  4.095   11.8   nan; ...
                           0.890  4.095  4.095   15.9   nan; ...
                           0.891  4.095  4.095    3.2   nan; ...
                           0.184  0.915  4.064   67.7   nan; ...
                           0.198  1.002  4.095   75.8   nan; ...
                           0.172  0.857  4.082   97.5   nan; ...
                           0.183  0.926  4.076   95.0   nan; ...
                           0.233  1.182  4.072   96.2   nan; ...
                           0.146  0.747  3.634  116.8   195; ...
                           0.134  0.681  3.405  272.8   109; ...
                           0.131  0.673  3.293  325.8   128; ...
                           0.133  0.678  3.396  330.0   127; ...
                           0.135  0.681  3.409  333.4   129; ...
                           0.135  0.681  3.426  333.2   128];
 
%.. use only one T < 103.795 value to test T out-of-bounds;
%.. use a v_r2 value of 2.8 to test C = 1/2.8 < 0.4204 out-of-bounds.
%.. add a range of T values for the other v_r ranges to test the latter.
%.. note that the surface shape renders the 2nd recordset result as a nan.
%
%.. the new recordset results are denoted as negative numbers to
%.. indicate that these values are not from the DPS.
proposed_unit_test_data = [0.440  4.095  4.095  105.4   -59; ...
                           0.380  4.095  4.095  241.9   nan; ...
                           0.320  4.095  4.095  374.2   -60; ...
                           0.184  0.915  4.064  105.4  -175; ...
                           0.198  1.002  4.095  241.9   -71; ...
                           0.172  0.857  4.082  374.2  -132; ...
                           0.183  0.926  4.076   84.6   nan; ...
                           0.233  2.800  4.072  250.2   nan; ...
                           0.146  0.747  3.634  116.8   195; ...
                           0.134  0.681  3.405  272.8   109; ...
                           0.131  0.673  3.293  325.8   128; ...
                           0.133  0.678  3.396  330.0   127; ...
                           0.135  0.681  2.000  333.4  -239; ...
                           0.135  0.681  1.000  333.2  -501];

%.. vertically concatenate original and proposed data 
unit_test_data = [original_unit_test_data; proposed_unit_test_data];

V_R1 = unit_test_data(:,1);
V_R2 = unit_test_data(:,2);
V_R3 = unit_test_data(:,3);
TEMP = unit_test_data(:,4);
chl_mmol_kg = unit_test_data(:,5);

nrecords = length(V_R1);
test_out = nan(nrecords,1);
for ii = 1:nrecords
    
    v_r1 = V_R1(ii);
    v_r2 = V_R2(ii);
    v_r3 = V_R3(ii);
    temp = TEMP(ii);
        
    %.. find the best resistivity range
    if     (v_r2  < 0.75)
        v_r = v_r3/5.0;
    elseif (v_r2 >= 0.75) && (v_r2 < 3.90)
        v_r = v_r2;
    else
        v_r = 5.0*v_r1;
    end
    
    %.. conductivity
    cond = 1.0/v_r;
    
    %.. set up [chl] vector
    Scurve = linspace(nanmin(Sdat(:)),nanmax(Sdat(:)),100);
    %.. constant temp vector
    Tcurve = Scurve*0 + temp;
    %.. interpolate a cond vector from temp and [chl] onto Cdat
    Ccurve = interp2(Tdat, Sdat, Cdat, Tcurve, Scurve);
    
    if isfinite(Ccurve)
        %.. interpolate conductivity msrmnt onto Ccurve to get [chl]
        S_mol_kg =interp1(Ccurve, Scurve, cond);
        Cl_mmol_kg = 1000.0 * S_mol_kg;
    else
        Cl_mmol_kg = nan;
    end
    
    test_out(ii) = round(Cl_mmol_kg);
end

%.. send comparison results to screen
disp('Original unit test:');
disp('First column, expected from DPS; 2nd column, calculated.')
disp([chl_mmol_kg(1:nrecords/2) test_out(1:nrecords/2)]);
disp(' ');
disp('Proposed unit test data:');
proposed_unit_test_data(:,end) = test_out(nrecords/2+1:end);
disp(proposed_unit_test_data);


