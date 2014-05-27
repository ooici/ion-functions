% Creates a test dataset using values provided by Sunburst (see
% SAMI_P0132_110414.txt and the example code in PH_OOI_2012-10-24.m). This code
% takes the first 5 values of that test dataset, replicates them 3 times
% yielding 3 sets of 15 test values and then computes a new pH using 3 sets of
% salinity values (30, 32 and 34) for those three sets.
%
% The final test values are compared to the final pH estimate and the coverted
% thermistor and battery values.
%
% C. Wingard 2014-05-27

% set the raw data and default calibration coefficients.
raw = [
    10	3480086699	2530	3288	1435	2484	2228	3289	1435	2490	2234	3286	1439	2485	2235	3287	1439	2487	2238	3289	1439	2487	2238	3290	1432	2486	2208	3289	1265	2484	1630	3283	943	2482	797	3289	663	2485	336	3286	512	2481	178	3282	451	2479	136	3288	449	2486	142	3285	478	2480	173	3290	530	2488	238	3288	600	2481	329	3287	675	2483	443	3286	756	2486	582	3288	836	2480	732	3288	913	2484	894	3291	983	2484	1047	3288	1043	2485	1199	3287	1103	2479	1341	3291	1156	2489	1469	3287	1198	2480	1577	3292	1236	2480	1675	3286	1266	2480	1760	3287	1295	2485	1834	0	3374	2528;
    10	3480086999	2530	3285	1509	2480	2341	3284	1510	2484	2341	3288	1507	2482	2341	3289	1506	2480	2339	3294	1513	2487	2345	3285	1498	2477	2324	3286	1276	2478	1885	3285	857	2477	1140	3287	545	2478	639	3289	397	2479	423	3282	346	2482	354	3289	350	2480	358	3287	387	2482	399	3287	446	2484	478	3293	522	2482	581	3289	607	2480	702	3288	695	2483	842	3289	789	2479	984	3285	879	2479	1132	3288	962	2486	1283	3286	1035	2480	1406	3291	1104	2476	1532	3281	1160	2479	1640	3289	1215	2482	1741	3289	1259	2481	1827	3290	1295	2481	1900	3288	1326	2477	1964	0	3374	2530;
    10	3480087299	2527	3291	1581	2480	2457	3290	1585	2481	2462	3289	1583	2481	2461	3290	1582	2486	2462	3290	1581	2480	2460	3286	1579	2477	2449	3286	1331	2478	1979	3285	895	2476	1218	3289	571	2481	701	3284	413	2474	472	3290	364	2477	400	3293	369	2481	405	3284	404	2477	447	3288	467	2476	528	3290	541	2485	632	3290	624	2478	757	3290	721	2477	898	3287	814	2479	1049	3288	900	2482	1191	3292	995	2481	1339	3289	1069	2478	1475	3289	1140	2479	1601	3288	1206	2480	1718	3291	1258	2482	1815	3291	1307	2483	1907	3291	1344	2482	1984	3289	1381	2480	2049	0	3373	2526;
    10	3480087599	2530	3284	1630	2477	2538	3289	1633	2479	2539	3288	1635	2480	2541	3286	1636	2476	2543	3288	1635	2479	2543	3290	1632	2482	2530	3284	1391	2480	2085	3286	915	2475	1252	3286	572	2483	704	3289	411	2479	471	3286	368	2481	409	3282	378	2480	416	3290	422	2483	469	3281	482	2478	551	3289	561	2479	663	3289	651	2479	791	3291	743	2480	934	3294	841	2484	1087	3287	932	2473	1235	3285	1020	2479	1387	3285	1100	2479	1527	3292	1173	2486	1656	3288	1238	2482	1776	3281	1298	2476	1880	3285	1345	2478	1971	3285	1388	2478	2050	3290	1425	2479	2122	0	3374	2527;
    10	3480087899	2526	3283	1705	2478	2658	3290	1706	2474	2654	3287	1706	2479	2654	3289	1707	2483	2660	3287	1710	2480	2657	3288	1702	2482	2644	3286	1449	2481	2174	3289	951	2485	1303	3289	589	2487	728	3287	432	2478	489	3286	382	2473	415	3283	388	2480	426	3291	436	2480	482	3286	503	2480	573	3290	585	2474	684	3288	675	2484	823	3285	772	2477	972	3289	878	2480	1136	3289	973	2483	1295	3291	1065	2482	1446	3286	1151	2489	1604	3287	1220	2481	1732	3289	1293	2477	1853	3286	1352	2486	1959	3291	1402	2484	2059	3289	1448	2482	2141	3294	1484	2485	2211	0	3373	2529;
];
raw = repmat(raw,3,1);

ea434 = 17533;
ea578 = 101;
eb434 = 2229;
eb578 = 38502;
ind_slope = 0.9698;
ind_offset = 0.2484;
Salinity = [ones(1,5)*30, ones(1,5)*32, ones(1,5)*35]';

% Thermistor Temperatures and Battery Voltages
% Rt = (raw(:,3) ./ (4096 - raw(:,3))) * 17400;
% InvT = 0.0010183 + 0.000241*(log(Rt)) + 0.00000015 * (log(Rt)).^3;
% TempK = 1 ./ InvT;
% TempC1 = TempK - 273.15;

Rt = (raw(:,114) ./ (4096 - raw(:,114))) * 17400;
InvT = 0.0010183 + 0.000241 * (log(Rt)) + 0.00000015 * (log(Rt)).^3;
TempK = 1 ./ InvT;
TempC2 = TempK - 273.15;
%TempFinal = (TempC1 + TempC2) / 2;
TempFinal = TempC2;
clear Rt InvT TempK TempC1 TempC2

VBatt = raw(:,113) * 15 / 4096;

% Blank calculations
Blank434A = raw(:,5) ./ raw(:,4);
Blank578A = raw(:,7) ./ raw(:,6);
Blank434B = raw(:,9) ./ raw(:,8);
Blank578B = raw(:,11) ./ raw(:,10);
Blank434C = raw(:,13) ./ raw(:,12);
Blank578C = raw(:,15) ./ raw(:,14);
Blank434D = raw(:,17) ./ raw(:,16);
Blank578D = raw(:,19) ./ raw(:,18);

blank434 = (Blank434A + Blank434B + Blank434C + Blank434D) / 4;
blank578 = (Blank578A + Blank578B + Blank578C + Blank578D) / 4;
clear Blank434A Blank434B Blank434C Blank434D
clear Blank578A Blank578B Blank578C Blank578D

% Light measurements
Ref434 = raw(:, 20:4:108);
I434 = raw(:, 21:4:109);
Ref578 = raw(:, 22:4:110);
I578 = raw(:, 23:4:111);

% Absorbance
A434 = -log10(I434 ./ Ref434);
A578 = -log10(I578 ./ Ref578);
clear Ref434 I434 Ref578 I578

A434blank = -log10(blank434);
A578blank = -log10(blank578);
clear blank434 blank578

Abs434 = A434 - repmat(A434blank,1,23);
Abs578 = A578 - repmat(A578blank,1,23);
R = (A578 - repmat(A578blank,1,23)) ./ (A434 - repmat(A434blank,1,23));
pKa = (1245.69 ./ (TempFinal + 273.15)) + 3.8275 + (0.0021 * (35 - Salinity));
clear A434 A434blank A578 A578blank

% Molar absorptivities and initial pH estimations
Ea434 = repmat(ea434 - (26 * (TempFinal - 24.788)),1,23);
Ea578 = repmat(ea578 + (TempFinal - 24.788),1,23);
Eb434 = repmat(eb434 + (12 * (TempFinal - 24.788)),1,23);
Eb578 = repmat(eb578 - (71 * (TempFinal - 24.788)),1,23);
e1 = Ea578 ./ Ea434;
e2 = Eb578 ./ Ea434;
e3 = Eb434 ./ Ea434;

V1 = R - e1;
V2 = e2 - R .* e3;

HI = ((Abs434 .* Eb578) - (Abs578 .* Eb434)) ./ ((Ea434 .* Eb578) - (Eb434 .* Ea578));
I = ((Abs578 .* Ea434) - (Abs434 .* Ea578)) ./ ((Ea434 .* Eb578) - (Eb434 .* Ea578));
IndConc = HI + I;
pointpH = real(repmat(pKa,1,23) + log10(V1 ./ V2));
clear Abs434 Abs578 R pKa Ea434 Eb434 Ea578 Eb578 e1 e2 e3 V1 V2 HI I

% **********************************************************************
% Determine the most linear region of points for seawater calculation
% Skip first 5 points
IndConca = IndConc(:, 6:end);
Y = pointpH(:, 6:end);
X = 1:23-5;

% create arrays for vectorized computations used in sum of squares below.
% reflows 1D and 2D arrays into 2D and 3D arrays, respectively, shifting
% each "row" of the arrays by one value, allowing the sum of square
% calculations to be computed in a vectorized fashion.
nRec = size(pointpH, 1);
step = 7;  % number of points to use
count = step + 1;
nPts = size(X, 2) - step;
x = zeros(nRec, nPts, count);
y = zeros(nRec, nPts, count);
for i = 1:nPts
    for j = 1:nRec
        x(j, i, :) = X(i:i+step);
        y(j, i, :) = Y(j, i:i+step);
    end %for
end %for
clear X nPts i j

% compute the range of best fitting points, using array multiplications to
% determine the best fit via the correlation coefficient.
sumx = sum(x, 3);
sumy = sum(y, 3);
sumxy = sum(x .* y, 3);
sumx2 = sum(x.^2, 3);
sumy2 = sum(y.^2, 3);
sumxx = sumx .* sumx;
sumyy = sumy .* sumy;
ssxy = sumxy - (sumx .* sumy) ./ count;
ssx = sumx2 - (sumxx ./ count);
ssy = sumy2 - (sumyy ./ count);
r2 = ssxy.^2 ./ (ssx .* ssy);
clear sumx sumy sumxy sumx2 sumy2 sumxx sumyy ssxy ssx ssy

% Range of seawater points to use
[~,strt] = max(r2,[],2);  % Find start point of the best fit using best R-squared
cutoff1 = strt;           % Start point
cutoff2 = strt + step;    % End point
IndConcS = zeros(nRec, count);
pointpHS = zeros(nRec, count);
for i = 1:nRec
    IndConcS(i,:) = IndConca(i, cutoff1(i):cutoff2(i));
    pointpHS(i,:) = Y(i, cutoff1(i):cutoff2(i));
end %for
clear nRec IndConca Y strt step cutoff1 cutoff2 i

% ************************* Final pH Calcs ********************************
sumx = sum(IndConcS, 2);
sumy = sum(pointpHS, 2);
sumxy = sum(pointpHS .* IndConcS, 2);
sumx2 = sum(IndConcS.^2, 2);
sumy2 = sum(pointpHS.^2, 2);
xbar = mean(IndConcS, 2);
ybar = mean(pointpHS, 2);
sumxx = sumx .* sumx;
sumyy = sumy .* sumy;
ssxy = sumxy - (sumx .* sumy) ./ count;
ssx = sumx2 - (sumxx ./ count);
ssy = sumy2 - (sumyy ./ count);
slope = ssxy ./ ssx;
ph = ybar - slope .* xbar;

% correct values for indicator impurities if pH greater than 8.2
m = ph >= 8.2;
ph(m) = ph(m) * ind_slope + ind_offset;
