function  [corrected_data] = acs_corr(data)

% Corrects ac9 for effects of temperature and salinity.

load TS.cor ;
% Temperature and salinity coefficients file.

Tref = 12.0;
% Reference temperature (temperature at time of calibration).

acs_cal = load('acs010_051006.cor' );
% This is an instrument-specific file supplied by the vendor.

c_wl = acs_cal(:,1)';
a_wl = acs_cal(:,3)';

acs_Wcal = [acs_cal(:,2)'; acs_cal(:,4)'];
% Transform calibration file to 2x84 (varies with acs) matrix; first
% row is a, second row is c.

% Find temperature and salinity coefficients at wavelengths matching
% acs.
t = zeros(2,length(a_wl);
for  k = 1:length(a_wl)
    [near_c, t(1,k)] = find_nearest(c_wl(k),TS(:,1));
    [near_a, t(2,k)] = find_nearest(a_wl(k),TS(:,1));
end
acs_c_tempcor = TS(t(:,1),2);
acs_a_tempcor = TS(t(:,2),2);
acs_c_salcor = TS(t(:,1),3);
acs_a_salcor = TS(t(:,2),4);

% Column assignments must be corrected for each dataset
temp = data(:,3); % CTD temperature 
sal = data(:,5); % CTD salinity column
% In this case, “data” is a merged CTD and ac-s data file

acs_c = data(:,12:95); % c (in ascending wavelengths)
acs_a = data(:,96:179); % a (in ascending wavelengths)

% Correct data to reference temp., 0 PSU, subtract Pure Water Offset.
for  z = 1:length(a_wl);
    acs_a_cor(:,z) = acs_a(:,z) - ((temp - Tref) * acs_a_tempcor(z)) ...
                     - (sal * acs_a_salcor(z)) - repmat(acs_Wcal(1,z),r,1);
    acs_c_cor(:,z) = acs_c(:,z) - ((temp - Tref) * acs_c_tempcor(z))...
                     - (sal * acs_c_salcor(z)) - repmat(acs_Wcal(2,z),r,1);
end

[corrected_data] = [acs_a_cor acs_c_cor];