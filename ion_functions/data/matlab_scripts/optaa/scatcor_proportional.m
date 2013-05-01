function  [scat_cor_apg] = scatcor_proportional(acs_a,acs_c,a_wl,c_wl, ...
                                                ref_wl,ag_ref_wl)
% Corrects for absorption tube scattering error using proportional
% method (Zaneveld et al. 1994, method 3).  ref_wl is 715.
% ag_ref_wl is the dissolved absorption at the reference wavelength
% (measured with a pre-filter).  It should be a column vector with
% the same length as acs_a and acs_c. Using this ensures any small
% absorption in the near IR is not subtracted. Make sure there are no
% negative ag values (use something like the routine below). If this
% measurement doesn't exist, insert 0.

% Find a(ref_wl) for scattering correction
ref_wl = 715
[near, ind_a_ref_wl] = find_nearest(ref_wl,a_wl);

% Ensure apg(ref_wl) is not negative
apg_ref_wl = acs_a(:,ind_a_ref_wl);
i = find(apg_ref_wl < 0);
apg_ref_wl(i) = 0;

% Proportional scattering correction to acs
% First interpolate c values to match a wavelengths.
acs_ci = pchip(c_wl,acs_c,a_wl);

% Find scattering error at ref_wl, ae_ref_wl
ae_ref_wl = apg_ref_wl - ag_ref_wl;

% Calculate apg using proportional method
 [scat_cor_apg] = acs_a - repmat(ae_ref_wl,1,length(a_wl)).*...
                  ((acs_ci - acs_a)./repmat(acs_ci(:,ind_a_ref_wl)...
                  - apg_ref_wl,1,length(a_wl)));
