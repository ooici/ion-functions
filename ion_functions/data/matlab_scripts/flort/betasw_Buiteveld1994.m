function[theta,betasw,bsw]=betasw_Buiteveld1994(lambda,sal,T);
%computes pure seawater scattering functions based on:
%Buiteveld et al. (1994). SPIE Ocean Optics XII, 2258:174-183.
%For refs for physical expressions below, see Buiteveld et al.
%coded by Michael Twardowski, 2005
%email: mtwardo@WET Labs, Inc2.com
%lambda in nm
%sal in Practical Salinity Units
%T in degC
%when sal=0, bsw is the total scattering of pure water
%for backscattering coefficients, divide total scattering by 2
theta_increment=0.01; %can vary this
theta=0:theta_increment:180;
rad=theta/180*pi;
k=1.38054e-23; %1.38054e-23 Boltzmann constant
depolar_ratio=0.051; %0.051 from expt data of Farinato and Roswell (1975)
n_wat=1.3247+3.3e3*lambda^-2-3.2e7*lambda^-4-2.5e-6*T2; %from Mcneil (1977)
%n_sw=1.3247+3.3e3*lambda^-2-3.2e7*lambda^-4-2.5e-6*T2+(5-2e-2*T)*4e-5*sal; %from Mcneil (1977)
%note that n_wat should be used instead of n_sw because the salinity
%adjustment below includes this effect
isothermal_compress=(5.062271-0.03179*T+0.000407*T2)*1e-10; %from expt data of Lepple and Millero
(1971)
%multiplication factor incorrectly reported by Buiteveld as 1e-11
%pressure derivative of refractive index
comp1=(-0.000156*lambda+1.5989)*1e-10;
comp2=(1.61857-0.005785*T)*1e-10;
n_pressure_derivative=(comp1*comp2)/1.5014e-010;
%***PURE WATER
%***Einstein-Smoluchowski Eqn
%better to use n_wat and assume Morel's salinity adjustment takes care of
%all of the salinity effect
beta90_wat=2*pi2*k*(273+T)*((lambda*10^-9)^-4)/isothermal_compress*n_wat2*...
n_pressure_derivative2*(6+6*depolar_ratio)/(6-7*depolar_ratio);
beta_wat=beta90_wat*(1+((cos(rad)).^2).*(1-depolar_ratio)/(1+depolar_ratio));
b_wat=8*pi/3*beta90_wat*(2+depolar_ratio)/(1+depolar_ratio);
%***SEAWATER
%30% enhancement for 37ppt seawater approximated from expt and theoretical data
%in Morel(1966,1974);
%probably about ±5% accuracy
betasw=beta_wat*(1+0.3*sal/37);
bsw=b_wat*(1+0.3*sal/37);
%Morel, A. 1966. Etude experimentale de la diffusion de la lumiere par l'eau,
%les solutions de chlorure de sodium, et l'eau de mer optiquement pures.
%J. Chim. Phys. 10:1359-1366.
%Morel, A. 1974. Optical properties of pure water and pure seawater., p.
%1-24. In: N.G.J.E. Steeman-Nielsen [ed.], Optical aspects of oceanography.
%Academic.
