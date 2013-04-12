function [uu,vv,ww] = ins2earth(u,v,w,H,tilt1,tilt2,adcpmode)
%
% USAGE
% [uu,vv,ww] = ins2earth(u,v,w,H,tilt1,tilt2,adcpmode)
%
%
% input
% u,v,w: velocity profile in instrument coordinates
% tilt1: measured pitch, degree
% tilt2: measured roll, degree
% adcpmode: 0 for downward looking ADCP, 1 for upward looking ADCP
%
% output
% uu,vv,ww: velocity profile earth coordinates
%
if adcpmode,
    R = tilt2+180;
else
    R = tilt2;
end

Rrad = deg2rad(R);
Hrad = deg2rad(H);
t1rad = deg2rad(tilt1);
t2rad = deg2rad(tilt2);

P = atan(tan(t1rad).*cos(t2rad));% rad
M1 = [cos(Hrad) sin(Hrad) 0; -sin(Hrad) cos(Hrad) 0;0 0 1];
M2 = [1 0 0;0 cos(P) -sin(P); 0 sin(P) cos(P)];
M3 = [cos(Rrad) 0 sin(Rrad); 0 1 0; -sin(Rrad) 0 cos(Rrad)];

vel = zeros(length(u),3);
for i=1:length(u),
    vel(i,:) = (M1 * M2 * M3 * [u(i); v(i); w(i)])';
end

uu = vel(:,1);
vv = vel(:,2);
ww = vel(:,3);

function radvalue = deg2rad(degvalue)
radvalue = degvalue./180*pi;
