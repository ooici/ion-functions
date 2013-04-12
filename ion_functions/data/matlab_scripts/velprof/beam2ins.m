function [u,v,w,e] = beam2ins(b1,b2,b3,b4)
%
% input:
% b1, b2, b3, b4: velocity profiles in beam coordinates
%
% output:
% velocity profiles (u,v,w,e) in instrument coordinates
%
%
theta = 20/180*pi;% 20 deg
a = 1./(2.*sin(theta));
b = 1./(4.*cos(theta));
c = 1;% convex (Long Ranger & Quartermaster)
d = a./sqrt(2);

u = c.*a.*(b1-b2);
v = c.*a.*(b4-b3);
w = b.*(b1+b2+b3+b4);
e = d.*(b1+b2-b3-b4);