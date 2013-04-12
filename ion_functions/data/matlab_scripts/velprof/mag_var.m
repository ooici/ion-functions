function [U,V] = mag_var(theta,u,v)
% magnetic variation correction
%
% input:
% theta (degree)
% u,v: horizontal velocity profile
%
%
% output:
% U,V
%
theta_rad = deg2rad(theta);

M = [cos(theta_rad) sin(theta_rad);-sin(theta_rad) cos(theta_rad)];

tmp = M * [u,v]';
U = tmp(1);
V = tmp(2);