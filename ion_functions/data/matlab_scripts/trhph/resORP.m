function [ORP] = resORP(ORP_V, offset, gain)
% This function computes ORP from raw Resistivity Probe measurements of ORP
% using a Pt-Ag/AgCl electrode pair. Initial inputs are the raw voltages for
% the ORP (ORP_V), offset, and gain

% For testing purposes, assume:
% offset = 2008;
% gain = 4.00;

ORP = ((ORP_V)*1000-offset)/gain;