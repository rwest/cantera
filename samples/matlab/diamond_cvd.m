%% DIAMOND_CVD - A CVD example simulating growth of a diamond film
% 
% This example computes the growth rate of a diamond film according to a
% simplified version of a particular published growth mechanism (see file
% diamond.yaml for details). Only the surface coverage equations are solved 
% here; the gas composition is fixed. (For an example of coupled gas-phase 
% and surface, see catalytic_combustion.py.) Atomic hydrogen plays an 
% important role in diamond CVD, and this example computes the growth rate 
% and surface coverages as a function of [H] at the surface for 
% fixed temperature and [CH3].
% 
% Requires: cantera >= 2.6.0, pandas >= 0.25.0, matplotlib >= 2.0
% Keywords: surface chemistry, kinetics
%

%% Initialization

help diamond_cvd

clear all
close all
cleanup
clc

t0 = cputime; % record the starting time

%% Operation Parameters

t = 1200.0;                     % surface temperature
p = 20.0 * oneatm;              % pressure

%% Create the gas object
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
% The gas phase will be taken from the definition of phase 'gas' in
% input file 'diamond.yaml'.

gas = Solution('diamond.yaml', 'gas');
gas.TP = {t, p};
x = gas.X;
ih = gas.speciesIndex('H');

%% Create the interface object

