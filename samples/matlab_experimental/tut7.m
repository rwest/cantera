%  Tutorial 7: Thermodynamic Properties
%
% Keywords: tutorial, thermodynamics

help tut7

LoadCantera;
clear all
close all
cleanup

% A variety of thermodynamic property methods are provided.
gas = Air
gas.TP = {800, oneatm}

% temperature, pressure, density
Temperature = gas.T
Pressure = gas.P
rho = gas.D
n = gas.molarDensity

% species non-dimensional properties
hrt = gas.enthalpies_RT            % vector of h_k/RT

% Unlike the legacy interface, each thermodynamic property is tied to the
% `basis`('molar' or 'mass') instead of being two separate properties. 

% mixture properties per mole
gas.basis = 'molar';
hmole = gas.H
umole = gas.U
smole = gas.S
gmole = gas.G

% mixture properties per unit mass
gas.basis = 'mass';
hmass = gas.H
umass = gas.U
smass = gas.S
gmass = gas.G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
cleanup
