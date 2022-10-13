% Tutorial 5:   Reaction information and rates
%
%    Topics:
%       - stoichiometric coefficients
%       - reaction rates of progress
%       - species production rates
%       - reaction equations
%       - equilibrium constants
%       - rate multipliers
%
% Keywords: tutorial, kinetics

help tut5
LoadCantera;
clear all
close all

g = Solution('gri30.yaml', 'gri30', 'None');
g.TPX = {1500, oneatm, ones(g.nSpecies, 1)};

% Methods are provided that compute many quantities of interest for
% kinetics. Some of these are:


% 1) Stoichiometric coefficients

nu_r   = g.stoichReactant    % reactant stoichiometric coefficient mstix
nu_p   = g.stoichProduct     % product stoichiometric coefficient mstix
nu_net = g.stoichNet         % net (product - reactant) stoichiometric
                             % coefficient mstix

% For any of these, the (k,i) matrix element is the stoichiometric
% coefficient of species k in reaction i. Since these coefficient
% matrices are very sparse, they are implemented as MATLAB sparse
% matrices.


% 2) Reaction rates of progress

% Methods ropForward, ropReverse, and ropNet return column vectors containing
% the forward, reverse, and net (forward - reverse) rates of
% progress, respectively, for all reactions.

qf = g.ropForward;
qr = g.ropReverse;
qn = g.ropNet;
rop = [qf, qr, qn]

% This plots the rates of progress
figure(1);
bar(rop);
legend('forward','reverse','net');

% 3) Species production rates

% Methods creationRates, destructionRates, and netProdRates return
% column vectors containing the creation, destruction, and net
% production (creation - destruction) rates, respectively, for all species.

cdot = g.creationRates;
ddot = g.destructionRates;
wdot = g.netProdRates;
rates = [cdot, ddot, wdot]

% This plots the production rates
figure(2);
bar(rates);
legend('creation','destruction','net');

% For comparison, the production rates may also be computed
% directly from the rates of progress and stoichiometric
% coefficients.

cdot2 = nu_p*qf' + nu_r*qr';
creation = [cdot; cdot2'; cdot - cdot2']

ddot2 = nu_r*qf' + nu_p*qr';
destruction = [ddot; ddot2'; ddot - ddot2']

wdot2 = nu_net * qn';
net = [wdot; wdot2'; wdot - wdot2']

% 4) Reaction equations

e8    = g.reactionEqn(8)             % equation for reaction 8
e1_10 = g.reactionEqns(1:10)          % equation for rxns 1 - 10
eqs   = g.reactionEqns                % all equations

% 5) Equilibrium constants

% The equilibrium constants are computed in concentration units,
% with concentrations in kmol/m^3.

kc = g.equilibriumConstants;
for i = 1:g.nReactions
   fprintf('%50s  %13.5g', eqs{i}, kc(i))
end

% 6) Multipliers

% For each reaction, a multiplier may be specified that is applied
% to the forward rate coefficient. By default, the multiplier is
% 1.0 for all reactions.

for i = 1:g.nReactions
   g.setMultiplier(i, 2*i);
   m = g.multiplier(i);
end

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
