% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% This code generates Figure 1.5 and Figure 1.6 in my thesis.
% See pg. 10-12 for details.
%
% Greg C. Zweigle
%
function [J] = current_calculation(Jsc,Ja,coeff,step,Ea,Ta,V,Jtry,Rs,Rsh)

K  = 1.381e-23;  % Boltzmann constant, J/K
q  = 1.602e-19;  % Electric charge, C

% Voltage dependent emitted photon flux.
b_u = coeff*Ea.^2 ./ (exp((Ea-q*(V+Jtry*Rs))/K/Ta) - 1);

% Current is sum of absorbed minus emitted.
Jdiode = q * (sum(b_u)*K*Ta/step) - Ja;

% Check the iterated current.
J = (Jsc - Jdiode - V / Rsh) / (1 + Rs/Rsh);
