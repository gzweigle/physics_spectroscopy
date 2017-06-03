% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% This code generates Figure 10.4 in my thesis.
%
% by Greg C. Zweigle
%
clear;

% Normalized range of energy for the horizontal axis.
E = [0:0.001:1];

% Equation 10.7, with data for blue curve in Figure 10.4.
k_squared = 1;
E1 = 0.25;
E2 = 0.75;
z1 = k_squared*log(abs((E1-E)./(E2-E)));

% Equation 10.7, with data for green curve in Figure 10.4.
k_squared = 2;
E1 = 0.25;
E2 = 0.75;
z2 = k_squared*log(abs((E1-E)./(E2-E)));

% Equation 10.7, with data for red curve in Figure 10.4.
k_squared = 1;
E1 = 0.45;
E2 = 0.55;
z3 = k_squared*log(abs((E1-E)./(E2-E)));

% Figure 10.4
plot(E,z1,E,z2,E,z3,'LineWidth',2);
xlabel('normalized energy, E');
ylabel('Z');
grid;
