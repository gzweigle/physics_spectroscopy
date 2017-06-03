% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% This code generates Figure 10.12 in my thesis.
% See page 168.
%
% by Greg C. Zweigle
%
clear

E = [0.4:0.001:.6 - 0.001];
a = 200;
k = 1000;
plot(E, (1/k)./(0.5-E), E, exp(-a*abs(0.5-E)).*sign(0.5-E), 'LineWidth', 2);
axis([0.4 0.6 -1 1]);
xlabel('Em');
grid;
