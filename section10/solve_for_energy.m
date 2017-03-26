% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% This code solves Equation 10.11 and generates Figure 10.8 in my thesis.
% See pages 149 - 155.
%
% Search for an Energy that satisfies Equation 10.11
% by locating the zero crossing. Here is the equation:
% coupling * F(Energy) - (Energy - E_phi) = 0
%
% by Greg C. Zweigle
%
function [Energy] = ...
solve_for_energy(coupling, E_phi, E1, E2, E_search_range, direction)

% Equation 10.11.
result = coupling * log(abs(E1 - E_search_range) ./ ...
abs(E2 - E_search_range)) - (E_search_range - E_phi);

% Find the zero crossing.
if direction == 0,
  first_neg = max(find(result >= 0));
else
  first_neg = min(find(result <= 0));
end;

% Return the energy that satisfies Equation 10.11.
Energy = E_search_range(first_neg);
