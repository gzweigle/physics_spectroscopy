% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% This function solves Equation 11.13 from the thesis multiple times, 
% once for each of the original energy levels.
%
% Energies in this function are in units of electronvolts (eV).
%
% by Greg C. Zweigle
%
function [shift_levels_ev Tprime] = energy_shift(M_val, wo_array_ev, ... 
fc_coeff, orig_levels_ev, TiO2_bands)

% Coefficient of exponential approximation term, see Figure 10.12 of thesis.
alpha = 100;

% The semiconductor effect. Coupling parameter K is set at 0.1eV.
% Explanation of this parameter value is given in the JPC paper.
log_curve = 0.1 * log(abs(TiO2_bands(1) - wo_array_ev) ./ ...
                      abs(TiO2_bands(2) - wo_array_ev));

% Loop over all of the M_val initial energy levels.
ref_Ep_ind = 1;

for ref_Ep = orig_levels_ev,

  % Loop over all levels not equal to the initial reference energy level.
  Tprime(ref_Ep_ind) = 0;
  subset_Ep_ind = 1;
  for subset_Ep = orig_levels_ev,

    % Approximation of fc_coeff / (ref_Ep - subset_Ep).
    % This is Equation 11.14 in the thesis.
    % Note that this equation includes "m = n", which is ok because
    % the approximation Equation 10.28 is equal to zero in this case.
    Tprime(ref_Ep_ind) = Tprime(ref_Ep_ind) + ...
      fc_coeff(subset_Ep_ind) * ...
      exp(-alpha * abs(ref_Ep - subset_Ep)) * sign(ref_Ep - subset_Ep);

    subset_Ep_ind = subset_Ep_ind + 1;

  end;

  % Find the intersection of the linear slope curve and the scaled
  % semiconductor effect. This is Equation 11.13 in the thesis.
  exp_curve = log_curve .* ...
    ((wo_array_ev - ref_Ep) * Tprime(ref_Ep_ind) + fc_coeff(ref_Ep_ind));
  min_exp_curve = abs(exp_curve - (wo_array_ev - ref_Ep));

  % This simple algorithm for finding the intersection isn't
  % the most robust approach. But it works for the cases considered.
  min_exp_inds = find(min_exp_curve <= min(min_exp_curve));

  % Find the shifted energy level. In the unlikely event that multiple
  % intersection are found, select the first min_exp_inds because it
  % would be closest to the intersection. This was tested to ensure
  % it never happened.
  shift_levels_ev(ref_Ep_ind) = wo_array_ev(min_exp_inds(1));

  ref_Ep_ind = ref_Ep_ind + 1;

end;
