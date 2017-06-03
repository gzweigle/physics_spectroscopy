% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% Reflects thesis (along with further post-graduation work.... because
% its an interesting problem, and there is room for further ideas).
%
% by Greg C. Zweigle
%
clear;

% Some physical constants
hbar = 6.626e-34 / 2 / pi;  % Js
speed_of_light = 299792458; % m / s
joule_per_ev = 1.602e-19;   % eV / J

% Energy unit conversions from eV to 1/nm and to 1/cm, respectively.
ev_per_nm = 2*pi * hbar * speed_of_light / joule_per_ev * 1e9;  % eV / nm^-1
ev_per_cm = ev_per_nm / 1e7;                                    % eV / cm^-1

% Model parameters:

% Titanium dioxide approximate conduction bands, eV. From Table 1 in JPC paper.
TiO2_bands(1) = 3.6;
TiO2_bands(2) = TiO2_bands(1) + 2.5;

% Franck-Condon normal mode frequencies. From Table 1 in the JPC paper.
omega_cm = [1525 1155 1005];  % Units are cm^-1.

% The free molecule transition frequency, omega sub eg, in units of cm^-1.
% The plots are in units of nanometers, so its helpful to assign the value
% in nanometers and show its conversion to cm^-1. Table 1 in the JPC paper.
weg_cm = 1e7 / 467;  % Units are cm^-1.

% This is the shift in equilibrium position between excited and ground states.
% See Equation 4.38 in thesis. Exact values were unknown, and they were
% selected to match the free spectrum.
delta = [1.30 0.4 0.4];  % Dimensionless.

% Use a Gaussian lineshape to compensate for unmodeled effects such as
% surface states and solvent interaction. See page. 185 of thesis.
% These variables are the variance parameter of the Gaussian. Mean is zero.
gaussian_free_variance = 650;  % Units are cm^-1.
gaussian_bound_variance = 750;  % Units are cm^-1.

% Terms in front of the absorption cross section.  See Equation 4.39 in
% thesis.  This is normalized out later because only interested in
% normalized spectrum for this work.  However, leave in the code as a
% reminder of the full equation.
% ugeo =  Evaluated absorption electronic transition moment at the
%         equilibrium geometry of the ground state (Condon approximation).
%         Ultimately, this value is irrelevant since it is normalized out later.
% gamma = Damping term. This is also used in Equation 4.39 denominator and in
%         that case it is not normalized out.
ugeo = 16.6;
gamma = 200;  % Units are cm^-1.
k_front = 4 * pi^2 * ugeo^2 / 3 / hbar / speed_of_light * gamma / pi;

% This is the number of excited states to consider. The first few dominate,
% and after that they have a rapidly diminishing effect on the spectrum.
% Making this larger significantly slows down the code, so for preliminary
% and fast results, set small. Then, make larger for final results.
total_quantum_numbers = 3;

% Experimental data:

% Experimental data for comparison against results modeled by theory.
% The values in this file were obtained by visual inspection of plot in paper.
load experimental_data.txt;
wo_experiment_nm = experimental_data(:,1);
spectrum_experiment_f = experimental_data(:,2) ./ ...
  max(abs(experimental_data(:,2)));
spectrum_experiment_b = experimental_data(:,3) ./ ...
  max(abs(experimental_data(:,3)));

% Theoretical calculations:

% The energy range of interest.
freq_interval = 10;  % Resolution of energy, in cm^-1.
wo_array_cm = [weg_cm - 5 * omega_cm:freq_interval:weg_cm + 15 * omega_cm(1)];
wo_array_ev = wo_array_cm * ev_per_cm;

% Setup the equation to span the relevant range of energies.
k_front = k_front * wo_array_cm .^2;

% Loops to calculate Equation 4.39 from thesis.
% There are three normal modes modeled.  One for loop for each.
fi = 1;
for v1 = 0:total_quantum_numbers,
  for v2 = 0:total_quantum_numbers,
    for v3 = 0:total_quantum_numbers,

      % Franck-Condon terms.
      delt_prod = prod(((delta.^2) / 2) .^ [v1 v2 v3]);
      fact_prod = prod(factorial([v1 v2 v3]));
      fc_coeff(fi) = delt_prod / fact_prod * exp(-sum(delta.^2)/2);

      % Denominator summation.
      orig_levels_cm(fi) = weg_cm + sum([v1 v2 v3] .* omega_cm);

      % Putting everything together to get Equation 4.39.
      spectrum_calc_f(fi,:) = k_front * fc_coeff(fi) ./ ...
      ((orig_levels_cm(fi) - wo_array_cm).^2 + gamma^2);

      % Increment the index, which reflects present combination
      % of quantum numbers v1, v2, and v3.
      fi = fi + 1;

    end;
  end;
end;

% Store the total number of excited frequencies: this is M in the thesis.
M_val = fi - 1;

% Gaussian model for spectral broadening.
% Using _f to indicate free molecule and _b to indicate bound molecule.
gaussian_f  = exp(-(wo_array_cm - mean(wo_array_cm)) .^2 / 2 / ...
  gaussian_free_variance^2);
gaussian_b = exp(-(wo_array_cm - mean(wo_array_cm)) .^2 / 2 / ...
  gaussian_bound_variance^2);

% Sum all of the individual terms to get a complete spectrum
% then convolve with the Gaussian linewidth.
spectrum_filtered_total_f = sum_and_convolve(spectrum_calc_f, gaussian_f);

% Compute the shifted energies due to molecular coupling to the semiconductor.
% This is where the central model from the thesis lives.
[shift_levels_ev Tprime] = energy_shift(M_val, wo_array_ev, fc_coeff, ...
  orig_levels_cm * ev_per_cm, TiO2_bands);

% Calculate Equation 4.39 again, but this time with the shifted energies
% predicted due to coupling between the molecule and semiconductor.
for fi = 1:length(shift_levels_ev),
  spectrum_calc_b(fi,:) = k_front .* fc_coeff(fi) ./ ...
  ((shift_levels_ev(fi)/ev_per_cm - wo_array_cm).^2 + gamma^2);
end;

% With shifted energies, sum all of the individual terms to get a
% complete spectrum then convolve with the Gaussian linewidth.
spectrum_filtered_total_b = sum_and_convolve(spectrum_calc_b, gaussian_b);

% Everything after this is plotting and displaying results:

% Axis and title for plots.  x-axis units are nanometers.
axis_ranges = [350 550 0 1.1];
molecule = 'CA11/ACOA';

% Plot all of the spectral comparisons.
plot_spectrum(ev_per_nm ./ wo_array_ev, ...
  spectrum_calc_f, spectrum_calc_b, ...
  spectrum_filtered_total_f, spectrum_filtered_total_b, ...
  wo_experiment_nm, spectrum_experiment_f, spectrum_experiment_b, ...
  ev_per_nm ./ TiO2_bands)

% It is interesting to see the impact of the shift on each
% individual Franck-Condon term.
figure(4);
subplot(1,1,1);
plot(ev_per_nm ./ (orig_levels_cm * ev_per_cm), Tprime, 'o', ...
     ev_per_nm ./ (orig_levels_cm * ev_per_cm), fc_coeff, 'x');
grid;
title('Tprime and Franck-Condon vs. original frequency levels');

% Find the peak of the original and coupled spectrum.
orig_peak_ind = find(spectrum_filtered_total_f >= ...
  max(spectrum_filtered_total_f));
shifted_peak_ind = find(spectrum_filtered_total_b >= ...
  max(spectrum_filtered_total_b));

% Convert from eV to nanometers.
orig_peak_nm = ev_per_nm ./ wo_array_ev(orig_peak_ind);
shifted_peak_nm = ev_per_nm ./ wo_array_ev(shifted_peak_ind);

% The difference in the peak locations is the amount of spectral shifting.
diff_peak_nm = orig_peak_nm - shifted_peak_nm;

% Print out the peaks and shift amounts.
dv = fprintf(...
'%s Free peak:%5.1f nm  Bound peak: %5.1f nm  Shift amount: %4.1f nm\n',...
molecule,  orig_peak_nm, shifted_peak_nm, diff_peak_nm);
