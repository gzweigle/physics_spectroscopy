% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% by Greg C. Zweigle
%
function plot_spectrum(wo_array_nm, ...
  spectrum_calc_f, spectrum_calc_b, ...
  spectrum_filtered_total_f, spectrum_filtered_total_b, ...
  wo_experiment_nm, spectrum_experiment_f, spectrum_experiment_b, ...
  TiO2_bands_nm)

% Axis and title for plots.  x-axis units are nanometers.
axis_ranges = [350 550 0 1.1];
molecule = 'CA11/ACOA';

% Plot the free spectrum along with underlying individual terms.
% Also show the semiconductor energies as vertical bars.
figure(1);
subplot(2,1,1);
plot(wo_array_nm, spectrum_filtered_total_f, '--', ...
  wo_array_nm, spectrum_calc_f/max(max(spectrum_calc_f)),...
  [1 1] * TiO2_bands_nm(1), [0 1.1], 'k', ...
  [1 1] * TiO2_bands_nm(2), [0 1.1], 'k', ...
  'LineWidth', 2);
xlabel('nm');
ylabel('normalized intensity');
dv = sprintf('%s free spectrum', molecule);
title(dv);
grid;
axis([axis_ranges]);

% Plot the bound spectrum along with underlying individual terms.
subplot(2,1,2);
plot(wo_array_nm, spectrum_filtered_total_b, '--', ...
  wo_array_nm, spectrum_calc_b/max(max(spectrum_calc_b)), ...
  [1 1] * TiO2_bands_nm(1), [0 1.1], 'k', ...
  [1 1] * TiO2_bands_nm(2), [0 1.1], 'k', ...
  'LineWidth', 2);
xlabel('nm');
ylabel('normalized intensity');
dv = sprintf('%s bound spectrum', molecule);
title(dv);
grid;
axis(axis_ranges);

% Plot the theoretical results for free vs. bound states.
figure(2);
subplot(2,1,1);
plot(wo_array_nm, spectrum_filtered_total_f, '--', ...
  wo_array_nm, spectrum_filtered_total_b, 'LineWidth', 2);
xlabel('nm');
ylabel('normalized intensity');
dv = sprintf('%s theory - free (dashed) vs. bound (solid)', molecule);
title(dv);
grid;
axis(axis_ranges);

% Plot the experimental results for free vs. bound states.
subplot(2,1,2);
plot(wo_experiment_nm, spectrum_experiment_f, '--', ...
  wo_experiment_nm, spectrum_experiment_b, 'LineWidth', 2);
xlabel('nm');
ylabel('normalized intensity');
dv = sprintf('%s experiment - free (dashed) vs. bound (solid)', molecule);
title(dv);
grid;
axis(axis_ranges);

% Direct comparison of theory vs. experimental results for free states.
figure(3);
subplot(2,1,1);
plot(wo_array_nm, spectrum_filtered_total_f, '--', ...
  wo_experiment_nm, spectrum_experiment_f, 'LineWidth', 2);
grid;
xlabel('nm');
ylabel('normalized intensity');
dv = sprintf('%s free - theory (dashed) vs. experiment (solid)', molecule);
title(dv);
axis(axis_ranges);

% Direct comparison of theory vs. experimental results for bound states.
subplot(2,1,2);
plot(wo_array_nm, spectrum_filtered_total_b, '--', ...
  wo_experiment_nm, spectrum_experiment_b, 'LineWidth', 2);
grid;
xlabel('nm');
ylabel('normalized intensity');
dv = sprintf('%s bound - theory (dashed) vs. experiment (solid)', molecule);
title(dv);
axis(axis_ranges);
