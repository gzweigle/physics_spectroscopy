% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% Summation and convolution happens twice, once for the decoupled energy
% levels and once for the coupled energy levels.
% So, break out into its own function.
%
% by Greg C. Zweigle
%
function spectrum_filtered_total = sum_and_convolve(spectrum, gaussian)

% Create the complete spectrum by adding up all of the individual contributors.
spectrum_summation = sum(spectrum);

% Convolve with the Gaussian
conv_spectrum = conv(gaussian, spectrum_summation);

% Shift by half of the length of the % Gaussian in order to recenter.
centered_spectrum = conv_spectrum(round(length(gaussian) / 2) : ...
  length(conv_spectrum));
spectrum_filtered_total = centered_spectrum(1:length(spectrum_summation));

% Normalize again.
spectrum_filtered_total = spectrum_filtered_total ./ ...
  max(spectrum_filtered_total);
