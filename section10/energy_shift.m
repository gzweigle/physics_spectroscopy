% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% This code solves Equation 10.11 and generates Figure 10.8 in my thesis.
% See pages 149 - 155.
%
% by Greg C. Zweigle
%
clear;

% Ground referenced lower semiconductor band energy.
E1 = 0.5;

% Ground referenced upper semiconductor band energy.
E2 = 1.0;

% E_phi is the energy of the perturbed excited molecular state, Equation 7.13.
% These equations define the step size for solving Equation 10.11 for E_phi.
E_phi_step = 0.01;
search_step = 1e-5;

% Break energy into segments to avoid solution discontinuities.
E_phi_range1 = [0.1:E_phi_step:E1-E_phi_step];
E_phi_range2 = [E1+E_phi_step:E_phi_step:(E1+E2)/2-E_phi_step];
E_phi_range3 = [(E1+E2)/2+E_phi_step:E_phi_step:E2-E_phi_step];
E_phi_range4 = [E2+E_phi_step:E_phi_step:1.4];
E_phi_range_all = [E_phi_range1 E_phi_range2 E_phi_range3 E_phi_range4];

% To show the trend, solve Equation 10.11 for three values of coupling.
coupling_range = [0.002 0.001 0.0005];
coupling_index = 1;
for coupling = coupling_range,

  % Solve for each value of E_phi (horizontal axis of plot).
  E_phi_index = 1;
  for E_phi = E_phi_range_all,

    % Handle two cases:
    % - Search starting at the mid-point and moving to higher energies.
    % - Search starting a low energies and moving toward the mid-point energy.
    if E_phi > (E1+E2)/2,
      [Energy] = solve_for_energy(coupling, E_phi, E1, E2,...
      [E_phi+search_step:search_step:10*E_phi], 1);
    else
      [Energy] = solve_for_energy(coupling, E_phi, E1, E2,...
      [-10*E_phi:search_step:E_phi+search_step], 0);
    end;

    % Save the results for plotting.
    E_save(coupling_index,E_phi_index) = Energy;

    E_phi_index = E_phi_index + 1;

  end;
 
  coupling_index = coupling_index + 1;

end;

% Plot the curves along with two vertical lines as references.
plot(E_phi_range_all, E_phi_range_all - E_save, ...
[.5 .5],[-0.015 0.015],'k',[1 1],[-0.015 0.015],'k','LineWidth',2);
xlabel('E\_phi');
ylabel('E\_phi - E');
grid;
