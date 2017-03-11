% Modeling the Spectroscopy of a Light Collecting Molecule Coupled to a
% Nanocrystalline Semiconductor
%
% Master of Science in Chemistry
% Washington State University
%
% This code generates Figure 6.2, 6.3, and 6.5 in my thesis.
% See page 89 - 93 for details.
%
% See Figure 6.3 for details of the geometry between two
% In the Matlab code, capital A and B represent alpha and beta, respectively.
%
% Greg C. Zweigle
%

clear;

% Convert degrees to radians
d2r = 2.0 * pi / 360.0;

% Electron distance from origin, normalized to nucleus distance.
R_electron = 0.5;

% Unit vector connecting the center point of the molecules.
% First term is along x-axis, second term is along y-axis.
Rhat = [1 0];

% Distance between the center point of the two molecules.
R_range = [4:25];

% Angle between electron centerline of electron b / nucleus B.
a_range = [0:5:85 95:5:180];

% Rotate the dipole connecting electron b / nucleus B.
alpha_index = 1; for alpha = a_range,

  % Vary the distance between the two molecules.
  R_index = 1; for R_molecules = R_range,

    % Rotate the dipole connecting electron a / nucleus A.
    T_index = 1; for theta = [0 45 90 180],

      % Vector between the molecules.
      R_v = R_molecules * Rhat;
      
      % Position vector of electron a.
      ra_v = R_electron * [cos(theta*d2r) sin(theta*d2r)];

      % Position vector of nucleus A.
      rA_v = -ra_v;

      % Position vector of electron b.
      rb_v = R_v + R_electron * [-cos((90-alpha)*d2r) sin((90-alpha)*d2r)];

      % Position vector of nucleus B.
      rB_v = R_v + R_electron * [cos((90-alpha)*d2r) -sin((90-alpha)*d2r)];

      % Vector pointing from electron a to nucleus A.
      raA_v = ra_v - rA_v;

      % Vector pointing from electron b to nucleus B.
      rbB_v = rb_v - rB_v;

      % Dipole-dipole potential.
      equation_6_29 = (raA_v * rbB_v' - ...
      3 * (raA_v * Rhat')*(rbB_v * Rhat')) / R_molecules^3;

      % For one particular distance, save numerator of Equation 6.18
      % so can plot as Figure 6.2.
      if R_molecules == 10,

        lone_dipole_plot(alpha_index,T_index) = equation_6_29;

        G_orient(alpha_index,T_index) = ...
	  2 * cos(theta*d2r) .* sin(alpha*d2r) + ...
              sin(theta*d2r) .* cos(alpha*d2r);

      end;

      % For one particular molecule aA orientation,
      % save dipole and Hamiltonian potentials.
      if theta == 90,
	
	% Equation 6.29 in Thesis.
        dipole(alpha_index,R_index) = equation_6_29;

	% Equation 6.28 in Thesis.
        hamilt(alpha_index,R_index) = ...
	  -1/sqrt((rb_v-rA_v)*(rb_v-rA_v)') + ...
          -1/sqrt((ra_v-rB_v)*(ra_v-rB_v)') + ...
           1/sqrt((rb_v-ra_v)*(rb_v-ra_v)') + ...
           1/sqrt((rB_v-rA_v)*(rB_v-rA_v)');

      end;

      T_index = T_index + 1;

    end;

  R_index = R_index + 1;

  end;

alpha_index = alpha_index + 1;

end;

% Calculate the ratio, Equation 6.30.
E = (dipole./hamilt - 1) * 100;

% Create Figure 6.4.
figure(1);
subplot(1,1,1);
meshc(R_range,a_range,E);
xlabel('R');
ylabel('alpha');
zlabel('Error');

% Plot Equation 6.28 and Equation 6.29 individually.  
figure(2);
subplot(2,1,1);
plot(a_range,hamilt,'LineWidth',2);
xlabel('alpha, degrees');
ylabel('Vh');
grid;
subplot(2,1,2);
plot(a_range,dipole,'LineWidth',2);
xlabel('alpha, degrees');
ylabel('Vdd');
grid;

% Create Figure 6.5.
figure(3);
subplot(1,1,1);
ind = find(R_range == 10);
plot(R_range(ind:length(R_range)),abs(E(1,ind:length(R_range))),...
'LineWidth',2);
xlabel('R');
ylabel('Error');
grid;

% Create Figure 6.2 two different ways. They match.
figure(4);
subplot(1,1,1);
plot(a_range,G_orient,'LineWidth',2);
grid;
xlabel('alpha, degrees');
ylabel('Orientation portion of G');
figure(5);
subplot(1,1,1);
plot(a_range,lone_dipole_plot*1000,'LineWidth',2);
grid;
xlabel('Check for correct orientation');
ylabel('alpha, degrees');
