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

clear;

% Physical constants.
K  = 1.381e-23;                % Boltzmann constant, J/K
h  = 6.262e-34;                % Plank constant, Js
c  = 299792458;                % Speed of light, m/s
Ta = 298.15;                   % Ambient temperature, K
Ts = 5760;                     % Solar temperature, K
q  = 1.602e-19;                % Electric charge, C

% The solid angle of received radiation.
Fs = pi*sin(0.26/360*2*pi)^2;
Fa = pi;

% Step size of integration.
step = 100;

% Operational energy range - 1.1 eV.
Eg = 1.1*q;

% Get the total power so can calculate the efficiency.
Et = K*Ts/step/10:K*Ts/step:100*K*Ts;  % J
b_t = (2*Fs/h^3/c^2)*Et.^2 ./ (exp(Et/K/Ts) - 1);
Pt = sum(b_t.*Et)*K*Ts/step;           % J/s

% Absorbed photon fluxes, start at the bandgap energy and increase to total.
Es = Eg:K*Ts/step:100*K*Ts; b_s = (2*Fs/h^3/c^2)*Es.^2 ./ (exp(Es/K/Ts) - 1);
Ea = Eg:K*Ta/step:100*K*Ta; b_a = (2*Fa/h^3/c^2)*Ea.^2 ./ (exp(Ea/K/Ta) - 1);

% Current due to ambient temperature.
Ja = q * sum(b_a*K*Ts/step);

% The short circuit current is based on the absorbed solar
% radiation and is independent of the applied voltage.
Jsc = q * (sum(b_s)*K*Ts/step);

% Loop over sample short circuit and series resistances.
r_index = 1;
for Rsh = [1e10 100 25] / 100 / 100, % This is Om-m^2
  for Rs = [0 0.2 0.4] / 100 / 100,  % This is Om-m^2

    % Loop over a range of externally applied voltages.
    V_range = [0.6:0.01:0.87];
    v_index = 1;
    for voltage = V_range,

      % Start with Rs = 0 and Rsh = "infinity".
      Jiter(r_index,v_index) = current_calculation(...
      Jsc, Ja, 2*Fa/h^3/c^2, step, Ea, Ta, voltage, 0, 0, 1e10);

      % Iterate to get the converged current value.
      for jind = 1:10,
        Jiter(r_index,v_index) = current_calculation(...
        Jsc, Ja, 2*Fa/h^3/c^2, step, Ea,
        Ta, voltage, Jiter(r_index,v_index), Rs, Rsh);
      end;

      v_index = v_index + 1;

    end;

    % To make the plot nice, remove any currents above the open circuit voltage.
    index = find(Jiter(r_index,:) < 0);
    Jiter(r_index,index) = zeros(1,length(index));

    % Find the peak power. Jiter is in units of A / m^2.
    peak_e(r_index) = max(100*V_range.*Jiter(r_index,:)/Pt);

    % The efficiency is the ratio of actual over maximum possible power.
    efficiency(r_index,:) = 100*V_range.*Jiter(r_index,:)/Pt;

    % Compute the fill factor.
    fill(r_index) = peak_e(r_index) * Pt / (Jiter(r_index,1) * ...
    V_range(min(find(Jiter(r_index,:)==0)))) / 100;

    r_index = r_index + 1;

  end;

end;

% Divide by 10 to convert A/m^2 to mA/cm^2.
figure(1);
plot(V_range,Jiter'/10);
grid;
xlabel('V: volts');
ylabel('J: mA/cm^2');

figure(2);
plot(V_range,efficiency);
grid;
xlabel('V: volts');
ylabel('Percent Efficiency');
