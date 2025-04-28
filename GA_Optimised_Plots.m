clc; clear all; close all;

% Parameters
frequency_Hz = 165e6;  
speedOfLight_m_s = 3e8;
wavelength_m = speedOfLight_m_s / frequency_Hz;
numDirectors = 4;

% Unconstrained Antenna Parameters
GAAntenna1 = [0.9091, 0.8909, 0.8000, 0.7884, 0.7897, 0.7871, 0.1818, 0.5777, 0.6729, 0.6628, 0.7380];                              

% Unconstrained Driven Element Configuration
d = dipoleFolded;
d.Length = GAAntenna1(2);
d.Width = cylinder2strip(0.003 * wavelength_m);
d.Spacing = d.Length / 60;

% Unconstrained Yagi-Uda Antenna Design
yagidesign = yagiUda;
yagidesign.Exciter = d;
yagidesign.NumDirectors = numDirectors;
yagidesign.ReflectorLength = GAAntenna1(1);
yagidesign.DirectorLength = GAAntenna1(3:6);
yagidesign.ReflectorSpacing = GAAntenna1(7);
yagidesign.DirectorSpacing = GAAntenna1(8:11);

% Constrained Antenna Parameters
GA_SLL_Cons_Ant = [0.9091, 0.8880, 0.8000, 0.7794, 0.7588, 0.7858, 0.2921, 0.4825, 0.7380, 0.7280, 0.6578];
                                 
% Constrained Driven Element Configuration
d = dipoleFolded;
d.Length = GA_SLL_Cons_Ant(2);
d.Width = cylinder2strip(0.003 * wavelength_m);
d.Spacing = d.Length / 60;

% Constrained Yagi-Uda Antenna Design
yagidesign_2 = yagiUda;
yagidesign_2.Exciter = d;
yagidesign_2.NumDirectors = numDirectors;
yagidesign_2.ReflectorLength = GA_SLL_Cons_Ant(1);
yagidesign_2.DirectorLength = GA_SLL_Cons_Ant(3:6);
yagidesign_2.ReflectorSpacing = GA_SLL_Cons_Ant(7);
yagidesign_2.DirectorSpacing = GA_SLL_Cons_Ant(8:11);

% Gain Calculations
gain = max(pattern(yagidesign, frequency_Hz, 0, 0:1:360));
gain_2 = max(pattern(yagidesign_2, frequency_Hz, 0, 0:1:360));

% Side Lobe Level (SLL) Calculations
lobeInfo = findLobes(polarpattern(pattern(yagidesign, frequency_Hz, 0, 0:1:360)));
sll = -lobeInfo.SLL;
lobeInfo_2 = findLobes(polarpattern(pattern(yagidesign_2, frequency_Hz, 0, 0:1:360)));
sll_2 = -lobeInfo_2.SLL;

% Display Results
disp('Unconstrained Antenna Configuration (meters):'); disp(GAAntenna1);
disp(['Unconstrained Antenna Gain (dBi): ', num2str(gain)]);
disp(['Unconstrained Antenna SLL (dB): ', num2str(sll)]);
disp(' ');
disp('Constrained Antenna Configuration (meters):'); disp(GA_SLL_Cons_Ant);
disp(['Constrained Antenna Gain (dBi): ', num2str(gain_2)]);
disp(['Constrained Antenna SLL (dB): ', num2str(sll_2)]);

% Visualizations

% Polar Radiation Plot
figure; pattern(yagidesign, frequency_Hz); title('Radiation Pattern');
figure; pattern(yagidesign_2, frequency_Hz);title('Radiation Pattern');

% E-plane Radiation Pattern (Azimuth = 0°)
theta_deg = 0:1:360;

% E-plane gain patterns for both antennas
Eplane_gain_unconstrained = pattern(yagidesign, frequency_Hz, 0, theta_deg);
Eplane_gain_constrained = pattern(yagidesign_2, frequency_Hz, 0, theta_deg);

% Convertion from degrees to radians for polarplot
theta_rad = deg2rad(theta_deg);

% Polar plot Creation
figure;
polarplot(theta_rad, Eplane_gain_unconstrained, ':b', 'LineWidth', 2); hold on;
polarplot(theta_rad, Eplane_gain_constrained, '-r', 'LineWidth', 2);
title('E-Plane Polar Radiation Pattern (Azimuth = 0°)');
legend('Unconstrained Antenna Design', 'Constrained Antenna Design', 'Location', 'northoutside');
rlim([-30 20]); 
thetaticks(0:30:360);
hold off;

% H-plane gain patterns for both antennas
Hplane_gain_unconstrained = pattern(yagidesign, frequency_Hz, 90, theta_deg);
Hplane_gain_constrained = pattern(yagidesign_2, frequency_Hz, 90, theta_deg);

% Polar plot Creation
figure;
polarplot(theta_rad, Hplane_gain_unconstrained, ':b', 'LineWidth', 2); hold on;
polarplot(theta_rad, Hplane_gain_constrained, '-r', 'LineWidth', 2);
title('H-Plane Polar Radiation Pattern (Azimuth = 90°)');
legend('Unconstrained Antenna Design', 'Constrained Antenna Design', 'Location', 'northoutside');
rlim([-30 20]); 
thetaticks(0:30:360);
hold off;

% Display Antenna Models
figure; yagidesign.Tilt = 90; yagidesign.TiltAxis = [0 1 0]; show(yagidesign);
title('Unconstrained Yagi-Uda Antenna Model');

figure; yagidesign_2.Tilt = 90; yagidesign_2.TiltAxis = [0 1 0]; show(yagidesign_2);
title('Constrained Yagi-Uda Antenna Model');

% 2D Polar Plots
theta = -180:1:180;
gain_dBi = pattern(yagidesign, frequency_Hz, 0, theta);
constrainedGain_dBi = pattern(yagidesign_2, frequency_Hz, 0, theta);

figure;
hold on; 
plot(theta, gain_dBi, ':b', 'LineWidth', 2); % Unconstrained Antenna (Blue)
plot(theta, constrainedGain_dBi, '-r', 'LineWidth', 2); % Constrained Antenna (Red)
grid on;
xlabel('Theta (degree)');
ylabel('Gain (dBi)');
title('Radiation Pattern of GA Antennas');
legend('Unconstrained Antenna Design', 'Constrained Antenna Design');
xlim([-180 180]); ylim([-30 20]); xticks(-180:45:180);
hold off; 
