clc; clear all; close all;

% Parameters
frequency_Hz = 165e6;  
speedOfLight_m_s = 3e8;
wavelength_m = speedOfLight_m_s / frequency_Hz;
numDirectors = 4;

% SA Antenna Parameters
SA_SLL_Cons_Ant = [0.9091, 0.8462, 0.8000, 0.7906, 0.7709, 0.7467, 0.2906, 0.5120, 0.6564, 0.6577, 0.7256];

% SA Driven Element Configuration
d = dipoleFolded;
d.Length = SA_SLL_Cons_Ant(2);
d.Width = cylinder2strip(0.003 * wavelength_m);
d.Spacing = d.Length / 60;

% GA Yagi-Uda Antenna Design
yagidesign = yagiUda;
yagidesign.Exciter = d;
yagidesign.NumDirectors = numDirectors;
yagidesign.ReflectorLength = SA_SLL_Cons_Ant(1);
yagidesign.DirectorLength = SA_SLL_Cons_Ant(3:6);
yagidesign.ReflectorSpacing = SA_SLL_Cons_Ant(7);
yagidesign.DirectorSpacing = SA_SLL_Cons_Ant(8:11);

% GA Antenna Parameters
GA_SLL_Cons_Ant = [0.9091, 0.8880, 0.8000, 0.7794, 0.7588, 0.7858, 0.2921, 0.4825, 0.7380, 0.7280, 0.6578];
                                 
% GA Design Driven Element Configuration
d = dipoleFolded;
d.Length = GA_SLL_Cons_Ant(2);
d.Width = cylinder2strip(0.003 * wavelength_m);
d.Spacing = d.Length / 60;

% GA Yagi-Uda Antenna Design
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
disp('SA-Optimised Antenna Configuration (meters):'); disp(SA_SLL_Cons_Ant);
disp(['SA-Optimised Antenna Gain (dBi): ', num2str(gain)]);
disp(['SA-Optimised Antenna SLL (dB): ', num2str(sll)]);
disp(' ');
disp('GA-Optimised Antenna Configuration (meters):'); disp(GA_SLL_Cons_Ant);
disp(['GA-Optimised Antenna Gain (dBi): ', num2str(gain_2)]);
disp(['GA-Optimised Antenna SLL (dB): ', num2str(sll_2)]);

% Visualizations

% Polar Radiation Plot
figure; pattern(yagidesign, frequency_Hz); title('Polar Radiation Pattern');
figure; pattern(yagidesign_2, frequency_Hz);title('Polar Radiation Pattern');

% E-plane Radiation Pattern (Azimuth = 0°)
theta_deg = 0:1:360;

% E-plane gain patterns for both antennas
Eplane_gain_unconstrained = pattern(yagidesign, frequency_Hz, 0, theta_deg);
Eplane_gain_constrained = pattern(yagidesign_2, frequency_Hz, 0, theta_deg);

% Convertion from degrees to radians for polarplot
theta_rad = deg2rad(theta_deg);

% E-plane Plots
figure;
polarplot(theta_rad, Eplane_gain_unconstrained, ':b', 'LineWidth', 2); hold on;
polarplot(theta_rad, Eplane_gain_constrained, '-r', 'LineWidth', 2);
title('E-plane Polar Radiation Pattern (Azimuth = 0°)');
legend('SA-Optimised Design', 'GA-Optimised Design', 'Location', 'northoutside');
rlim([-40 20]); 
thetaticks(0:30:360);
hold off;

% H-plane gain patterns for both antennas
Hplane_gain_unconstrained = pattern(yagidesign, frequency_Hz, 90, theta_deg);
Hplane_gain_constrained = pattern(yagidesign_2, frequency_Hz, 90, theta_deg);

% Polar plot Creation
figure;
polarplot(theta_rad, Hplane_gain_unconstrained, ':b', 'LineWidth', 2); hold on;
polarplot(theta_rad, Hplane_gain_constrained, '-r', 'LineWidth', 2);
title('H-plane Polar Radiation Pattern (Azimuth = 90°)');
legend('SA-Optimised Design', 'GA-Optimised Design', 'Location', 'northoutside');
rlim([-40 20]); 
thetaticks(0:30:360);
hold off;

% Display Antenna Models
figure; yagidesign.Tilt = 90; yagidesign.TiltAxis = [0 1 0]; show(yagidesign);
title('SA-Optimised Yagi-Uda Antenna Model');

figure; yagidesign_2.Tilt = 90; yagidesign_2.TiltAxis = [0 1 0]; show(yagidesign_2);
title('GA-Optimised Yagi-Uda Antenna Model');

% 2D Polar Plots
theta = -180:1:180;
gain_dBi = pattern(yagidesign, frequency_Hz, 0, theta);
constrainedGain_dBi = pattern(yagidesign_2, frequency_Hz, 0, theta);

figure;
hold on; 
plot(theta, gain_dBi, ':b', 'LineWidth', 2); 
plot(theta, constrainedGain_dBi, '-r', 'LineWidth', 2); 
grid on;
xlabel('Theta (degree)');
ylabel('Gain (dBi)');
title('Cartesian Radiation Patterns of SA and GA Antennas');
legend('SA-Optimised Design', 'GA-Optimised Design');
xlim([-180 180]); ylim([-30 20]); xticks(-180:45:180);
hold off; 
