clc; clear all; close all;

% Parameters
frequency_Hz = 165e6;  
speedOfLight_m_s = 3e8;
wavelength_m = speedOfLight_m_s / frequency_Hz;
numDirectors = 4;

% Antenna Parameters
initialAntenna = wavelength_m * [0.5, 0.47, 0.406, 0.406, 0.406, 0.406, 0.2, 0.34, 0.34, 0.34, 0.34];
                                 

% Driven Element Configuration
d = dipoleFolded;
d.Length = initialAntenna(2);
d.Width = cylinder2strip(0.003 * wavelength_m);
d.Spacing = d.Length / 60;
d.Tilt = 0;
d.TiltAxis = "X";

% Yagi-Uda Antenna Design
yagidesign = yagiUda;
yagidesign.Exciter = d;
yagidesign.NumDirectors = numDirectors;
yagidesign.ReflectorLength = initialAntenna(1);
yagidesign.DirectorLength = initialAntenna(3:6);
yagidesign.ReflectorSpacing = initialAntenna(7);
yagidesign.DirectorSpacing = initialAntenna(8:11);

% Gain Calculation
gain = max(pattern(yagidesign, frequency_Hz, 0, 0:1:360));

% Side Lobe Level (SLL) Calculation
lobeInfo = findLobes(polarpattern(pattern(yagidesign, frequency_Hz, 0, 0:1:360)));
sll = -lobeInfo.SLL;

% Display Results
disp('Initial Antenna Configuration (meters):'); disp(initialAntenna);
disp(['Antenna Gain (dBi): ', num2str(gain)]);
disp(['Initial Antenna SLL (dB): ', num2str(sll)]);

% Visualization

figure; pattern(yagidesign, frequency_Hz); title('3D Radiation Pattern');

figure; pattern(yagidesign, frequency_Hz, 0, 0:1:360); title('E-plane Pattern (Azimuth = 0°)');

figure; pattern(yagidesign, frequency_Hz, 90, 0:1:360); title('H-plane Pattern (Elevation = 90°)');

figure; yagidesign.Tilt = 90; yagidesign.TiltAxis = [0 1 0]; show(yagidesign);
title('Initial Yagi-Uda Antenna Model');

% 2D Polar Plot
theta = -180:1:180;
gain_dBi = pattern(yagidesign, frequency_Hz, 0, theta);

figure;
plot(theta, gain_dBi, 'r', 'LineWidth', 2);
grid on;

xlabel('Theta (degree)'); ylabel('Gain (dBi)');
title('2D Radiation Pattern');
xlim([-180 180]); ylim([-40 20]); xticks(-180:45:180);