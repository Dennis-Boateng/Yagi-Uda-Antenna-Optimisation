% Simulated Annealing for Yagi-Uda Antenna Optimisation
clc; clear all;

% Simulation Parameters
frequency_Hz = 165e6;  % Operating frequency
maxSLL = -15; % Max allowable side lobe level (SLL)
initialTemp = 100;     % Initial temperature
coolingRate = 0.999;    % Cooling factor per iteration
minTemp = 1;        % Stop when temperature is low

% Define antenna parameter ranges (in meters)
speedOfLight_m_s = 3e8;
wavelength_m = speedOfLight_m_s / frequency_Hz;
reflectorRange = [0.5, 0.55] * wavelength_m;
drivenRange = [0.45, 0.49] * wavelength_m;
directorRange = [0.35, 0.44] * wavelength_m;
spacingRange = [0.1, 0.45] * wavelength_m;
elementRadius = 0.003 * wavelength_m;

numElements = 6;  % Reflector + Exciter + 4 Directors
numDirectors = 4;

tic;  % Start timing

% Initialize antenna parameters with progressive director decay
currentSolution = [
    rand() * (reflectorRange(2) - reflectorRange(1)) + reflectorRange(1), ...
    rand() * (drivenRange(2) - drivenRange(1)) + drivenRange(1), ... 
    rand() * (directorRange(2) - directorRange(1)) + directorRange(1), ... 
];

% Apply decay logic for remaining directors
for j = 4:6
    decayFactor = 0.95 + rand() * (0.98 - 0.95);  % Slight decay variation
    value = decayFactor * currentSolution(1, j - 1);
    currentSolution(1, j) = max(directorRange(1), value);  % Clamp to min
end

for j = 7:11
    currentSolution(1, j) = spacingRange(1) + (rand() * (spacingRange(2) - spacingRange(1)));
end

[currentFitness, gain, sll] = evaluateFitnessAntennaToolbox(currentSolution, frequency_Hz, maxSLL, numDirectors, elementRadius);

bestSolution = currentSolution;
bestFitness = currentFitness;
temp = initialTemp;
iteration = 1;

fprintf('Iteration %d - Temp: %.4f - Best Fitness: %.2f\n', iteration, temp, bestFitness);

goalFitness = 16;  % Define target fitness
fitnessHistory = [];
sllHistory = [];
gainHistory = [];
sllHistory(1) = sll;
gainHistory(1) = gain;
fitnessHistory(1) = currentFitness;

% Start Simulated Annealing Process
while temp > minTemp 
    
    % Generate a new solution with a small perturbation
    k = 0.03;  % Scaling factor for perturbation
    alpha = 0.7;  % Temperature sensitivity exponent

    perturbedSolution = currentSolution + k * (temp^alpha) * (rand(size(currentSolution)) - 0.5);

    % Ensure Reflector is longer than Driven Element
    reflectorLength = max(reflectorRange(1), min(reflectorRange(2), perturbedSolution(1)));
    drivenLength = max(drivenRange(1), min(drivenRange(2), perturbedSolution(2)));
    

    if reflectorLength <= drivenLength
        reflectorLength = drivenLength * (1.05 + rand() * 0.03); % Ensure reflector > driven with small buffer
    end

    newSolution(1) = reflectorLength; % Assign corrected reflector length
    newSolution(2) = drivenLength;    % Assign corrected driven length

    % Ensure the first director does not get stuck at minimum value
    bufferFactor = 1.05;  % Buffer to prevent stagnation
    firstDirector = max(directorRange(1) * bufferFactor, min(directorRange(2), perturbedSolution(3)));
    newSolution(3) = firstDirector;  % First director, allowing freedom to move

    for j = 4:6  % Remaining directors
        decayFactor = 0.95 + (0.98 - 0.95) * k * (temp^alpha); % Adaptive decay with temp
        newSolution(j) = max(directorRange(1), decayFactor * newSolution(j-1)); % Maintain decreasing trend
    end

    % Ensure new values remain within valid ranges
    rangeMin = [reflectorRange(1), drivenRange(1), repmat(directorRange(1), 1, 4), repmat(spacingRange(1), 1, 5)];
    rangeMax = [reflectorRange(2), drivenRange(2), repmat(directorRange(2), 1, 4), repmat(spacingRange(2), 1, 5)];

    newSolution = min(max(perturbedSolution, rangeMin), rangeMax);
    
    % Evaluate new solution
    [newFitness, gain, sll] = evaluateFitnessAntennaToolbox(newSolution, frequency_Hz, maxSLL, numDirectors, elementRadius);
    
    % Store values for plotting
    sllHistory(end+1) = sll;
    gainHistory(end+1) = gain;

    % Acceptance criterion (Sigmoid Function)
    delta = newFitness - currentFitness;
    acceptProb = 1 / (1 + exp(-delta / temp));

    if newFitness > currentFitness || rand() < acceptProb
        currentSolution = newSolution;
        currentFitness = newFitness;
    end
    
    % Update best solution if needed
    if currentFitness > bestFitness
        bestFitness = currentFitness;
        bestSolution = currentSolution;
    end

    fitnessHistory(end+1) = bestFitness;
    
    % Cool down temperature
    temp = temp * coolingRate;
    iteration = iteration + 1;
    
    % Display progress
    fprintf('Iteration %d - Temp: %.4f - Best Fitness: %.2f\n', iteration, temp, bestFitness);
     
    % Check if goal fitness is met
    if newFitness >= goalFitness
        break;  % Exit loop when goal is reached
    end
end

elapsedTime = toc;  % Stop timing

% Display final optimized antenna configuration
disp('Optimized Antenna Configuration (meters):');
disp(bestSolution);
disp(['Best Fitness: ', num2str(bestFitness)]);
disp(['Best SLL (dB): ', num2str(sll)]);
disp(['Time to reach goal: ', num2str(elapsedTime), ' seconds']);

d = dipoleFolded;
d.Length = bestSolution(2);
d.Width = cylinder2strip(elementRadius);
d.Spacing = d.Length/60;

yagidesign = yagiUda;
yagidesign.Exciter = d;
yagidesign.NumDirectors = numDirectors;
yagidesign.ReflectorLength = bestSolution(1);
yagidesign.DirectorLength = bestSolution(3:6);
yagidesign.ReflectorSpacing = bestSolution(7);
yagidesign.DirectorSpacing = bestSolution(8:11);

pattern(yagidesign,frequency_Hz);

figure;
pattern(yagidesign,frequency_Hz,0,0:1:360);

figure;
pattern(yagidesign,frequency_Hz,90,0:1:360);

theta = -180:1:180;
gain_dBi = pattern(yagidesign, frequency_Hz, 0, theta); 

% Plot
figure;
plot(theta, gain_dBi, 'r', 'LineWidth', 2);  % Red line
grid on;

% Axis labels and title
xlabel('Theta (degree)');
ylabel('Gain (dBi)');
title('Radiation Pattern');

% Set axis limits for clarity (adjust as needed)
xlim([-180 180]);
ylim([-40 20]);
xticks(-180:45:180);

% Plot Fitness vs Iterations
figure;
plot(1:length(fitnessHistory), fitnessHistory, '-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Fitness Score');
title('Fitness Convergence Performance');
grid on;

figure;
plot(1:length(gainHistory), gainHistory, '-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Gain (dBi)');
title('Gain Convergence Performance');
grid on;

% Plot SLL vs Iterations
figure;
plot(1:length(sllHistory), sllHistory, '-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Side Lobe Level (SLL) in dB');
title('SLL vs. Iterations');
grid on;

figure;
scatter(gainHistory, sllHistory, 'o', 'filled');
xlabel('Gain (dBi)');
ylabel('Side Lobe Level (SLL) in dB');
title('SLL vs. Gain');
grid on;


% --- Supporting Function Definitions ---

function [fitness, gain, sll] = evaluateFitnessAntennaToolbox(antennaParameters, frequency_Hz, maxSLL, numDirectors, elementRadius)
    
    
    % Extract antenna parameters
    reflectorLength = antennaParameters(1);
    drivenLength = antennaParameters(2);
    directorLength = antennaParameters(3:6);
    spacings = antennaParameters(7:end);

    try
        % Define Yagi elements
        d = dipoleFolded;
        d.Length = drivenLength;
        d.Width = cylinder2strip(elementRadius);
        d.Spacing = d.Length/60;

        yagidesign = yagiUda;
        yagidesign.Exciter = d;
        yagidesign.NumDirectors = numDirectors;
        yagidesign.ReflectorLength = reflectorLength;
        yagidesign.DirectorLength = directorLength(1:4);
        yagidesign.ReflectorSpacing = spacings(1); 
        yagidesign.DirectorSpacing = spacings(2:5);

        % Perform antenna pattern analysis
        [patternData, azimuth, elevation] = pattern(yagidesign, frequency_Hz);

        % Calculate Gain
        gain = max(patternData(:)); % Maximum gain from the radiation pattern

        % Calculate Side Lobe Level (SLL)
        D = pattern(yagidesign,frequency_Hz,0,0:1:360);
        p = polarpattern(D);
    
        lobeInfo = findLobes(p);
        sll = -lobeInfo.SLL;

        % Fitness Calculation (High Gain, Low SLL)
        penaltySLL = max(0, sll - maxSLL);
        fitness = gain - (0.5 * penaltySLL);

    catch err
        disp(['Error in fitness evaluation: ', err.message]);
        fitness = -inf; % Assign very low fitness if evaluation fails
    end
end







