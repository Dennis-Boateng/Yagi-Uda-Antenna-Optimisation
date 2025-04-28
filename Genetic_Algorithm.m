% Genetic Algorithm for Yagi-Uda Antenna Optimisation
clc; clear all;

% Parameters
populationSize = 100; % Number of individuals
numElements = 6; % Number of antenna elements
numGenes = 11; % Total number of genes (lengths + spacings)
frequency_Hz = 165e6; % Operating frequency
speedOfLight_m_s = 3e8;
maxSLL = -15; % Maximum Side Lobe Level (dB)
bitsPerGene = 7; % Bits per gene
mutationRate = 0.1; 
numGenerations = 100;
numDirectors = 4;

% Wavelength
wavelength_m = speedOfLight_m_s / frequency_Hz;

% Parameter Ranges (in wavelengths)
reflectorLengthRange = [0.5, 0.55];
drivenLengthRange = [0.45, 0.49];
directorLengthRange = [0.35, 0.44];
spaceLengthRange = [0.1, 0.45];

% Convert ranges to meters
reflectorLengthRange_m = reflectorLengthRange * wavelength_m;
drivenLengthRange_m = drivenLengthRange * wavelength_m;
directorLengthRange_m = directorLengthRange * wavelength_m;
spaceLengthRange_m = spaceLengthRange * wavelength_m;
elementRadius = 0.003 * wavelength_m;

tic;  % Start timing

% Initialize Population
population_binary = initializePopulationBinary(populationSize, numGenes, bitsPerGene, ...
    reflectorLengthRange_m, directorLengthRange_m, drivenLengthRange_m, spaceLengthRange_m);

goalFitness = 16;  % Define target fitness
fitnessHistory = [];
sllHistory = [];
gainHistory = [];

% --- Main GA Loop ---
for generation = 1:numGenerations

    % Decode the binary population
    population = decodePopulation(population_binary, numGenes, bitsPerGene, ...
        reflectorLengthRange_m, drivenLengthRange_m, directorLengthRange_m, spaceLengthRange_m);
    
    % Evaluate fitness
    [fitness, bestSLL, bestGain] = evaluateFitnessAntennaToolbox(population, frequency_Hz, maxSLL, numDirectors, elementRadius);
    fitnessHistory(end+1) = max(fitness);
    sllHistory(end+1) = bestSLL;
    gainHistory(end+1) = bestGain;

    % Check if goal fitness achieved
    if max(fitness) >= goalFitness
        fprintf('Goal fitness reached at generation %d.\n', generation);
        break;
    end

    % Elitism
    eliteCount = floor(0.1 * populationSize); % Top 10% are elites
    [sortedFitness, sortedIndices] = sort(fitness, 'descend');
    elitePopulationBinary = population_binary(sortedIndices(1:eliteCount), :);

    % Selection
    rouletteIndices = rouletteWheelSelection(fitness, sortedIndices, eliteCount, populationSize - eliteCount);
    selectedNonEliteBinary = population_binary(rouletteIndices, :);

    % Crossover
    offspringPopulationBinary = singlePointCrossover(selectedNonEliteBinary);

    % Mutation
    offspringPopulationBinary = mutatePopulation(offspringPopulationBinary, mutationRate);

    % Create new population
    population_binary = [elitePopulationBinary; offspringPopulationBinary];

    % Display progress
    bestFitness = max(fitness);
    fprintf('Generation %d: Best Fitness = %.2f\n', generation, bestFitness);

    % Check if goal fitness is met
    if bestFitness >= goalFitness
        break;  
    end
end

% Final Results
finalPopulation = decodePopulation(population_binary, numGenes, bitsPerGene,...
    reflectorLengthRange_m, drivenLengthRange_m, directorLengthRange_m, spaceLengthRange_m);
[fitness, bestSLL, bestGain] = evaluateFitnessAntennaToolbox(finalPopulation, frequency_Hz, maxSLL, numDirectors, elementRadius);
[bestFitness, bestIndex] = max(fitness);
bestAntenna = finalPopulation(bestIndex, :);

elapsedTime = toc;  % Stop timing

% Display final results
disp('Best Antenna Configuration (meters):');
disp(bestAntenna);
disp(['Best Fitness: ', num2str(bestFitness)]);
disp(['SLL (dB): ', num2str(bestSLL)]);
disp(['Time to reach goal: ', num2str(elapsedTime), ' seconds']);

d = dipoleFolded;
d.Length = bestAntenna(2);
d.Width = cylinder2strip(elementRadius);
d.Spacing = d.Length/60;

yagidesign = yagiUda;
yagidesign.Exciter = d;
yagidesign.NumDirectors = numDirectors;
yagidesign.ReflectorLength = bestAntenna(1);
yagidesign.DirectorLength = bestAntenna(3:6);
yagidesign.ReflectorSpacing = bestAntenna(7);
yagidesign.DirectorSpacing = bestAntenna(8:11);

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

% Plot Best Antenna Gain vs Iterations
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
title('SLL Convergence');
grid on;

figure;
scatter(gainHistory, sllHistory, 'o', 'filled');
xlabel('Gain (dBi)');
ylabel('Side Lobe Level (SLL) in dB');
title('SLL vs. Gain');
grid on;


% --- Supporting Function Definitions ---

function population_binary = initializePopulationBinary(populationSize, numGenes, bitsPerGene, ...
    reflectorLengthRange_m, directorLengthRange_m, drivenLengthRange_m, spaceLengthRange_m)
    
    population = zeros(populationSize, numGenes);
    population_binary = zeros(populationSize, numGenes * bitsPerGene);

    for i = 1:populationSize
        % Assign random lengths and spacings within respective ranges
        population(i, 1) = reflectorLengthRange_m(1) + (rand() * (reflectorLengthRange_m(2) - reflectorLengthRange_m(1)));
        population(i, 2) = drivenLengthRange_m(1) + (rand() * (drivenLengthRange_m(2) - drivenLengthRange_m(1)));

        population(i, 3) = directorLengthRange_m(1) + rand() * (directorLengthRange_m(2) - directorLengthRange_m(1));

        for j = 4:6
            decayFactor = 0.95 + rand() * (0.98 - 0.95);
            value = decayFactor * population(i, j - 1);
            population(i, j) = max(directorLengthRange_m(1), value);  % Clamp to min
        end

        for j = 7:11
            population(i, j) = spaceLengthRange_m(1) + (rand() * (spaceLengthRange_m(2) - spaceLengthRange_m(1)));
        end
        
        % Convert to binary representation
        for j = 1:numGenes
            if j == 1
                geneRange = reflectorLengthRange_m;
            elseif j == 2
                geneRange = drivenLengthRange_m;
            elseif j >= 3 && j <= 6
                geneRange = directorLengthRange_m;
            else
                geneRange = spaceLengthRange_m;
            end
            scaledValue = (population(i, j) - geneRange(1)) / (geneRange(2) - geneRange(1));
            quantizedValue = round(scaledValue * (2^bitsPerGene - 1));
            quantizedValue = min(quantizedValue, 2^bitsPerGene - 1);  % Clip to valid max
            binaryString = dec2bin(quantizedValue, bitsPerGene);
            population_binary(i, (j - 1) * bitsPerGene + 1 : j * bitsPerGene) = binaryString - '0';
        end
    end
end

function population = decodePopulation(population_binary, numGenes, bitsPerGene,...
    reflectorLengthRange_m, drivenLengthRange_m, directorLengthRange_m, spaceLengthRange_m)

    populationSize = size(population_binary, 1);
    population = zeros(populationSize, numGenes);

    for i = 1:populationSize
        for j = 1:numGenes
            if j == 1
                geneRange = reflectorLengthRange_m;
            elseif j == 2
                geneRange = drivenLengthRange_m;
            elseif j >= 3 && j <= 6
                geneRange = directorLengthRange_m;
            else
                geneRange = spaceLengthRange_m;
            end
            startIndex = round((j - 1) * bitsPerGene + 1);
            endIndex = round(startIndex + bitsPerGene - 1);
            binaryString = char(population_binary(i, startIndex:endIndex) + '0');
            quantizedValue = bin2dec(binaryString);
            scaledValue = quantizedValue / (2^bitsPerGene - 1);
            population(i, j) = geneRange(1) + (scaledValue * (geneRange(2) - geneRange(1)));
        end
    end
end

function [fitness, bestSLL, bestGain] = evaluateFitnessAntennaToolbox(population, frequency_Hz, maxSLL, numDirectors, elementRadius)
  
    % Evaluate fitness using MATLAB's Antenna Toolbox
    populationSize = size(population, 1);
    fitness = zeros(populationSize, 1);
    bestGain = 0;  % Initialize to worst possible gain
    bestSLL = 20;    % Initialize to worst possible SLL

    for i = 1:populationSize
        % Extract antenna parameters
        antennaParameters = population(i, :);
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
            fitness(i) = gain - (0.5 * penaltySLL);

            % **Separate Best Gain Selection**
            if gain > bestGain
                bestGain = gain;
            end

            % **Separate Best SLL Selection**
            if sll < bestSLL
                bestSLL = sll;
            end

        catch err
            disp(['Error in fitness evaluation: ', err.message]);
            fitness(i) = -inf; % Assign very low fitness if evaluation fails
        end
    end
end

function selectedIndices = rouletteWheelSelection(fitness, sortedIndices, eliteCount, numSelections)
    % Exclude elite individuals from the fitness array
    nonEliteFitness = fitness(sortedIndices(eliteCount + 1:end)); 
    
    % Compute selection probabilities based on non-elite fitness
    totalFitness = sum(nonEliteFitness); 
    probabilities = nonEliteFitness / totalFitness;
    cumulativeProbabilities = cumsum(probabilities);

    selectedIndices = zeros(numSelections, 1);
    
    for i = 1:numSelections
        randomNumber = rand();
        selectedIndices(i) = find(cumulativeProbabilities >= randomNumber, 1, 'first');
    end
end

function offspringPopulationBinary = singlePointCrossover(selectedNonEliteBinary)
    populationSize = size(selectedNonEliteBinary, 1);
    chromosomeLength = size(selectedNonEliteBinary, 2);
    offspringPopulationBinary = selectedNonEliteBinary;

    for i = 1:2:populationSize - 1
        if i + 1 <= populationSize
            crossoverPoint = randi(chromosomeLength - 1);
            offspringPopulationBinary(i, :) = [selectedNonEliteBinary(i, 1:crossoverPoint), ...
                                               selectedNonEliteBinary(i + 1, crossoverPoint + 1:end)];
            offspringPopulationBinary(i + 1, :) = [selectedNonEliteBinary(i + 1, 1:crossoverPoint), ...
                                                   selectedNonEliteBinary(i, crossoverPoint + 1:end)];
        end
    end

    % Handle last individual separately if Non Elite population size is odd
    if mod(size(selectedNonEliteBinary, 1), 2) ~= 0
        randomIndex = randi(size(selectedNonEliteBinary, 1) - 1); % Select a random individual
        crossoverPoint = randi(chromosomeLength - 1);
        offspringPopulationBinary(end, :) = [selectedNonEliteBinary(end, 1:crossoverPoint), ...
                                             selectedNonEliteBinary(randomIndex, crossoverPoint + 1:end)];
    end
end

function offspringPopulationBinary = mutatePopulation(offspringPopulationBinary, mutationRate)
    populationSize = size(offspringPopulationBinary, 1);
    chromosomeLength = size(offspringPopulationBinary, 2);

    for i = 1:populationSize
        if rand() < mutationRate
            bitToFlip = randi(chromosomeLength);
            offspringPopulationBinary(i, bitToFlip) = 1 - offspringPopulationBinary(i, bitToFlip);
        end
    end
end


