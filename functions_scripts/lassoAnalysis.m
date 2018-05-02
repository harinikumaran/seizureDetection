%% Lasso Analysis
% Determine feature weights through LASS0 regression
% this has the advantage that many feature weights will be 0, which means
% these features a not informative and can be turned of inside the chip to
% potentially save power

%% Load data
matrixpath = '/Matrices/Study_005_matrix.mat';
load(matrixpath)

%% Generate feature matrix and target vector
weightSeizure = matrix.Matrix{1}.CostSensitivity;
weightInterictal = matrix.Matrix{1}.CostFalseAlarmRate;

nSeizures = matrix.NumSeizures;
nIntericalSamples = sum(~matrix.Matrix{1}.getBitIsSeizure);

% remove seizure areas for now
interical = double(matrix.FlagMatrix);
interical = reshape(interical, matrix.NumElements, matrix.FeatureLength)';
interical(matrix.Matrix{1}.getBitIsSeizure, :) = [];

% compact every seizure into one sample so that there is only on sample per 
% seizure
% either 1 if that seizure was detected or 0 otherwise
detected = reshape(matrix.BitSeizureDetectedMatrix, matrix.NumElements, matrix.NumSeizures)';

% add compacted seizures at the end
X = [sqrt(weightInterictal)*interical; detected*sqrt(weightSeizure)];

% generate target vector
y = [zeros(nIntericalSamples, 1); sqrt(weightSeizure)*ones(nSeizures, 1)]; 

%% Run LASSO regression
% calculate lasso weights for different values of lambda
% this takes a long time to run
[B, fitInfo] = lasso(X, y);

%% Analyse results
% for every value of lambda, calculate the area under curve

nLambdas = size(B, 2);
md = MatrixDetection(matrix, 'Name', 'Lasso');

sensitivities = [];
falseAlarmRates = [];
area1 = zeros(nLambdas, 1);
area2 = zeros(nLambdas, 1);
area5 = zeros(nLambdas, 1);
for k = 1:nLambdas
    fprintf('*** Loop index %d - Lambda = %.10f\n', k, fitInfo.Lambda(k))
    fprintf('*** Number of zero elements: %d\n', fitInfo.DF(k))
    weights = B(:, k);
    md.WeightMatrix(:) = weights;
    md.calculate();
    sensitivities = [sensitivities md.CharacteristicCurve.Sensitivity]; %#ok<*AGROW>
    falseAlarmRates = [falseAlarmRates md.CharacteristicCurve.FalseAlarmRate];
    area1(k) = md.CharacteristicCurve.calculateArea(1); % area under curve up to 1 false alarm per house
    area2(k) = md.CharacteristicCurve.calculateArea(2); % area under curve up to 2 false alarm per house
    area5(k) = md.CharacteristicCurve.calculateArea(5); % area under curve up to 5 false alarm per house
end % for k

%% Plot results

f1 = figure(1);

plot(fitInfo.DF, area1)
xlabel('Number of nonzero coefficients')
ylabel('Area under Curve')
title('Lasso Linear Regression Analysis')
grid on
xlim([0 size(B, 1)])
f1.Color = 'w';
increaseSize(f1);