%% Latency Analysis

%% Load data
matrixpath = '/Users/paul/Google Drive/Microchip_Biosignal_Computation/Seizure_Data/Matrices/Study_005_matrix.mat';
load(matrixpath)

% create MatrixDetection object
mdOLS = MatrixDetection_LinearRegression(matrix, 'Name', 'OLS', 'WeightSeizure', 50);

% set function handle to least squares solution
mdOLS.CalculationFunctionHandle = @(X, y) (X'*X + 0.01*eye(size(X'*X)))\X'*y;
mdOLS.calculateWeights()
mdOLS.calculate()

%% Calculate latencies for different threshold values

thresholds = mdOLS.CharacteristicCurve.ThresholdFactor;
latency = zeros(mdOLS.NumSeizures, length(thresholds));
sensitivity = mdOLS.CharacteristicCurve.Sensitivity;
falseAlarmRate = mdOLS.CharacteristicCurve.FalseAlarmRate;

for k = 1:length(thresholds)
    mdOLS.Threshold = thresholds(k);
    mdOLS.weightedDetection;
    latency(:, k) = mdOLS.calculateDetectionLatency();
end % for

%% Plot results

indices = [12, 19, 34, 41];
fig = figure();

labels = cell(length(indices), 1);
for k = 1:length(labels)
    labels{k} = [num2str(sensitivity(indices(k)), '%.2f') '; ' num2str(falseAlarmRate(indices(k)), '%.2f')];
end % for
boxplot(latency(:, indices), 'Labels',labels)
hold on
plot(mean(latency(:, indices), 1, 'omitnan'), 'k+', 'MarkerSize', 20)
hold off
ylabel('Latency [seconds]')
xlabel('Sensitivity; False Alarm Rate')
legend('Mean')
grid on
ylim([0 50])
increaseSize(fig)
export_fig('latency', '-pdf')
export_fig('latency', '-png', '-m2')