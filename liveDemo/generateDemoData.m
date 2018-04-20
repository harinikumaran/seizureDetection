% generate data for liveDemo

%% Load Relevant data

measurementPath = '/Users/paul/Google Drive/Microchip_Biosignal_Computation/Seizure_Data/Study_005/Study_005_channel1.mat';
load(measurementPath)

matrixPath = '/Users/paul/Google Drive/Microchip_Biosignal_Computation/Seizure_Data/Matrices/Study_005_matrix.mat';
load(matrixPath)

%% Generate Features
measurement.downsample(2);

fs = measurement.Fs;
indices_meas = measurement.SeizureStart(2)-10*fs:measurement.SeizureEnd(2)+20*fs;

measurement.preprocess();

lineLength = LineLength(measurement);
lineLength.calculate();

nonlinearEnergy = NonlinearEnergy(measurement);
nonlinearEnergy.calculate();

alphaBand = SpectralBandPower(measurement, 'AlphaBand', [14 32]);
alphaBand.calculate();

betaBand = SpectralBandPower(measurement, 'BetaBand', [8 12]);
betaBand.calculate();

thetaBand = SpectralBandPower(measurement, 'ThetaBand', [4 8]);
thetaBand.calculate();

mdOLS = MatrixDetection_LinearRegression(matrix, 'Name', 'OLS');
mdOLS.CalculationFunctionHandle = @(X, y) (X'*X + 0.01*eye(size(X'*X)))\X'*y;
mdOLS.calculateWeights;
mdOLS.calculate()

%% Extract data for demo
stepSize = lineLength.StepSize;
data = zeros(length(indices_meas), 8);
isSeizure = measurement.isSeizure(indices_meas);
data(:, 1) = measurement.Data(indices_meas);
data(:, 2) = measurement.Data(indices_meas);
data(isSeizure, 1) = NaN;
data(~isSeizure, 2) = NaN;

for k = 1:length(indices_meas)
    index = floor(indices_meas(k)/stepSize)+1;
    data(k, 3) = lineLength.Value(index);
    data(k, 4) = nonlinearEnergy.Value(index);
    data(k, 5) = alphaBand.Value(index);
    data(k, 6) = betaBand.Value(index);
    data(k, 7) = thetaBand.Value(index);
    data(k, 8) = mdOLS.Value(index); 
end % for k

data(:, 3) = data(:, 3)/max(data(:, 3));
data(:, 4) = data(:, 4)/max(data(:, 4));
data(:, 5) = data(:, 5)/max(data(:, 5));
data(:, 6) = data(:, 6)/max(data(:, 6));
data(:, 7) = data(:, 7)/max(data(:, 7));
data(:, 8) = data(:, 8)/100;

%% Generate data for characteristic curves data
curveData_matrix = [mdOLS.CharacteristicCurve.ThresholdFactor'/100;  
    mdOLS.CharacteristicCurve.Sensitivity';
    mdOLS.CharacteristicCurve.FalseAlarmRate']'; 
   
curveData_bestFeature = [matrix.bestFeature().CharacteristicCurve.ThresholdFactor'/100;
    matrix.bestFeature().CharacteristicCurve.Sensitivity'; 
    matrix.bestFeature().CharacteristicCurve.FalseAlarmRate']';

%% Save data

save('demoData', 'data', 'curveData_matrix', 'curveData_bestFeature', 'fs')