%% Main file
% shows the functionality of the seizureDetection classes and functions

%% Download iEEG measurement from iEEG.org

channel = 4; % channel number to dowlaod
loginName = 'login_name';
loginBin = 'ieeglogin.bin';
studyName = 'Study 005';

% download the iEEG data and store it in a measurement object
measurement = IEEG_getData(studyName, channel, loginName, loginBin); 

% downsample measurement
measurement.downsample(2)

% do a deep copy of the measuremnt object
newMeasurent = measurement.copy(); % do a deep copy

% trim the measurement to keep the two first hours only
newMeasurent.cut(1, 2*3600*measurement.Fs)

% save the measurement object
savename = 'dummyName';
save(savename, 'measurement', '-v7.3')

% plot the measurement
% can be very slow for large measurements
measurement.plot()

%% Download data in a loop and save it

numChannels = 16;
loginName = 'login_name';
loginBin = 'ieeglogin.bin';
studyName = 'Study 005';
savename_truncated = 'Study_005_channel';
for channel = 1:numChannels
    measurement = IEEG_getData(studyName, channel, loginName, loginBin); 
    name = [savename_truncated num2str(channel)];
    save(name, 'measurement', '-v7.3')
end

%% Preprocess the measurement

% customize filter order and cutoff frequency
measurement.FilterOrder = 6; 
measurement.FilterFc1 = 1;
measurement.FilterFc2 = 70;

% filter data, original data will be replaced
measurement.preprocess(); 

%% Calculate features

% lineLength feature
lineLength = LineLength(measurement); % use default values for parameters
lineLength.calculate();
lineLength.calculateOperatingCharacteristicCurve();
lineLength.findBestThreshold();

% also available: customize parameters, see FeatureCalculation object help
% for more details
lineLength = LineLength(measurement, ...
    'StepSize_seconds', 0.5, 'WindowLength_seconds', 4);

% nonlinear energy
nonlinearEnergy = NonlinearEnergy(measurement);
nonlinearEnergy.calculate();
nonlinearEnergy.calculateOperatingCharacteristicCurve();
nonlinearEnergy.findBestThreshold();

% see SpectralBandPower object for more details
thetaBand = SpectralBandPower(measurement, 'ThetaBand', [4 8]);
alphaBand = SpectralBandPower(measurement, 'AlphaBand', [14 32]);
betaBand = SpectralBandPower(measurement, 'BetaBand', [8 12]);

% acces sensitivity & false alarm rate
% sensitivity = lineLength.Sensitivity;
% falseAlarmRate = lineLength.FalseAlarmRate;

% plot feature
lineLength.plot('Threshold')
% with more options and more information displayed on plot
% lineLength.plot('Measurement', 'Threshold', 'Baseline', 'Mark')

% plot sensitivity false alarm rate curve
[h, ax] = lineLength.CharacteristicCurve.plot('Legend String'); % create new figure
lineLength.CharacteristicCurve.addCurve(h, 'Legend String 2'); % add a plot to existing figure

%% Generate feature matrix

% call script
generateMatrix

% or just load it
matrixPath = 'Study_005_matrix.mat';
load(matrixPath)

% plot one channel or one feature of the matrix
matrix.plot('feature', 4, 'Threshold', 'Bare', 'Mark')
matrix.plot('channel', 1, 'Threshold', 'Mark')

%% Matrix detection
mdOLS = MatrixDetection_LinearRegression(matrix, 'Name', 'OLS', 'WeightSeizure', 50);

% set function handle to least squares solution
mdOLS.CalculationFunctionHandle = @(X, y) (X'*X + 0.01*eye(size(X'*X)))\X'*y;
mdOLS.calculateWeights()
mdOLS.calculate()

% display area under curve (up to 1 false alarm per hour)
mdOLS.CharacteristicCurve.calculateArea(1)

% display the weighted sum
mdOLS.plotWeightedDetection()
mdOLS.plotWeightedDetection('Threshold', 'Mark')

% plot characteristic curve of weighted detection together with the best
% feature
[h, ax] = mdOLS.CharacteristicCurve.plot();
matrix.bestFeature().CharacteristicCurve.addCurve(h);

% get detection latency
mdOLS.calculateDetectionLatency();

%% CrossValidate Matrix
numFolds = 4;
[curves_training, curves_validation] = MatrixCrossValidation(matrix, 15, numFolds);
save('curves', 'curves_training', 'curves_validation', '-v7.3')

% compare training and validation curves (plot them)
area_validation = mean(cellfun(@(x) x.calculateArea, curves_validation));
area_training = mean(cellfun(@(x) x.calculateArea, curves_training));

sensitivity_training = zeros(size(curves_training{1}.Sensitivity));
sensitivity_validation =  zeros(size(curves_validation{1}.Sensitivity));
falseAlarmRate_training = zeros(size(curves_training{1}.FalseAlarmRate));
falseAlarmRate_validation =  zeros(size(curves_validation{1}.FalseAlarmRate));
for k = 1:numFolds
    sensitivity_training = sensitivity_training + curves_training{k}.Sensitivity/numFolds;
    sensitivity_validation =  sensitivity_validation + curves_validation{k}.Sensitivity/numFolds;
    falseAlarmRate_training = falseAlarmRate_training + curves_training{k}.FalseAlarmRate/numFolds;
    falseAlarmRate_validation = falseAlarmRate_validation + curves_validation{k}.FalseAlarmRate/numFolds;
end % for k
fig = figure();
plot(falseAlarmRate_training, sensitivity_training, 'k')
hold on
plot(falseAlarmRate_validation, sensitivity_validation, 'r')
hold off
xlabel('False Alarm Rate [false alarms/hour]')
ylabel('Sensitivity')
legend('Training', 'Validation')
xlim([0 1])
increaseSize(fig)
formatplot(fig)

%% Live Demo
% run the live demo
load('./liveDemo/demoData.mat')
liveDemo(data, fs, curveData_matrix, curveData_bestFeature)

%% Lasso Analysis
% script
lassoAnalysis