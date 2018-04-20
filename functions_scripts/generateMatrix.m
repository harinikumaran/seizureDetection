% generate matrix script
% generate feature matrix from measurement objects
%% Initialize

savename = 'Matrix';

% files should be named sequencially by channel
% truncated path name, the only thing that should be missing is the channel
% number (and the .mat extension) - it will be added in the loop
filename_truncated = '/Users/paul/Google Drive/Microchip_Biosignal_Computation/Seizure_Data/Study_005/Study_005_channel';

% number of channels
numChannels = 16;

% Features
% names cell array is used to initialize FeatureMatrix object. It should
% contain the names of the different features (must match FeatureName
% property of each feature object)
names = {'Line Length', 'Nonlinear Energy', 'Power', 'Theta Band', 'Alpha Band', 'Beta Band'};
numFeatures = length(names);

% create empty FeatureMatrix object
matrix = FeatureMatrix(numChannels, names);

% create waitbar
lengthWaitbar = numChannels*numFeatures;
loopCounter = 1;
wb = waitbar(0, ['Progress: ' num2str(loopCounter/lengthWaitbar*100, '%.2f') '%'], 'Name', 'Calculating...');
totalTimer = tic;

%% Loop through all channels to create matrix
for channel = 1:numChannels
    
    % load file
    filename = [filename_truncated num2str(channel) '.mat'];
    load(filename);
    % measurement.FilterOrder = 6;
    % measurement.FilterFc1 = 4; % Hz
    % measurement.FilterFc2 = 70; % Hz
    
    % preprocess measurement (Bandpass filter)
    measurement.preprocess();
    
    % for each feature, calculate feature value, find best threshold and
    % save element in FeatureMatrix object
    for featureNumber = 1:numFeatures
        switch featureNumber
            case 1
                feature = LineLength(measurement);
            case 2
                feature = Power(measurement);
            case 3
                feature = NonlinearEnergy(measurement);
            case 4
                feature = SpectralBandPower(measurement, 'ThetaBand', [4 8]); % Theta Band
            case 5
                feature = SpectralBandPower(measurement, 'AlphaBand', [14 32]); % Alpha Band
            case 6
                feature = SpectralBandPower(measurement, 'BetaBand', [8 12]); % Beta Band                
        end % switch
               
        feature.calculate();
        feature.calculateOperatingCharacteristicCurve()
        feature.findBestThreshold();
        
        % save in feature matrix
        matrix.setElement(channel, feature);
        
        % update Waitbar
        loopCounter = loopCounter + 1;
        waitbar(loopCounter/lengthWaitbar, wb, [num2str(loopCounter/lengthWaitbar*100, '%.2f') '%, ' num2str(toc(totalTimer)/60, '%.2f') ' min. elapsed'])
    end % for featureNumber
    
    fprintf('*** Done with channel %d, saving... *** \n', channel)
    save(savename, 'matrix', '-v7.3') % save progress
end % for channel

close(wb)
fprintf('Done in %.2f minutes.\n', toc(totalTimer)/60)