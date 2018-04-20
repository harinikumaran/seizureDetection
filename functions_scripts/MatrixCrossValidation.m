function [curves_training, curves_validation] = MatrixCrossValidation(...
    matrix, chunkSize_minutes, numFolds)
% trains the matrix detection algorithms on parts of the data and evaluates
% it on another part of the data
% returns a cell array of OperatingCharacteristicCurve objects, one for 
% every split of the data (there are numFolds splits in total)

%% Generate chunks of data
% cut the data into chunks, these chunks will be shuffle to perform the
% cross validation
dummyFeature = matrix.Matrix{1};
indCell = cell(0);
index = 1;
chunkSize = (dummyFeature.MeasurementFs/dummyFeature.StepSize)*60*chunkSize_minutes;

while index + chunkSize - 1 < dummyFeature.Length
    endIndex = index + chunkSize - 1;
    while dummyFeature.isseizure(endIndex) && index + 100 - 1 <= dummyFeature.Length
        % do not cut seizures, increase chunk size until seizure ends
        endIndex = endIndex + 100;
    end % while
    
    indCell{end+1} = index:endIndex;
    index = endIndex + 1;
end % while

%% Do Cross validation
numChunks = length(indCell);
randomCellIndices = randperm(numChunks); % shuffle the chunk indices
ratioTraining = 1-1/numFolds;
curves_training = cell(numFolds, 1);
curves_validation = cell(numFolds, 1);

for k = 1:numFolds
    
    % Generate indices that will be used to shuffle the data
    trainingIndices = [];
    validationIndices = [];
    for i = 1:floor(numChunks*ratioTraining)
        trainingIndices = [trainingIndices indCell{randomCellIndices(i)}]; %#ok<*AGROW>
    end % for k
    for i = floor(numChunks*ratioTraining)+1:numChunks
        validationIndices = [validationIndices indCell{randomCellIndices(i)}];
    end
    
    % Shift randomCellIndices so next iteration will take different data
    % for training (this cycles through the data)
    randomCellIndices = circshift(randomCellIndices, floor(numChunks*(1-ratioTraining)));
    
    %% Split Dataset   

    % create validation matrix with shuffled data
    matrix_validation = copy(matrix);
    matrix_validation.cut(validationIndices)
    
    % create training matrix with shuffled data
    matrix_training = copy(matrix);
    matrix_training.cut(trainingIndices)
    
    %% Do training on training matrix
    matrix_training.calculateFeatureCharacteristicCurve();
    cellfun(@findBestThreshold, matrix_training.Matrix)
    matrix_training.fillMatrices();
    matrix_training.fillAreaMatrix();
    
    % matrix detection
    mdOLS_training = MatrixDetection_LinearRegression(matrix_training, 'Name', 'Ridge Training');
    mdOLS_training.CalculationFunctionHandle = @(X, y) (X'*X + 0.01*eye(size(X'*X)))\X'*y;
    mdOLS_training.calculateWeights;
    mdOLS_training.calculate()
    
    %% Do  Validation
    for i = 1:matrix_training.NumElements
        % set threshold of individual features to training threshold
        matrix_validation.Matrix{i}.ThresholdBaselineFactor = matrix_training.Matrix{i}.ThresholdBaselineFactor;
    end % for k
    cellfun(@evaluate, matrix_validation.Matrix)
    matrix_validation.fillMatrices();
    
    % matrix detection
    md_validation = MatrixDetection(matrix_validation, 'Name', 'Validation');
    md_validation.WeightMatrix = mdOLS_training.WeightMatrix; % set weights to training weights
    md_validation.calculate()
   
    %% Get results
    [h, ~] = mdOLS_training.CharacteristicCurve.plot('OLS Training');
    md_validation.CharacteristicCurve.addCurve(h, 'Validation')
    
    curves_training{k} = copy(mdOLS_training.CharacteristicCurve);
    curves_validation{k} = copy(md_validation.CharacteristicCurve);
    
    fprintf('** Iteration %d done.\n', k)
    fprintf('Training area under curve: %.4f\n', mdOLS_training.CharacteristicCurve.calculateArea(1))
    fprintf('Validation area under curve: %.4f\n\n', md_validation.CharacteristicCurve.calculateArea(1))
end