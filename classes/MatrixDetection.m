classdef MatrixDetection < handle
    % MATRIXDETECTION is a class that handles weighted sum detection of
    % seizures
    % The idea is to look at several channel and features at once. These
    % are stored in the FeatureMatrix object. Each of the entries of the
    % matrix will be either zero (value below threshold) or 1 (value above
    % threshold). This object allows to attribute a weight to each element
    % of the matrix and then multiply the value (0 or 1) of each matrix 
    % element to obtain a weighted sum. This weighted sum can be compared
    % to a threshold to detect seizures.
    %
    % Some parts of this object are similar to the FeatureCalculation
    % object, becauser in a sense they are similar, expect that the
    % threshold here is a constant (no baseline) and the value is a 
    % weighted sum
    %
    % The construct has one required input parameter which is a
    % FeatureMatrix object
    % Several Name-Value pairs are possible for customization:
    %   -CostSensitivity: for determining the best threshold, higher values
    %       mean that detecting seizures is more important than avoiding 
    %       false alarms
    %   -CostFalseAlarmRate: for determining the best threshold, higher
    %       values mean that less false alarms are desired 
    %   -ThresholdFactorValues: the possible values for the threshold.
    %       Depending on how the weights are calculated, the weighted sum 
    %       could be negative or only positive and it could have a 
    %       different range.

    
    properties
        Name
        
        WeightMatrix % matrix of weights associated with each feature
        
        Value % weighted sum value
        Threshold % threshold
        LinkedMatrix % linked FeatureMatrix object
        
        CostSensitivity % factor for calculating cost when finding best threshold
        CostFalseAlarmRate % factor for calculating cost when finding best threshold
        
        ThresholdFactorValues % for calculating the operating characteristic curve
    end
    
    properties (SetAccess = protected)
        Ranking % a ranking the features by their respective weight 
        
        CharacteristicCurve % holds an OperatingCharacteristicCurve object
        
        NumElements % number of elements in the matrix
        
        Sensitivity % Percentage of detected seizures
        Specificity % Percentage of time above threshold outside of seizure times
        FalseAlarms % Number of false alarms
        FalseAlarmRate % Number of false alarms per hour
        ThresholdCrossings % Number of times the feature value exceeded the threshold
        MissedSeizures % Numeration of seizure that were not detected
        BitSeizureDetected % Logical array, index i is true when seizure i is detected by this feature, false otherwise 
        ThresholdCrossedLocation % Indices where the feature value crossed the threshold, in terms of feature samples
        FalseAlarmLocation % Indices of false alarms, in terms of feature samples
        
        FeatureLength % Total number of calculated values
        SeizureStart % In terms of feature samples
        SeizureEnd % In terms of feature samples
        NumSeizures % Total number of seizures
        SeizureHoldTime % duration of a seizure for evaluation purpose in feature samples
    end
    
    methods
        % constructor
        function obj = MatrixDetection(varargin)
            % first input argument should be a FeatureMatrix object
            
            % parse input
            defaultCostSensitivity = 50;
            defaultCostFalseAlarmRate = 1;
            defaultName = 'Generic Matrix Detection';
            defaultThresholdFactorValues = 0:0.05:100;
            
            p = inputParser;
            p.KeepUnmatched = true;
            addRequired(p, 'matrix', @(x) validateattributes(x,{'FeatureMatrix'},{'nonempty'},mfilename,'matrix',1));
            addParameter(p,'CostSensitivity',defaultCostSensitivity, @isscalar);
            addParameter(p,'CostFalseAlarmRate',defaultCostFalseAlarmRate, @isscalar);
            addParameter(p, 'ThresholdFactorValues', defaultThresholdFactorValues, @isnumeric)
            addParameter(p, 'Name', defaultName, @ischar);
            parse(p,varargin{:});
            matrix = p.Results.matrix;
            
            obj.LinkedMatrix = matrix;
            
            obj.FeatureLength = matrix.FeatureLength;
            obj.SeizureStart = matrix.Matrix{1}.SeizureStart;
            obj.SeizureEnd = matrix.Matrix{1}.SeizureEnd;
            obj.NumSeizures = matrix.NumSeizures;
            obj.SeizureHoldTime = matrix.Matrix{1}.SeizureHoldTime;
            
            obj.NumElements = matrix.NumElements;
            
            obj.WeightMatrix = zeros(size(matrix.Matrix));

            obj.CostSensitivity = p.Results.CostSensitivity;
            obj.CostFalseAlarmRate = p.Results.CostFalseAlarmRate;
            
            obj.Name = p.Results.Name;
            
            obj.ThresholdFactorValues = p.Results.ThresholdFactorValues;
        end % constructor
        
        function sobj = saveobj(obj)
            % this functions executes when the object is saved
            % otherwise, the whole feature matrix will be saved
            warning('The linked feature matrix will not be saved! Saving this object is not that usefull...')
            obj.LinkedMatrix = [];
            try
                obj.Ranking = rmfield(obj.Ranking, 'Feature');
            catch
                % nothing
            end % try
            sobj = obj;           
        end % function sobj
        
        function calculate(obj)
            % calculates the rank of the features, the weighted sum and the
            % operating characteristic curve

            if isempty(obj.WeightMatrix)
                warning('Weights not calculated. Not doing anything.')
                return
            elseif all(obj.WeightMatrix(:) == 0)
                warning('All weights are 0.')
            end
            
            % rank each feature
            obj.rankFeatures();
            
            % calculate value
            obj.calculateValue(); 
            
            % calculate operating characteristic curve
            obj.calculateOperatingCharacteristicCurve();
                     
        end % function calculate
                
        function weightedDetection(obj, varargin)
            % evaluate the weighted sum detection, calculate sensitivity,
            % false alarm rate and other metrics
            
            if isempty(obj.WeightMatrix)
                error('WeightMatrix not yet calculated. Run calculateWeights function first.')
            end % if
            
            printResults = false;
            if nargin > 1
                if any(strcmp(varargin, 'Print'))
                    printResults = true;
                end
            end % if nargin > 1
            
            % initialize stuff
            feature = obj.LinkedMatrix.Matrix{1};
            numSeizures = obj.NumSeizures;        
            bitDetected = false(numSeizures, 1);
            isSeizure = false(feature.Length, 1);
            
            % mark seizure as detected if at least one value is above
            % threshold
            for k = 1:numSeizures
                indexSeizure = obj.SeizureStart(k):obj.SeizureEnd(k);
                if any(obj.Value(indexSeizure) > obj.Threshold)
                    bitDetected(k) = true;
                end % if
                
                % also, mark where the seizures are
                isSeizure(indexSeizure) = true;
            end % for k
            sensitivity = sum(bitDetected)/numSeizures;   
            
            % find locations where the threshold is crossed
            aboveThreshold = obj.Value > obj.Threshold;
            
            thresholdCrossed = [diff(aboveThreshold) == 1; false];  % is true when the threshold was crossed between this and the previous sample
            thresholdCrossed(isSeizure) = false;
            
            thresholdCrossings = sum(thresholdCrossed);
            thresholdCrossedLocation = find(thresholdCrossed);
            
            if ~isempty(thresholdCrossedLocation)
                % only do it if the threshold was crossed at least once 
                
                % if multiple threshold crossings occur within a short time
                % (specified by seizureHoldTime) of each other, only count
                % one false alarm
                
                seizureHoldTime = obj.SeizureHoldTime;
                
                loc = thresholdCrossedLocation;
                indexIgnore = false(size(thresholdCrossedLocation));
                k = 1;
                while ~isempty(k)
                    
                    loc = loc - loc(k);
                    indexIgnore = indexIgnore | (loc > 0 & loc <= seizureHoldTime);
                    
                    k = find(loc > seizureHoldTime, 1, 'first');
                    
                end % while k
                
                falseAlarmLocation = thresholdCrossedLocation(~indexIgnore);
            else
                % no false alarms
                falseAlarmLocation = [];                
            end % if ~isempty
            
            % get number of false alarms and their location
            falseAlarms = length(falseAlarmLocation);
            falseAlarmRate = falseAlarms/(feature.Duration/3600);
            
            % Store missed seizures and false alarms location
            missedSeizures = find(~bitDetected);
            bitSeizureDetected = bitDetected;
            
            % Print result to command line
            if printResults
                fprintf('*** Results for weighted detection:\n')
                fprintf('Sensitivity: %3.2f%%\n', sensitivity*100)
                fprintf('Number of false alarms: %d\n', falseAlarms)
                fprintf('False alarm rate: %.6f\n', falseAlarmRate)
            end % if printResults
            
            % store results in object
            obj.Sensitivity = sensitivity;
            obj.ThresholdCrossings = thresholdCrossings;
            obj.ThresholdCrossedLocation = thresholdCrossedLocation;
            obj.FalseAlarms = falseAlarms;
            obj.FalseAlarmRate = falseAlarmRate;
            obj.FalseAlarmLocation = falseAlarmLocation;
            obj.MissedSeizures = missedSeizures;
            obj.BitSeizureDetected = bitSeizureDetected;
            
        end % function weightedDetection
        
        function plotWeightedDetection(obj, varargin)
            % plot the weighted sum
            % possible plot customization options are
            % 'Mark' and 'Threshold'
            
            % check for additional input
            plotThreshold = false;
            bitMark = false;
            if nargin > 1
                if any(strcmp(varargin, 'Mark'))
                    bitMark = true;
                end
                if any(strcmp(varargin, 'Threshold'))
                    plotThreshold = true;
                end
            end % if nargin > 1
            
            h = figure( 'units', 'pixels', ...
                'Name', 'Weighted Score', ...
                'Tag', 'Weighted Score', ...
                'Visible', 'off'); % plotting while not visible is faster
            
            ax = axes;
            
            feature = obj.LinkedMatrix.Matrix{1};
            seizureStart = obj.SeizureStart;
            seizureEnd = obj.SeizureEnd;
            
            % generate time vector
            Fs = feature.MeasurementFs/feature.StepSize;
            Ts = 1/Fs;
            time_feature = 0:Ts:(feature.Length-1)*Ts;
            
            % plot samples where no seizure occurs first
            % replace seizure samples with NaNs
            tempData = obj.Value;
            tempTime = time_feature;
            for k = 1:length(seizureStart)
                tempTime(seizureStart(k):seizureEnd(k)) = nan;
                tempData(seizureStart(k):seizureEnd(k)) = nan;
            end % for k
            l(1) = plot(ax, tempTime, tempData);
            
            % now plot seizure samples in red so they stand out
            tempData = nan(size(tempData));
            tempTime = nan(size(tempTime));
            for k = 1:length(seizureStart)
                tempTime(seizureStart(k):seizureEnd(k)) = time_feature(seizureStart(k):seizureEnd(k));
                tempData(seizureStart(k):seizureEnd(k)) = obj.Value(seizureStart(k):seizureEnd(k));
            end % for k
            hold(ax, 'on')
            l(2) = plot(ax, tempTime, tempData, 'r');
            hold(ax, 'off')
            legendstr = {'Value', 'Seizure'};
            
            if plotThreshold
                hold(ax, 'on')
                l(end+1) = plot(ax, time_feature([1, end]), obj.Threshold([1, 1]));
                hold(ax, 'off')
                legendstr{end+1} = 'Threshold';
            end % if
            
            if bitMark
                hold(ax, 'on')
                % mark false alarms
                falseAlarmLocation = obj.FalseAlarmLocation;
                missedSeizures = obj.MissedSeizures;
                if ~isempty(falseAlarmLocation)
                    l(end+1) = plot(ax, time_feature(falseAlarmLocation), ...
                        obj.Value(falseAlarmLocation), 'kx', 'MarkerSize', 10);
                    legendstr{end+1} = 'False Alarm';
                end % if
                if ~isempty(missedSeizures)
                    indexMissed = seizureStart(missedSeizures);
                    l(end+1) = plot(ax, time_feature(indexMissed), ...
                        obj.Value(indexMissed), 'dk', 'MarkerSize', 10);
                    legendstr{end+1} = 'Missed Seizure';
                end % if
                hold(ax, 'off')
            end % if bitMark
            
            ylabel(ax, 'Weighted score')
            xlabel(ax, 'Time [s]')
            legend(l, legendstr)
            
            formatplot(h);
            
            h.Visible = 'on';
        end % function plotWeightedDetection
        
        function findBestThreshold(obj)
            % find the best threshold through cost function minimization
            % first, the operating characteristic curve is calculated, it
            % contains the sensitivity and false alarm rates for different
            % threshold values (sweep of all possible values)
            
            costSensitivity = obj.CostSensitivity;
            costFalseAlarmRate = obj.CostFalseAlarmRate;
            
            if isempty(obj.CharacteristicCurve)
                warning('Operating Characteristic Curve not calculated yet. Doing it now.')
                obj.calculateOperatingCharacteristicCurve;
            end % if
            
            sensitivity = obj.CharacteristicCurve.Sensitivity;
            falseAlarmRate = obj.CharacteristicCurve.FalseAlarmRate;
            
            % calculate the cost
            cost = costSensitivity*(1-sensitivity).^2 ...
                + costFalseAlarmRate*(falseAlarmRate).^2;
            
            [~, index] = min(cost);
            
            bestFactor = obj.CharacteristicCurve.ThresholdFactor(index);
            
            obj.Threshold = bestFactor;
            obj.weightedDetection('Print')
            
        end % function findBestThreshold
                        
        function calculateValue(obj)
            % calculates the WeightedDetectionScore
            % the score is calculated for every feature sample by
            % taking the sum over all feature of the product
            % bitAboveThreshold*score of every feature.
            % This product will either be 0 if the feature value is not
            % above the threshold for that sample OR it will be the score
            % associated with that feature if the feature value is above
            % the threshold for that sample
            % This calculation can be done very quickly by matrix
            % multiplication, this is why it only takes one line to compute
            % it
            % This function can be overwritten by children class that
            % implement another way of calculating the value
            
            % - obj.WeightMatrix(:) turns the score matrix into a column
            % vector of size [numElements x 1]
            % - reshape(reshape(obj.AboveThresholdMatrix, obj.NumElements, obj.FeatureLength)
            % turn the 3D-matrix that is obj.AboveThresholdMatrix into a
            % 2D-matrix of size [numElements x featureLength] where every
            % row corresponds to one feature
            % - we transpose this matrix so that is has size [featureLength
            % x numElements] and can be (matrix-)multiplied with the score,
            % i.e. the resuls will be a [featureLength x 1] matrix/vector
            obj.Value = transpose(reshape(obj.LinkedMatrix.FlagMatrix, obj.NumElements, obj.FeatureLength))*obj.WeightMatrix(:);
            
            % normalize to values between 0 and 100
            obj.Value = obj.Value*100/max(obj.Value);
            
        end % function
        
        function rankFeatures(obj)
            % rank features according to their weight
            
            % obj.FeatureScore(:) turns the score matrix into a columns vector
            [sortedWeight, indexRunning] = sort(obj.WeightMatrix(:), 'descend');
            
            [indexRow, indexCol] = ind2sub(size(obj.WeightMatrix), indexRunning);
            
            try
                obj.Ranking = struct('Weight', num2cell(sortedWeight), ...
                    'Sensitivity', num2cell(obj.LinkedMatrix.SensitivityMatrix(indexRunning)), ...
                    'FalseAlarmRate', num2cell(obj.LinkedMatrix.FalseAlarmRateMatrix(indexRunning)), ...
                    'Specifiticy', num2cell(obj.LinkedMatrix.SpecificityMatrix(indexRunning)), ...
                    'AreaUnderCurve', num2cell(obj.LinkedMatrix.AreaMatrix(indexRunning)), ...
                    'Channel', num2cell(indexRow), ...
                    'Feature', obj.LinkedMatrix.Matrix(indexRunning), ...
                    'IndexCol', num2cell(indexCol), ...
                    'IndexRunning', num2cell(indexRunning));
            catch
                obj.Ranking = struct('Weight', num2cell(sortedWeight), ...
                    'Sensitivity', num2cell(obj.LinkedMatrix.SensitivityMatrix(indexRunning)), ...
                    'FalseAlarmRate', num2cell(obj.LinkedMatrix.FalseAlarmRateMatrix(indexRunning)), ...
                    'Specifiticy', num2cell(obj.LinkedMatrix.SpecificityMatrix(indexRunning)), ...
                    'Channel', num2cell(indexRow), ...
                    'Feature', obj.LinkedMatrix.Matrix(indexRunning), ...
                    'IndexCol', num2cell(indexCol), ...
                    'IndexRunning', num2cell(indexRunning));
            end % try
        end % function rankFeatures
        
        function calculateOperatingCharacteristicCurve(obj)
            % calculate the receiving operating characteristive curve of
            % this feature by sweeping through multiple threshold values
            % and calculating the sensitivity and false alarm rate for each
            % of them
            
            timeElapsed = tic;
            fprintf('Calculating operating characteristics...\n')
            
            thresholdFactorValues = obj.ThresholdFactorValues;
            sensitivity = zeros(size(thresholdFactorValues));
            falseAlarmRate = zeros(size(thresholdFactorValues));
          
            k = 1;
            while (1)
                obj.Threshold = thresholdFactorValues(k);
                obj.weightedDetection;
                sensitivity(k) = obj.Sensitivity;
                falseAlarmRate(k) = obj.FalseAlarmRate;
                
                if (sensitivity(k) == 0 && falseAlarmRate(k) == 0 ...
                        && sensitivity(k-1) == 0 && falseAlarmRate(k-1) == 0) ...
                        || k == length(sensitivity)
                    
                    % stop when sensitivity and false alarm rate are 0 for
                    % more than two steps
                    
                    break;
                else
                    k = k + 1;
                end % if
            end % while
            
            sensitivity = sensitivity(1:k);
            falseAlarmRate = falseAlarmRate(1:k);
            thresholdFactorValues = thresholdFactorValues(1:k);
            
            obj.CharacteristicCurve = OperatingCharacteristicCurve(obj.Name, ...
                sensitivity, falseAlarmRate, thresholdFactorValues);
            
            fprintf('%s operating characteristic calculated in %.2f seconds.\n', obj.Name, toc(timeElapsed));
            
            obj.findBestThreshold();
        end % function calculateOperatingCharacteristicCurve
        
    end % methods
    
end

