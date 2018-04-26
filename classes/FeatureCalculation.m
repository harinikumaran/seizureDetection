classdef (Abstract) FeatureCalculation < matlab.mixin.Copyable
    % FEATURECALCULATION is a parent class for all individual feature
    % detection classes
    %
    % Each feature is computed over a sliding window. The length of the
    % window is 'WindowLength' original samples longs and gets updated
    % every 'StepSize' samples. The total number of computed values is
    % stored in the property 'Length'. Note: One feature value is also
    % refered to as one feature sample in the various comments in the code.
    %
    % Some more notes:
    % - The baseline is calculate by averaging a certain feature values
    %   over a certain (long) window. However, that window is not directly
    %   next the current value, rather it is lagging behind by a certain
    %   amount of time (for example, by default, the baseline is calculated
    %   over a 3 minutes window that is lagging 2 minutes behind). The
    %   baseline is not updated very often minimize the amount of registers
    %   needed to store values in hardware
    % - When detecting false alarms, only one false alarm can occour every
    %   amount of time specified by the SeizureHoldTime property. This
    %   means that if a false alarm occurs, all other false alarms
    %   happenings that amount of time after the first one are ignored for
    %   evaluation purposes. This is to avoid a high amount of false alarms
    %   resulting from the feature values oscillating about the treshold.
    %   This is also how it would work on a chip, since having multiple
    %   seizures right after each other is not realistic.
    % - A seizure counts as detected if at least one feature value is above
    %   the threshold during the seizure
    % - This class is a children of matlab.mixin. Copyable so that it can be
    %   deep-copied if need by calling copy(feature).
    % 
    % Required input for the constructor:
    %   -measurement: a measurement object
    %
    % Possible Name-Value pairs accepted by the constructor for
    % customization are (for default values, see constructor):
    %   -StepSize_seconds: by how much the sliding window is moved every
    %       step (in seconds)
    %   -WindowLength_seconds: how long the sliding window is (in seconds)
    %   -BaselineAverageTime_seconds: over what time window to average the
    %       feature value to calculate the baseline (in seconds)
    %   -AverageWindowShift_seconds: the amount of time the baseline is
    %       lagging behind the feature value (in seconds)
    %   -BaselineRefreshRate_seconds: how often the baseline is refreshed
    %   -SeizureHoldTime_seconds: see note above
    %   -CostSensitivity: for determining the best threshold, higher values
    %       mean that detecting seizures is more important than avoiding 
    %       false alarms
    %   -CostFalseAlarmRate: for determining the best threshold, higher
    %       values mean that less false alarms are desired 
    %   -ThresholdBaselineFactor: the threshold is calculated as
    %       ThresholdBaselinefactor*Baseline
        
    properties
        Value % Feature value
        
        CharacteristicCurve % hold an OperatingCharacteristicCurve object
        
        StudyName % Study name of original measurement
        StudyChannelNumber % Channel number of original measurement
        StudyChannelName
        
        BaselineWindowLength % In terms of feature step size
        BaselineShift % In terms of feature step size
        BaselineRefreshrate % In terms of feature step size
        
        ThresholdBaselineFactor % factor by which the baseline is multiplied to obtain the threshold
        
        LinkedMeasurement % a handle (reference/link) to the measurement object that is being analysed
        
        PlotHandle % a handle to the plot of the feature (if there is one)
        AxesHandle % a handle array containing the axes handles of the plot
        
        CostSensitivity % factor for calculating cost of not detecting seizures when finding best threshold
        CostFalseAlarmRate % factor for calculating cost of false alarms when finding best threshold
        SeizureHoldTime % duration of a seizure for evaluation purpose in feature samples
        
    end % properties
    
    properties (SetAccess = protected) % Read-Only properties
        Baseline % Feature value baseline, average of feature value over a longer window
        
        StepSize % In terms of measurement samples
        WindowLength % In terms of measurement samples
        Length % Total number of calculated values
        
        SeizureStart % In terms of feature samples
        SeizureEnd % In terms of feature samples
        NumSeizures % Total number of seizures
        MeasurementFs % Sampling rate of linked measurement, useful is case the LinkedMeasurement is deleted
        Duration % Duration of linked measurement, useful is case the LinkedMeasurement is deleted
        
        % Properties used for or calculated by feature evaluation
        Sensitivity % Percentage of detected seizures
        Specificity % Percentage of time above threshold outside of seizure times
        FalseAlarms % Number of false alarms
        ThresholdCrossings % Number of times the feature value exceeded the threshold
        FalseAlarmRate % Number of false alarms per hour
        MissedSeizures % Enumeration of seizure that were not detected
        BitSeizureDetected % Logical array, index i is true when seizure i is detected by this feature, false otherwise 
        ThresholdCrossedLocation % Indices where the feature value crossed the threshold, in terms of feature samples
        FalseAlarmLocation % Indices of false alarms, in terms of feature samples
        
    end % properties (SetAcces = protected)
    
    properties (Abstract)
        FeatureName
    end % properties (Abstract)
    
    methods
        % constructor
        function obj = FeatureCalculation(measurement, varargin)
            
            % parse input
            p = inputParser;
            p.KeepUnmatched = true;
            
            % define default values
            defaultStepSize_seconds = 0.2; 
            defaultWindowLength_seconds = 1;
            defaultBaselineAverageTime_seconds = 180;
            defaultAverageWindowShift_seconds = 120;
            defaultBaselineRefreshRate_seconds = 30;
            defaultCostSensitivity = 50;
            defaultCostFalseAlarmRate = 1;
            defaultSeizureHoldTime_seconds = 60;
            defaultThresholdBaselineFactor = 5;

            addRequired(p, 'measurement', @(x) validateattributes(x,{'Measurement'},{'nonempty'},mfilename,'measurement',1));
            addParameter(p, 'StepSize_seconds', defaultStepSize_seconds, @isscalar);
            addParameter(p, 'WindowLength_seconds', defaultWindowLength_seconds, @isscalar);
            addParameter(p, 'BaselineAverageTime_seconds', defaultBaselineAverageTime_seconds, @isscalar);
            addParameter(p, 'AverageWindowShift_seconds', defaultAverageWindowShift_seconds, @isscalar);
            addParameter(p, 'BaselineRefreshRate_seconds', defaultBaselineRefreshRate_seconds, @isscalar);
            addParameter(p, 'SeizureHoldTime_seconds', defaultSeizureHoldTime_seconds, @isscalar);
            addParameter(p, 'CostSensitivity', defaultCostSensitivity, @isscalar);
            addParameter(p, 'CostFalseAlarmRate', defaultCostFalseAlarmRate, @isscalar);
            addParameter(p, 'ThresholdBaselineFactor', defaultThresholdBaselineFactor, @isscalar);
            
            % check inputs
            parse(p,measurement, varargin{:});
            
            % calculate object properties
            
            % store a link to measurement object
            measurement = p.Results.measurement;
            
            obj.LinkedMeasurement = measurement;
            
            obj.StepSize = floor(measurement.Fs*p.Results.StepSize_seconds);
            obj.WindowLength = floor(measurement.Fs*p.Results.WindowLength_seconds);
            if rem(obj.WindowLength, obj.StepSize) ~= 0
                error('Window length has to be a multiple of step size.')
            end % if
            
            obj.Length = length(1:obj.StepSize:measurement.Length);
            
            downsamplingFactor = obj.StepSize;
            
            obj.SeizureStart = floor(obj.LinkedMeasurement.SeizureStart/downsamplingFactor);
            obj.SeizureEnd = ceil(obj.LinkedMeasurement.SeizureEnd/downsamplingFactor);
            obj.NumSeizures = length(measurement.SeizureStart);
            Fs = measurement.Fs;
            obj.MeasurementFs = Fs;
            obj.Duration = obj.LinkedMeasurement.Duration;
           
            obj.BaselineWindowLength = floor(p.Results.BaselineAverageTime_seconds*Fs/obj.StepSize);
            obj.BaselineShift = floor(p.Results.AverageWindowShift_seconds*Fs/obj.StepSize);
            obj.BaselineRefreshrate = floor(p.Results.BaselineRefreshRate_seconds*Fs/obj.StepSize);
            
            obj.SeizureHoldTime = Fs/obj.StepSize*p.Results.SeizureHoldTime_seconds;
            
            obj.CostSensitivity = p.Results.CostSensitivity;
            obj.CostFalseAlarmRate = p.Results.CostFalseAlarmRate;
            
            obj.StudyName = measurement.StudyName;
            obj.StudyChannelNumber = measurement.StudyChannelNumber;
            obj.StudyChannelName = measurement.StudyChannelName;
            
            obj.ThresholdBaselineFactor = p.Results.ThresholdBaselineFactor;
            
        end % constructor
        
        function calculate(obj)
            % calculate feature value
            
            if isempty(obj.LinkedMeasurement)
                warning('Object %s has no measurement linked to it. No calculations done.', obj.FeatureName)
                return
            end % if
            
            disp(['Calculating ' obj.FeatureName '...'])
            
            calcTimer = tic;
            
            % initialize vectors
            obj.Value = zeros(obj.Length, 1);
            
            % do feature specific calculations
            obj.calculateFeatureValue;
            
            % replace NaNs entries with their nearest non NaN
            % neighbor
            if any(isnan(obj.Value))
                x = transpose(1:obj.Length);
                obj.Value = interp1(x(~isnan(obj.Value)), ...
                    obj.Value(~isnan(obj.Value)), x, 'nearest', 'extrap');
            end % if
            
            % calculate baseline
            obj.calculateBaseline()
            
            disp([obj.FeatureName ' calculated in ' num2str(toc(calcTimer), '%.2f') ' seconds.'])
            
        end % function calculate
        
        function evaluate(obj, varargin)
            % evaluate detection algorithm with current threshold factor
            % calculates the sensitivity and false alarm rate as well as
            % some other usefull metrics
            % optional input argument 'Print' can be passed to print out
            % the results of the detection in the command line
            
            % check if feature was calculated before call to this function
            if isempty(obj.Value)
                error(['The value of ' obj.FeatureName ' has not been calculated yet.'])
            end % if
            
            % check input arguments for 'Print' option
            printResults = false;
            if nargin > 1
                if any(strcmp(varargin, 'Print'))
                    printResults = true;
                end
            end % if nargin > 1
            
            % calculate sensitivity:
            % initialize some values
            numSeizures = length(obj.SeizureStart);
            bitDetected = false(numSeizures, 1);
            isSeizure = false(obj.Length, 1);
            
            % mark seizure as detected if at least one value is above
            % threshold
            for k = 1:numSeizures
                indexSeizure = obj.SeizureStart(k):obj.SeizureEnd(k);
                if any(obj.Value(indexSeizure) > obj.Threshold(indexSeizure))
                    bitDetected(k) = true;
                end % if
                
                % also, mark where the seizures are
                isSeizure(indexSeizure) = true;
            end % for k
            
            % now do false alarms:
            % calculate where the feature value exceeds the threshold
            aboveThreshold = obj.Value > obj.Threshold;
            
            % calculate where the feature value crosses the threshold
            % thresholdCrossed is true when the threshold was crossed between this and the previous sample
            thresholdCrossed = [diff(aboveThreshold) == 1; false];  
            % set thresholdCrossed to false in seizure areas since if the
            % threshold is crossed during a seizure, it should not count as
            % a false alarm
            thresholdCrossed(isSeizure) = false;
            
            % number of times the threshold was crossed
            obj.ThresholdCrossings = sum(thresholdCrossed);
            % indices where the threshold was crossed
            obj.ThresholdCrossedLocation = find(thresholdCrossed);
            
            % calculate false alarm locations
            % one false alarm can only occur once in a certain amount of
            % time (see note in class description)
            falseAlarmLocation = obj.calculateFalseAlarmLocation(obj.ThresholdCrossedLocation); 
            
            % now we have the actual false alarms
            obj.FalseAlarms = length(falseAlarmLocation);
            obj.FalseAlarmRate = obj.FalseAlarms/(obj.Duration/3600);
            obj.FalseAlarmLocation = falseAlarmLocation;           
            
            % Calculate sensitivity
            obj.Sensitivity = sum(bitDetected)/numSeizures;
            
             % Store missed seizures location
            obj.MissedSeizures = find(~bitDetected);
            obj.BitSeizureDetected = bitDetected;
            
            % Calculate specificity
            aboveThreshold(isSeizure) = false; % only count samples that are above threshold outside of where seizures occur
            obj.Specificity = (length(aboveThreshold) - sum(aboveThreshold))/length(aboveThreshold);
  
            % Print result to command line
            if printResults
                fprintf('*** Results for %s:\n', obj.FeatureName)
                fprintf('Sensitivity: %3.2f%%\n', obj.Sensitivity*100)
                fprintf('Specificity: %3.6f%%\n', obj.Specificity*100)
                fprintf('Number of false alarms: %d\n', obj.FalseAlarms)
                fprintf('False alarm rate: %.6f\n', obj.FalseAlarmRate)
            end % if printResults
        end % function evaluateFeature
        
        function calculateBaseline(obj)
            % calculate the baseline
            
            baseline = zeros(obj.Length, 1);
            refreshRate = obj.BaselineRefreshrate;
            baselineShift = obj.BaselineShift;
            baselineWindowLength = obj.BaselineWindowLength;
            value = obj.Value;
            
            for k = 1 + baselineShift + baselineWindowLength:refreshRate:obj.Length
                baseline(k-refreshRate:k) = mean(value(...
                    k - baselineShift - baselineWindowLength + 1 : k - baselineShift));
            end % for k
            
            baseline(1:baselineShift + baselineWindowLength) = baseline(baselineShift + baselineWindowLength+1);
            baseline(end-refreshRate:end) = baseline(end - refreshRate - 1);
            
            obj.Baseline = baseline;
            
        end % function calculateBaseline
        
        function findBestThreshold(obj)
            % find the best threshold through cost function minimization
            % first, the operating characteristic curve is calculated, it
            % contains the sensitivity and false alarm rates for different
            % threshold values (sweep of all possible values)
            
            costSensitivity = obj.CostSensitivity;
            costFalseAlarmRate = obj.CostFalseAlarmRate;
            
            if isempty(obj.CharacteristicCurve)
                warning('Operating Characteristic Curve is empty. Run function calculateOperatingCharacteristicCurve first.')
            end % if
            
            sensitivity = obj.CharacteristicCurve.Sensitivity;
            falseAlarmRate = obj.CharacteristicCurve.FalseAlarmRate;
            
            % calculate the cost
            cost = costSensitivity*(1-sensitivity).^2 ...
                + costFalseAlarmRate*(falseAlarmRate).^2;
            
            % find the index that minimizes the cost
            [~, index] = min(cost);
            
            bestFactor = obj.CharacteristicCurve.ThresholdFactor(index);
            
            obj.ThresholdBaselineFactor = bestFactor;
            obj.evaluate('Print')
            
        end % function findBestThreshold
        
        function plot(obj, varargin)
            % plot the feature value
            % possible input arguments are
            % 'Measurement': also plot the linked measurement in an
            %       additional subplot
            % 'Threshold': also plot the threshold
            % 'Baseline': also plot the baseline
            % 'Mark': mark missed seizures and false alarms
            % 'Bare': removes legends and ticks and other annotations
            
            plotTimer = tic;
            
            disp(['Plotting ' obj.FeatureName '...'])
            
            % numbers of subplots
            numSubplots = 1;
            
            % check if measurement data should be plotted
            plotMeasurement = false;
            if nargin > 1
                if any(strcmp(varargin, 'Measurement'))
                    plotMeasurement = true;
                    numSubplots = numSubplots + 1;
                end
            end % if nargin > 1
            
            % create figure
            h = figure( 'units', 'pixels', ...
                'Name', obj.FeatureName, ...
                'Position', get(0,'ScreenSize'), ... % make it fullscreen
                'Tag', obj.FeatureName, ...
                'Visible', 'off'); % plotting while not visible is faster
            
            % create subplots
            for k = 1:numSubplots
                ax(k) = subplot(numSubplots, 1, k); %#ok<AGROW>
            end % for k
            
            obj.AxesHandle = ax;
            obj.PlotHandle = h;
            
            index = 1;
            % plot measurement
            if plotMeasurement
                obj.plotMeasurement(ax(index), varargin{:})
                index = index + 1;
            end % if plot Measurement
            
            % Plot feature
            obj.plotFeature(ax(index), varargin{:})
            
            % format the plot
            formatplot(h);
            
            % label xaxes of last subplot
            xlabel(ax(index), 'Time [s]')
            
            % make plot visible
            obj.PlotHandle.Visible = 'on';
            
            disp(['Plotting was done in ' num2str(toc(plotTimer), '%.2f') ' seconds.'])
            
        end % function plot
        
        function plotFeature(obj, ax, varargin)
            % plot only the feature value in the axes specified by input
            % argument ax. Possible optional input arguments are:
            % 'Threshold': also plot the threshold
            % 'Baseline': also plot the baseline
            % 'Mark': mark missed seizures and false alarms
            % 'Bare': removes legends and ticks and other annotations
            
            % check input argument
            plotThreshold = false;
            plotBaseline = false;
            bitBare = false;
            bitMark = false;
            if nargin > 2
                if any(strcmp(varargin, 'Threshold'))
                    plotThreshold = true;
                end
                if any(strcmp(varargin, 'Baseline'))
                    plotBaseline = true;
                end
                if any(strcmp(varargin, 'Bare'))
                    bitBare = true;
                end
                if any(strcmp(varargin, 'Mark'))
                    bitMark = true;
                end
            end % if nargin > 2
            
            % plot feature value, with seizures being in a different color
            [l(1), l(2)] = obj.plotSeizureHighlight(ax);
            legendstr = {obj.FeatureName, 'Seizure'};
            
            hold(ax, 'on')
            % plot Baseline and Threshold if requested
            Ts = 1/obj.MeasurementFs;
            time_feature = 0:obj.StepSize*Ts:(obj.Length-1)*obj.StepSize*Ts;
            
            if plotBaseline
                l(end+1) = plot(ax, time_feature, obj.Baseline, 'Color', [0.47 0.67 0.19]);
                legendstr{end + 1} = 'Baseline';
            end % if plotBaseline
            if plotThreshold
                l(end+1) = plot(ax, time_feature, obj.Threshold, 'Color', [0.93 0.69 0.13]);
                legendstr{end+1} = 'Threshold';
                uistack(l(end), 'bottom')
            end % if plotThreshold
            hold(ax, 'off')
            
            % mark missed seizures and false alarms if requested
            if bitMark
                ltemp = obj.markFalseAlarms(ax);
                if ~isempty(ltemp)
                    l(end+1) = ltemp;
                    legendstr{end+1} = 'False Alarms';
                end % if
                
                ltemp = obj.markMissedSeizures(ax);
                if ~isempty(ltemp)
                    l(end+1) = ltemp;
                    legendstr{end+1} = 'Missed Seizures';
                end % if
            end % if
            
            if ~bitBare
                legend(l, legendstr)
                ylabel(ax, obj.FeatureName)
            else
                ax.YTickLabel = [];
            end % if ~bitBare
        end % function plotFeature
        
        function calculateOperatingCharacteristicCurve(obj)
            % calculate the receiving operating characteristive curve of
            % this feature by sweeping through multiple threshold values
            % and calculating the sensitivity and false alarm rate for each
            % of them
            
            timeElapsed = tic;
            fprintf('Calculating operating characteristics...\n')
            
            % smaller stepsize for lower thresholdFactorValues because
            % there are more changes in that area
            thresholdFactorValues = [1:0.1:10 10.2:0.2:100];
            
            sensitivity = zeros(size(thresholdFactorValues));
            falseAlarmRate = zeros(size(thresholdFactorValues));
            
            % sweep through thresholds
            k = 1;
            while (1)
                obj.ThresholdBaselineFactor = thresholdFactorValues(k);
                obj.evaluate;
                sensitivity(k) = obj.Sensitivity;
                falseAlarmRate(k) = obj.FalseAlarmRate;
                
                if (sensitivity(k) == 0 && falseAlarmRate(k) == 0 ...
                        && sensitivity(k-1) == 0 && falseAlarmRate(k-1) == 0) ...
                        || k == length(sensitivity)
                        
                    % stop when sensitivity and false alarm rate are 0 for
                    % more than two steps or if thresholds have been
                    % calculated                   
                    break;
                else
                    k = k + 1;
                end % if
            end % while
            
            sensitivity = sensitivity(1:k);
            falseAlarmRate = falseAlarmRate(1:k);
            thresholdFactorValues = thresholdFactorValues(1:k);
            
            % create OperatingCharacteristicCurve object with results
            obj.CharacteristicCurve = OperatingCharacteristicCurve(obj.FeatureName, ...
                sensitivity, falseAlarmRate, thresholdFactorValues);
            
            fprintf('%s operating characteristic calculated in %.2f seconds.\n', obj.FeatureName, toc(timeElapsed));
            
        end % function calculateOperatingCharacteristicCurve
        
        function isSeizure = getBitIsSeizure(obj)
            % return a logical array of the same length as the value 
            % that is true at index i if that index corresponds to a 
            % seizure
            
            isSeizure = false(obj.Length, 1);
            for k = 1:obj.NumSeizures
                index = obj.SeizureStart(k):obj.SeizureEnd(k);
                isSeizure(index) = true;
            end % for k
            
        end % function
        
        function bool = isseizure(obj, index)
            % return true is input index corresponds to a seizure index and
            % false otherwise
            
            if index < 1 || index > obj.Length
                error('Index %d out of bound. Min index is 1 and max index is %d.', index, obj.Length)
            end
            
            bool = false;
            for k = 1:obj.NumSeizures
                if index >= obj.SeizureStart(k) && index <= obj.SeizureEnd(k)
                    bool = true;
                    break
                end
            end % for k
        end
        
        function bitAboveThreshold = BitAboveThreshold(obj, index)
            % return a logical vector where index i is true if the feature
            % value at index i is above the threshold
            % second input parameter is optional
            % if no index is provided, return full vector, otherwise return
            % only the specified index/indices
            
            if nargin < 2
                bitAboveThreshold = obj.Value > obj.Threshold;
            else
                if index(1) < 1 || index(end) > obj.Length
                    error('Invalid index, index should be between 1 and %d.', obj.Length)
                else
                    % input is good
                    bitAboveThreshold = obj.Value(index) > obj.Threshold(index);
                end % if
            end % if nargin
        end % function BitAboveThreshold
        
        function threshold = Threshold(obj, index)
            % return the threshold as baseline*constant
            % second input parameter is optional
            % if no index is provided, return full vector, otherwise return
            % only the specified index/indices
            % The threshold is not stored as a object property in order to
            % reduce file size. Since it can be computed very easily, this
            % should not have a big impact on speed
            
            if nargin < 2
                threshold = obj.Baseline.*obj.ThresholdBaselineFactor;
            else
                if index(1) < 1 || index(end) > obj.Length
                    error('Invalid index, index should be between 1 and %d.', obj.Length)
                else
                    % input is good
                    threshold = obj.Baseline(index).*obj.ThresholdBaselineFactor;
                end % if
            end % if nargin
            
        end % function Threshold
        
        function cut(obj, indices)
            % cuts the feature. Removes all values and seizures that are
            % not contained in the input argument indices from the feature 
            % and adjusts object properties like seizureStart and 
            % seizureEnd accordingly. 
            
            validateattributes(indices, {'numeric'}, {'vector'}, mfilename, 'indices', 1)
            
            % do a bunch of checks to see if cutting is valid (not cutting
            % through a seizure)
            if min(indices) < 1 || max(indices) > obj.Length
                error('Invalid indices.')
            end % if
            
            containsStart = ismember(obj.SeizureStart, indices);
            containsEnd = ismember(obj.SeizureEnd, indices);
            
            if any(containsStart ~= containsEnd)
                error('A seizure is getting cut. Please do not do this. Nothing is getting cut.')
            end

            % do the  cutting
            bitIsSeizure = obj.getBitIsSeizure();
            bitIsSeizure = bitIsSeizure(indices);
            bitIsSeizure = diff(bitIsSeizure);
            
            obj.SeizureStart = find(bitIsSeizure == 1);
            obj.SeizureEnd = find(bitIsSeizure == -1);      
            obj.NumSeizures = length(obj.SeizureStart);
            
            % cut values
            obj.Baseline = obj.Baseline(indices);
            obj.Value = obj.Value(indices);
            obj.Length = length(obj.Value);
            obj.Duration = obj.Length*(obj.StepSize/obj.MeasurementFs);
            
            % reset stuff that is not valid anymore
            obj.LinkedMeasurement = [];
            obj.CharacteristicCurve = [];
            obj.Sensitivity = [];
            obj.Specificity = [];
            obj.FalseAlarms = [];
            obj.ThresholdCrossedLocation = [];
            obj.FalseAlarmRate = [];
            obj.MissedSeizures = [];
            obj.BitSeizureDetected = [];
            obj.ThresholdCrossedLocation = [];
            obj.FalseAlarmLocation = [];

        end % function cut
    end % methods
    
    methods (Access = protected)
        
        function [l1, l2] = plotSeizureHighlight(obj, ax)
            % plot feature value with the parts where a seizure 
            % occurs in a different color
            
            Ts = 1/obj.MeasurementFs;
            time_feature = 0:obj.StepSize*Ts:(obj.Length-1)*obj.StepSize*Ts;
            
            % plot samples where no seizure occured
            % replace seizures samples with NaN
            tempData = obj.Value;
            tempTime = time_feature;
            for k = 1:length(obj.SeizureStart)
                tempTime(obj.SeizureStart(k):obj.SeizureEnd(k)) = nan;
                tempData(obj.SeizureStart(k):obj.SeizureEnd(k)) = nan;
            end % for k
            
            hold(ax, 'on')
            l1 = plot(ax, tempTime, tempData);
            
            % plot seizures samples in red
            tempData = nan(size(tempData));
            tempTime = nan(size(tempTime));
            for k = 1:length(obj.SeizureStart)
                tempTime(obj.SeizureStart(k):obj.SeizureEnd(k)) = ...
                    time_feature(obj.SeizureStart(k):obj.SeizureEnd(k));
                tempData(obj.SeizureStart(k):obj.SeizureEnd(k)) = ...
                    obj.Value(obj.SeizureStart(k):obj.SeizureEnd(k));
            end % for k
            
            l2 = plot(ax, tempTime, tempData, 'r');
            hold(ax, 'off')
            
        end % function seizureHighlighPlot
        
        function plotMeasurement(obj, ax, varargin)
            % plots the measurement
            % indizes where seizure occurs are in a different color
            
            bitBare = false;
            if nargin > 3
                if any(strcmp(varargin, 'Bare'))
                    bitBare = true;
                end % if any
            end % if
            
            measurement = obj.LinkedMeasurement;
            
            tempData = measurement.Data;
            time_measurement = 0:measurement.Ts:(measurement.Length-1)*measurement.Ts;
            tempTime = time_measurement;
            
            for k = 1:length(measurement.SeizureStart)
                tempTime(measurement.SeizureStart(k):measurement.SeizureEnd(k)) = nan;
                tempData(measurement.SeizureStart(k):measurement.SeizureEnd(k)) = nan;
            end % for k
            l = plot(ax, tempTime, tempData);
            legendstr = {'iEEG Data'};
            
            if ~isempty(measurement.SeizureStart)
                % plot seizures in a different color
                tempData = nan(size(measurement.Data));
                tempTime = nan(size(time_measurement));
                
                for k = 1:length(measurement.SeizureStart)
                    tempTime(measurement.SeizureStart(k):measurement.SeizureEnd(k)) = ...
                        time_measurement(measurement.SeizureStart(k):measurement.SeizureEnd(k));
                    tempData(measurement.SeizureStart(k):measurement.SeizureEnd(k)) = ...
                        measurement.Data(measurement.SeizureStart(k):measurement.SeizureEnd(k));
                end % for k
                hold(ax, 'on')
                l(end+1) = plot(ax, tempTime, tempData, 'r');
                hold(ax, 'off')
                legendstr{end+1} = 'Seizure';
            end % if ~isempty
            
            if ~bitBare
                legend(l, legendstr);
                ylabel('EEG Data')
            end % if
        end % function plotMeasurement
        
        function l = markFalseAlarms(obj, ax)
            % mark false alarms on plot
            
            if isempty(obj.FalseAlarms)
                warning([obj.FeatureName ' has no false alarms.'])
                l = [];
            else
                Ts = 1/obj.MeasurementFs;
                time_feature = 0:obj.StepSize*Ts:(obj.Length-1)*obj.StepSize*Ts;
                value = obj.Threshold;
                
                hold(ax, 'on')
                l = plot(ax, time_feature(obj.FalseAlarmLocation), ...
                    value(obj.FalseAlarmLocation), 'kx', 'MarkerSize', 10);
                hold(ax, 'off')
            end
        end % function markFalseAlarms
        
        function l = markMissedSeizures(obj, ax)
            % mark missed seizures on plot
            
            if isempty(obj.MissedSeizures)
                warning([obj.FeatureName ' has no missed seizures.'])
                l = [];
            else
                % plot them
                Ts = 1/obj.MeasurementFs;
                time_feature = 0:obj.StepSize*Ts:(obj.Length-1)*obj.StepSize*Ts;
                indexMissed = obj.SeizureStart(obj.MissedSeizures);
                
                hold(ax, 'on')
                l = plot(ax, time_feature(indexMissed), obj.Value(indexMissed), 'dk', 'MarkerSize', 10);
                hold(ax, 'off')
            end % if
        end % function markMissedSeizures
        
        function falseAlarmLocation = calculateFalseAlarmLocation(obj, thresholdCrossedLocation)
            % calculate location of false alarms. A false alarm can only
            % occur a certain amount of samples after a previous false
            % alarm, not before
            
            if ~isempty(thresholdCrossedLocation)
                % only do it if the threshold was crossed at least once 
                
                % if multiple threshold crossings occur within a short time
                % (specified by seizureHoldTime) of each other, only count
                % one false alarm (the first one)
                
                loc = thresholdCrossedLocation;
                indexIgnore = false(size(thresholdCrossedLocation));
                k = 1;
                while ~isempty(k)
                    
                    loc = loc - loc(k);
                    indexIgnore = indexIgnore | (loc > 0 & loc <= obj.SeizureHoldTime);
                    
                    k = find(loc > obj.SeizureHoldTime, 1, 'first');
                    
                end % while k
                
                falseAlarmLocation = thresholdCrossedLocation(~indexIgnore);
            else
                % no false alarms
                falseAlarmLocation = [];                
            end % if ~isempty
        end % function calculateFalseAlarmLocation
    end % methods (Access = protected)
    
    methods (Abstract, Access = protected)
        calculateFeatureValue(obj)
    end % methods (Abstract, Access = protected)
end

