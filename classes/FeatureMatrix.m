classdef FeatureMatrix < matlab.mixin.Copyable
    % FEATUREMATRIX is an object that hold a matrix of features. For every
    % measurement, there are a certain number of channels and features
    % (line length, power, etc...). This gives a gives a #channels x
    % #features number of different signals to evaluate. Some of these
    % signals (from now on called features, i.e one feature is for example
    % line length of channel 2) give better results than other.
    %
    % All features are stored inside the Matrix property. The elements of
    % the matrix have to be set after initialization through the setElement
    % setRow or addRow function.
    %
    % The constructor has 2 required input parameters, numChannels and
    % featureNames, numChannel being of type double and
    % featureNames being a cell array of feature name string. These strings
    % have to match the name of the feature objects exactly.
    
    properties
        
    end
    
    properties (SetAccess = protected) % Read-Only properties
        
        NumChannels % Number of channels
        NumFeatures % Number of features
        FeatureNames % Cell array of feature names in order in which they appear
        FeatureLength % Length of one feature
        
        NumSeizures % Number of seizures per feature (same for all elements of the matrix since they based on the same measurement)
        NumElements % Number of elements in the matrix
        
        Matrix % Cell matrix of feature objects
        SensitivityMatrix % Sensitivity of each feature
        SpecificityMatrix % Specificity of each feature
        FalseAlarmRateMatrix % FAR for each feature
        AreaMatrix % Area under false alarm rate sensitivity curve
        
        FlagMatrix % 3D array of vectors
        BitSeizureDetectedMatrix %[numChannels x numFeatures x numSeizures], index i,j,k is true when feature j on channel i detect seizure k
        
    end %  properties (SetAccess = protected) % Read-Only properties
    
    methods
        % constructor
        function obj = FeatureMatrix(numChannels, featureNames)
            % constructor has 2 required input parameters, numChannels and
            % featureNames (see class description for more detail)
            
            obj.FeatureNames = featureNames;
            obj.NumFeatures = length(featureNames);
            obj.NumChannels = numChannels;
            obj.Matrix = cell(obj.NumChannels, obj.NumFeatures);
            obj.NumElements = numel(obj.Matrix);
        end % constructor
        
        function setElement(obj, indexChannel, feature)
            % set a matrix element, the feature object name has to match
            % one of the strings inside the FeatureNames property of the
            % matrix
            
            % check input
            if ~isobject(feature) || ~isprop(feature, 'FeatureName')
                error('Wrong input.')
            elseif indexChannel > obj.NumChannels
                error('Channel index too big. The Matrix has %d channels', obj.NumChannels)
            end
            
            indexFeature = find(strcmp(obj.FeatureNames, feature.FeatureName));
            if isempty(indexFeature)
                error('Input feature name does not match with any columns name.')
            end % if
            
            % all elements from the matrix must come from the same
            % dataset, i.e. the same amount of seizures
            if isempty(obj.NumSeizures)
                obj.NumSeizures = feature.NumSeizures;
                obj.FeatureLength = feature.Length;
            elseif obj.NumSeizures ~= feature.NumSeizures
                error('The number of seizures in the input does not match current matrix elements.')
            elseif obj.FeatureLength ~= feature.Length
                error('The input feature has different length than the other elements of the matrix.')
            end % if
            
            % all good, save feature in matrix
            feature.LinkedMeasurement = [];
            obj.Matrix{indexChannel, indexFeature} = copy(feature);
            
        end % function setElement
        
        function fillAreaMatrix(obj)
            % fill AreaMatrix property
            
            % upperBound for the false alarm rate when calculating the area under curve
            % area under the false alarm rate-sensitivity curve will the
            % calculate betwenn false alarm rate of 0 and upperBound
            upperBound = 1;
            
            obj.AreaMatrix = zeros(size(obj.Matrix));
            
            % fill matrices with results from every feature
            for row = 1:obj.NumChannels
                for col = 1:obj.NumFeatures
                    try
                        obj.AreaMatrix(row, col) = obj.Matrix{row, col}.CharacteristicCurve.calculateArea(upperBound);
                    catch
                        warning('Could not calculate area under curve of feature %s.', obj.Matrix{row, col}.Name)
                    end % try
                end % for col
            end % for row
        end % function fillAreaMatrix(obj)
        
        function valueMatrix = getValueMatrix(obj, dimString)
            % return a matrix containing the values of the single feature
            % optional input argument dimString can be set to '2D' and the
            % returned matrix will have dimensions [featureLengh x
            % numTotalFeatures], otherwise it will be [numChannel X
            % numFeaturesPerChannel x  featureLength]
            
            
            % 3D Matrix [channel x feature x sample]
            valueMatrix = zeros(obj.NumChannels, obj.NumFeatures, obj.FeatureLength);
            
            % fill matrices with results from every feature
            for row = 1:obj.NumChannels
                for col = 1:obj.NumFeatures
                    
                    % above threshold matrix
                    valueMatrix(row, col, :) = obj.Matrix{row, col}.Value;
                    
                end % for col
            end % for row
            
            if nargin == 2 && strcmp(dimString, '2D')
                valueMatrix = transpose(reshape(valueMatrix, obj.NumElements, obj.FeatureLength));
            end
            
        end % function getFlagMatrix
        
        function fillMatrices(obj)
            % fill the SensitivityMatrix, FalseAlarmRateMatrix and
            % SpecificityMatrix properties
            
            obj.SensitivityMatrix = zeros(size(obj.Matrix));
            obj.FalseAlarmRateMatrix = zeros(size(obj.Matrix));
            obj.SpecificityMatrix = zeros(size(obj.Matrix));
            
            % 3D Matrices [channel x feature x sample]
            obj.FlagMatrix = false(obj.NumChannels, obj.NumFeatures, obj.FeatureLength);
            obj.BitSeizureDetectedMatrix = true(obj.NumChannels, obj.NumFeatures, obj.NumSeizures);
            
            % fill matrices with results from every feature
            for row = 1:obj.NumChannels
                for col = 1:obj.NumFeatures
                    
                    obj.SensitivityMatrix(row, col) = obj.Matrix{row, col}.Sensitivity;
                    obj.SpecificityMatrix(row, col) = obj.Matrix{row, col}.Specificity;
                    obj.FalseAlarmRateMatrix(row, col) = obj.Matrix{row, col}.FalseAlarmRate;
                    
                    % above threshold matrix
                    obj.FlagMatrix(row, col, :) = obj.Matrix{row, col}.BitAboveThreshold;
                    
                    % missed seizures matrix
                    obj.BitSeizureDetectedMatrix(row, col, obj.Matrix{row, col}.MissedSeizures) = false;
                    
                end % for col
            end % for row
        end % function fillMatrices
        
        function feature = bestFeature(obj, varargin)
            % return the best feature or if an index k is provided,
            % the kth best feature
            % the 'best' feature is determined by the area under curve
            
            if nargin == 2
                k = varargin;
            else
                k = 1;
            end
            
            if isempty(obj.AreaMatrix)
                obj.fillAreaMatrix();
            end % if
            
            [~, I] = sort(obj.AreaMatrix(:), 'descend'); % sort areas from best to worst
            
            featureIndex = I(k); % get index of kth best area
            
            feature = obj.Matrix{featureIndex}; % return that feature
        end % function bestFeature
        
        function plot(obj, strWhat, index, varargin)
            %   plot('feature', 2, 'Threshold', 'Bare', 'Mark')
            %   plot('channel', 13, 'Threshold', 'Mark')
            %
            % plot features, either plot all features for one channel or
            % all channels for one feature
            % first input argument (besides obj) should be either 'channel'
            % or 'feature', second should be the index of the
            % feature/channel to plot. Further input arguments can be
            % entered to costumize how the plot looks. Possibilities:
            % 'Threshold', 'Bare', 'Mark', 'Baseline'
            
            if strcmp(strWhat, 'channel')
                if index > obj.NumChannels
                    error('Channel index out of bounds.')
                end % if
                indexToPlot = index:obj.NumChannels:obj.NumElements;
            elseif strcmp(strWhat, 'feature')
                if index > obj.NumFeatures
                    error('Feature index out of bounds.')
                end % if
                startIndex = 1+(obj.NumChannels*(index-1));
                indexToPlot = startIndex:1:startIndex+obj.NumChannels-1;
            else
                error('Invalid input. First input argument should either be "channel" or "feature".')
            end % if
            
            bitBare = false;
            if nargin > 3
                if any(strcmp(varargin, 'Bare'))
                    bitBare = true;
                end % if
            end % if
            
            h = figure( 'units', 'pixels', ...
                'Name', 'Feature Matrix', ...
                'Position', get(0,'ScreenSize'), ... % make it fullscreen
                'Tag', 'Feature Matrix', ...
                'Visible', 'off'); % plotting while not visible is faster
            
            numSubplots = length(indexToPlot);
            
            % call plot function from feature object
            for k = 1:numSubplots
                ax(k) = subplot(numSubplots, 1, k); %#ok<AGROW>
                obj.Matrix{indexToPlot(k)}.plotFeature(ax(k), varargin{:});
            end % for
            
            % format the plot
            if bitBare
                formatplot(h, 'TopMargin', 0.01, 'LeftMargin', 0.01, 'BottomMargin', 0.01, 'RightMargin', 0.01, 'Gap', 0.005);
                ax(end).XTickLabel = [];
            else
                formatplot(h, 'TopMargin', 0.01, 'Gap', 0.002, 'LegendOutside', true);
                xlabel(ax(end), 'Time [s]');
            end
            
            h.Visible = 'on';
            
        end % function plot
        
        function deleteRows(obj, indexRow)
            % delete rows with index indexRow
            
            if indexRow(1) < 1 || indexRow(end) > obj.NumChannels
                error('Invalid row index.')
            end % if
            
            obj.Matrix(indexRow, :) = [];
            obj.NumChannels = size(obj.Matrix, 1);
            obj.NumElements = numel(obj.Matrix);
            
        end % function deleteRow
        
        function deleteFeature(obj, nameFeature)
            % delete a feature with name nameFeature, the whole column will
            % be deleted!
            
            indexCol = find(strcmp(obj.FeatureNames, nameFeature));
            if isempty(indexCol)
                error('Invalid feature name.')
            end % if
            
            obj.Matrix(:, indexCol) = [];
            obj.FeatureNames(indexCol) = [];
            obj.NumFeatures = length(obj.FeatureNames);
            obj.NumElements = numel(obj.Matrix);
            
            fprint('Succefully deleted feature %s.\n', nameFeature)
            
        end % function deleteRow
        
        function addRow(obj, insertCell, indexInsert)
            % adds a row to the matrix. indexInsert parameter is optional,
            % if it is not specified insert it at the end
            
            % check input
            if ~all(size(insertCell) == [1, obj.NumFeatures])
                error('The row to insert should be a 1 x %d cell array', obj.NumFeatures)
            elseif ~all(strcmp(obj.FeatureNames, cellfun(@(x) x.FeatureName, insertCell, 'UniformOutput', false)))
                error('Input does not match feature names and/or order in Matrix.')
            end % if
            
            if nargin == 2
                % insert at end
                
                obj.Matrix = [obj.Matrix; insertCell];
                indexInsert = obj.NumChannels + 1;
            else
                % insert it at indexRow
                
                if indexInsert > obj.NumChannels + 1
                    error('indexRow insert is too big.')
                end % if
                
                obj.Matrix = [obj.Matrix(1:indexInsert-1, :); insertCell; obj.Matrix(indexInsert:end, :)];
            end % if
            
            fprintf('Row succesfully inserted at row index %d\n', indexInsert)
            
            obj.NumChannels = size(obj.Matrix, 1);
            obj.NumElements = numel(obj.Matrix);
            
        end % function addRow
        
        function setRow(obj, insertCell, index)
            % set a row of the matrix
            
            % check input
            if ~all(size(insertCell) == [1, obj.NumFeatures])
                error('The row to insert should be a 1 x %d cell array', obj.NumFeatures)
            elseif ~all(strcmp(obj.FeatureNames, cellfun(@(x) x.FeatureName, insertCell, 'UniformOutput', false)))
                error('Input does not match feature names and/or order in Matrix.')
            end % if
            
            
            obj.Matrix(index, :) = insertCell;
            
            fprintf('Row %d succesfully set.\n', index)
            
        end % function setRow
        
        function cut(obj, indices)
            % removes data from the matrix, only keeps indices specified by
            % input argument
            % this function removes all indices not contained in the input
            % argument indices from all features inside the matrix and
            % adjusts the properties of the individual features and of the
            % matrix
            
            % call cut function in every feature
            for k = 1:obj.NumElements
                obj.Matrix{k}.cut(indices);
            end % for k
            
            obj.FeatureLength = obj.Matrix{1}.Length;
            obj.NumSeizures = obj.Matrix{1}.NumSeizures;
            
            % delete all matrices since they are not valid anymore
            obj.SensitivityMatrix = [];
            obj.SpecificityMatrix = [];
            obj.FalseAlarmRateMatrix = [];
            obj.AreaMatrix = [];
            obj.FlagMatrix = [];
            obj.BitSeizureDetectedMatrix = [];   
        end % function cut
        
        function calculateFeatureCharacteristicCurve(obj)
            % calcute the operating characteristic curve of all features
            
            % create a wait bar because this could take a while
            totalTimer = tic();
            wb = waitbar(0, 'Starting calculations...', 'Name', 'Calculating...');
            
            for k = 1:obj.NumElements
                obj.Matrix{k}.calculateOperatingCharacteristicCurve;
                waitbar(k/obj.NumElements, wb, [num2str(k/obj.NumElements*100, '%.2f') '%, ' num2str(toc(totalTimer)/60, '%.2f') ' min. elapsed'])
            end % for k
            
            close(wb)
            fprintf('Calculated curves of all features in %.2f seconds.\n', toc(totalTimer))
            
        end % function calculateElementCharacteristicCurve
        
    end % methods
    
    methods (Access = protected)
        function cpObj = copyElement(obj)
            disp('Copying feature matrix.')
            
            cpObj = copyElement@matlab.mixin.Copyable(obj); % shallow copy of all properties
            
            for k = 1:obj.NumElements
                cpObj.Matrix{k} = copy(obj.Matrix{k});
            end % for
        end % function copyElement
    end % methods (Access = protected)
end

