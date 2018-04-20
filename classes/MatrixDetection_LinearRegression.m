classdef MatrixDetection_LinearRegression < MatrixDetection
    % Calculate the weights through linear regression
    %
    % The CalculationFunctionHandle property should be a function handle to
    % calculate the weights
    % For ridge regression with a small regularisation parameter this could
    % for example be
    %   @(X, y) (X'*X + 0.01*eye(size(X'*X)))\X'*y
    %
    % Note: for the purpose of calculating the weights, every seizure is
    % compacted to a single value which is either 1 if the feature detected
    % that seizure and 0 otherwise
    
    
    
    properties
        WeightSeizure % Datapoint weight for seizure datapoints
        WeightInterictal % Datapoint weight for interictal datapoint
        
        CalculationFunctionHandle % function handle
    end
    
    methods
        function obj = MatrixDetection_LinearRegression(varargin)
            obj = obj@MatrixDetection(varargin{:});
            % first input argument must be matrix, second input argument
            % must be function handle
            
            % parse input
            defaultWeightSeizure = obj.CostSensitivity;
            defaultWeightInterictal = obj.CostFalseAlarmRate;
            
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'WeightSeizure', defaultWeightSeizure, @isscalar);
            addParameter(p,'WeightInterictal',defaultWeightInterictal, @isscalar);

            parse(p,varargin{2:end}); % first input parameter was matrix and is required by superclass
            
            obj.WeightSeizure = p.Results.WeightSeizure;
            obj.WeightInterictal = p.Results.WeightInterictal;
            
        end % constructor   
        
        function obj = set.CalculationFunctionHandle(obj, funcHandle)
            % setter function for CalculationFunctionHandle property
            validateattributes(funcHandle,{'function_handle'}, {'nonempty'}, mfilename)
            obj.CalculationFunctionHandle = funcHandle;           
        end % function setCalculationFunctionHandle
  
        function calculateWeights(obj)
            % calculate the feature weights
            
            matrix = obj.LinkedMatrix;
            
            % get number of seizures and number of interictal samples
            nSeizures = matrix.NumSeizures;
            nIntericalSamples = sum(~matrix.Matrix{1}.getBitIsSeizure);
            
            % first remove all seizure areas from the data
            interical = double(matrix.FlagMatrix);
            interical = reshape(interical, matrix.NumElements, matrix.FeatureLength)';
            interical(matrix.Matrix{1}.getBitIsSeizure, :) = [];
            
            detected = reshape(matrix.BitSeizureDetectedMatrix, matrix.NumElements, matrix.NumSeizures)';
            
            % now put the seizure back at the end of the data. Every
            % seizure will only correspond to 1 data point which is either
            % 1 if the seizure is detect or 0 otherwise
            X = [sqrt(obj.WeightInterictal)*interical; detected*sqrt(obj.WeightSeizure)];
            y = [zeros(nIntericalSamples, 1); sqrt(obj.WeightSeizure)*ones(nSeizures, 1)]; % target vector
            
            % calculate the weights
            w = obj.CalculationFunctionHandle(X, y);        
            
            % store the weights
            obj.WeightMatrix = zeros(size(matrix.Matrix));
            obj.WeightMatrix(1:matrix.NumElements) = w;
            
        end      
    end % methods
        
end % classdef

