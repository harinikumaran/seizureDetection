classdef LineLength < FeatureCalculation
    % LINELENGTH is a class to analysis the line length feature
    %
    % The detection threshold is given by: k*baseline + c
    
    properties
        FeatureName = 'Line Length'
    end % properties
    
    methods 
        % constructor
        function obj = LineLength(measurement, varargin)
            obj = obj@FeatureCalculation(measurement, varargin{:});
        end % function LineLength
    end % methods
    
    methods (Access = protected)
        function calculateFeatureValue(obj)
            
            measurement = obj.LinkedMeasurement;
            data = measurement.Data;
            data(measurement.isNaN) = NaN; % replace original NaNs with NaNs
            
            temp = movsum(abs(diff(data)), [obj.WindowLength-1 0])/obj.WindowLength;
            obj.Value = temp(1:obj.StepSize:end);
            
        end % function calculateFeatureValue
    end % methods (Access = protected)
    
end % classdef

