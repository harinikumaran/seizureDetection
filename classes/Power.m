classdef Power < FeatureCalculation
    
    properties
        FeatureName = 'Power'
    end
    
    methods
        % constructor
        function obj = Power(measurement, varargin)
            obj = obj@FeatureCalculation(measurement, varargin{:});
        end % function Power
    end % methods
    methods (Access = protected)
        function calculateFeatureValue(obj)
            
            measurement = obj.LinkedMeasurement;
            data = measurement.Data;
            data(measurement.isNaN) = NaN; % replace original NaNs with NaNs
            
            temp = movmean(data.^2, [obj.WindowLength-1 0]);            
            obj.Value = temp(1:obj.StepSize:end);
            
        end % function calculateFeatureValue
    end % methods (Access = protected)
end

