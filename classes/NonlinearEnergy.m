classdef NonlinearEnergy < FeatureCalculation
    %NONLINEARENERGY is a class to calculate and evaluate the
    %nonlinear energy feature
    
    properties
        FeatureName = 'Nonlinear Energy'
    end
    
    methods
        % constructor
        function obj = NonlinearEnergy(measurement, varargin)
            obj = obj@FeatureCalculation(measurement, varargin{:});
        end % function FeatureNonlinearEnergy        
    end % methods
    
    methods (Access = protected)
        function calculateFeatureValue(obj)
            
            measurement = obj.LinkedMeasurement;
            data = measurement.Data;
            data(measurement.isNaN) = NaN; % replace original NaNs with NaNs
             
            val = data.^2 - circshift(data, 1).*circshift(data, -1);
            val([1 end]) = 0;
            temp = movmean(val, [obj.WindowLength-1 0]);
            
            obj.Value = temp(1:obj.StepSize:end);
            
        end % function calculateFeatureValue
    end % methods (Access = protected)
end

