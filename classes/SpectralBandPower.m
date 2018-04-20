classdef SpectralBandPower < FeatureCalculation
    % SPECTRALBANDPOWER is an object to calculate the spectral band power
    % feature
    %
    % Contrary to the other FeatureCalculation children, it has two more
    % required input parameters in addition to the measurement object:
    %   -bandName: the name of the frequency band
    %   -frequencyBounds: the frequency boundaries of the band
    %
    % Example:
    %   thetaBand = SpectralBandPower(measurement, 'ThetaBand', [4 8]);
    
    properties
        FeatureName
        FrequencyBounds % [lowerLimit upperLimit] in Hz
        
        FilterOrder
    end
    
    methods
        % constructor
        function obj = SpectralBandPower(measurement, bandName, frequencyBounds,  varargin)
            obj = obj@FeatureCalculation(measurement, varargin{:});

            % parse input
            p = inputParser;
            p.KeepUnmatched = true;
            
            % define default values            
            defaultFilterOrder = 6;
            
            addRequired(p, 'bandName', @(x) validateattributes(x,{'char'},{'scalartext'},mfilename,'bandName',2));
            addRequired(p, 'frequencyBounds', @(x) validateattributes(x,{'numeric'},{'numel', 2},mfilename,'frequencyBounds',3));
            addParameter(p, 'FilterOrder', defaultFilterOrder, @isscalar);
            
            parse(p, bandName, frequencyBounds, varargin{:});
            
            obj.FeatureName = p.Results.bandName;
            obj.FrequencyBounds = p.Results.frequencyBounds;
            obj.FilterOrder = p.Results.FilterOrder;
        end % constructor
    end % methods
    
    methods (Access = protected)
        function calculateFeatureValue(obj)
            
            measurement = obj.LinkedMeasurement;
            Fs = measurement.Fs;
            
            % filter the data using a 10th order butterworth filter
            N   = obj.FilterOrder;   % Order
            Fc1 = obj.FrequencyBounds(1);  % First Cutoff Frequency
            Fc2 = obj.FrequencyBounds(2);  % Second Cutoff Frequency
            
            h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
            Hd = design(h, 'butter');
            
            dataFiltered = filter(Hd, measurement.Data);

            dataFiltered(measurement.isNaN) = NaN; % replace original NaNs with NaNs
            
            temp = movmean(dataFiltered.^2, [obj.WindowLength-1 0]);            
            obj.Value = temp(1:obj.StepSize:end);
                     
        end % function calculateFeatureValue
    end % methods
end

