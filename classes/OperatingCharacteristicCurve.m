classdef OperatingCharacteristicCurve < matlab.mixin.Copyable
    %OPERATINGCHARACTERISTICCURVE hold data relating to the operating
    %characterictic curve
    %
    % The operating characteristic curve here is a curve with the false
    % alarm rate on the x-axis and the sensitivity on the y-axis
    %
    % This class also contains a function calculateArea to calculate the 
    % area under the curve 
    
    properties
        Name % Name for legend when plotting
    end
    
    properties (SetAccess = protected)
        Sensitivity % vector of achieved sensitivities (y-axis)
        FalseAlarmRate % vector ofcorresponding false alarm rate (x-axis)
        ThresholdFactor % vector of threshold factors corresponding to the obtained values
        Area % scalar area under the curve 
    end
    
    methods
        % contructor
        function obj = OperatingCharacteristicCurve(name, sensitivity, falseAlarmRate, threshold)
            
            if ~isequal(size(sensitivity), size(falseAlarmRate))
                error('Inputs do not have same size.')
            end % if
                       
            obj.Sensitivity = sensitivity(:);
            obj.FalseAlarmRate = falseAlarmRate(:);
            obj.ThresholdFactor = threshold(:);
            obj.Name = name;
            
        end % function OperatingCharacteristicCurve
        
        function [h, ax] = plot(obj, displayName)
            % plot the curve on a new plot, returns the figure and axes
            % handle
            
            if nargin == 1
                displayName = obj.Name;
            end % if
            
            h = figure( 'units', 'pixels', ...
                'Name', 'Operating Characteristic Curve', ...
                'Tag', 'Operating Characteristic Curve');
            
            ax = axes;
            
            % use stairs instead of plot to better represent the discrete
            % natur of the data
            % plot(ax, obj.FalseAlarmRate, obj.Sensitivity, 'DisplayName', displayName);
            stairs(ax, obj.FalseAlarmRate, obj.Sensitivity, 'DisplayName', displayName);
            
            xlabel('False Alarm Rate [false alarms/hour]')
            ylabel('Sensitivity')
            legend('location', 'southeast')
            title(ax, 'Sensitivity vs. false alarm rate plot')
            
            formatplot(h);
            ax.XLim = [0 5];
            
        end % function plot
        
        function addCurve(obj, handle, displayName)
            % adds a curve to an already existing plot
            % - second input parameter 'displayName' is optional
            % - first input parameter 'handle' can be either a handle to the
            % axis or a handle to the figure to which the new curve is to
            % be added
            
            if isa(handle, 'matlab.graphics.axis.Axes')
                % input handle is axis object
                
                ax = handle;
                
            elseif isa(handle, 'matlab.ui.Figure')
                % input handle is figure handle
                
                % find all children of figure that are of type Axes
                ax = findall(handle, 'Type', 'Axes');
                
                if ~isscalar(ax)
                    % there are more than one object of type Axes in this
                    % figure, something is wrong
                    
                    error('Found multiple axes within the figure. No plotting done.')
                    
                end % if ~isscalar
                
            else
                % unknown handle
                error('Bad input.')
            end % if
            
            if nargin < 3
                displayName = obj.Name;
            end % if nargin < 3
            
            hold(ax, 'on')
            % plot(ax, obj.FalseAlarmRate, obj.Sensitivity, 'DisplayName', displayName)
            stairs(ax, obj.FalseAlarmRate, obj.Sensitivity, 'DisplayName', displayName)
            hold(ax, 'off')
            
        end % function addCurve
        
        function area = calculateArea(obj, upperBound)
            % calculate the area under the curve
            %
            % upperbound gives the upper integration limit,
            % the calculated area will be the integral from x = 0 to 
            % x = upperBound (the x-axis is the false alarm rate)
            
            if nargin == 1
                upperBound = 1;
            end
            
            index = find(obj.FalseAlarmRate <= upperBound);
            falseAlarmRate = obj.FalseAlarmRate(index);
            sensitivity = obj.Sensitivity(index);
            
            if isempty(sensitivity)
                area = 0;
                obj.Area = 0;
                return
            end % if isempty
            
            [falseAlarmRate, indexSort] = sort(falseAlarmRate);
            sensitivity = sensitivity(indexSort);
            
            area = trapz([falseAlarmRate; upperBound], ...
                [sensitivity; sensitivity(end)]);
            
            obj.Area = area;
            
        end % function calculateArea
    end
    
end

