classdef Measurement < matlab.mixin.Copyable
    % MEASUREMENT is a class used to store the EEG recording data
    %
    % Preprocesing is done through a Butterworth Bandpass filter
    %
    % The class is a children of 'matlab.mixin.Copyable' so that it can be
    % deep copied by calling the copy function
    %
    % Required input arguments for the constructor are:
    %   data: the recording data
    %   fs: sampling frequency in Hz
    %   seizureStart: indices where seizures begin (in samples)
    %   seizureEnd: indices where seizures end (in samples)
    %
    % Parameters passed as Name-Value pair to the constructor are:
    %   StudyName: the name of the study [char]
    %   StudyChannel: the channel of the measurement [double]
    %   StudyChannelName: sometimes channels have a specific name [char]
    
    
    properties
        FilterOrder = 6;
        FilterFc1 = 1; % Hz
        FilterFc2 = 70; % Hz
        
        StudyName % Name of associated study
        StudyChannelNumber % Channel number in associated study
        StudyChannelName % Name of the channel
    end
    
    properties (SetAccess = protected) % Read-Only properties
        Data % EEG or iEEG data
        Ts % Sampling time of the data
        Fs % Sampling rate/frequency of the data
        Length % Total number of samples
        NumSeizures % Total number of seizures
        SeizureStart % In terms of samples
        SeizureEnd % In terms of samples
        SeizureDuration % In seconds
        Duration % In seconds
        isNaN % Logic vector, true if the original sample was NaN
    end % properties (SetAcces = protected)
    
    methods
        % constructor
        function obj = Measurement(varargin)
            % obj = Measurement(data, Fs, seizureStart, seizureEnd, 'Name', 'Value')
            
            defaultStudyName = 'generic study';
            defaultStudyChannelNumber = -1;
            defaultStudyChannelName = '';
            
            p = inputParser;
            p.KeepUnmatched = false;
            addRequired(p, 'data', @isnumeric);
            addRequired(p, 'Fs', @isscalar);
            addRequired(p, 'seizureStart', @isnumeric);
            addRequired(p, 'seizureEnd', @isnumeric);
            addParameter(p,'StudyName', defaultStudyName, @ischar);
            addParameter(p,'StudyChannelName', defaultStudyChannelName, @ischar);
            addParameter(p,'StudyChannelNumber', defaultStudyChannelNumber, @isscalar);
            parse(p,varargin{:});
                        
            obj.Data = p.Results.data(:); % force column vector
            obj.Fs = p.Results.Fs;
            obj.Ts = 1/obj.Fs;
            obj.Length = length(obj.Data);
            obj.Duration = (obj.Length-1)*obj.Ts;
            obj.SeizureStart = p.Results.seizureStart(:);
            obj.SeizureEnd = p.Results.seizureEnd(:);
            obj.NumSeizures = length(obj.SeizureStart);

            obj.StudyName = p.Results.StudyName;            
            obj.StudyChannelNumber = p.Results.StudyChannelNumber;
            
            % calculate seizure duration
            obj.SeizureDuration = (obj.SeizureEnd - obj.SeizureStart)*obj.Ts;
            
            obj.details()
            
            checkData(obj)
        end % constructor
                
        function preprocess(obj)
            % preprocess data
            
            disp('Preprocessing data...')
            
            preprocessTimer = tic;
            
            % Create Filter            
            h  = fdesign.bandpass('N,F3dB1,F3dB2', obj.FilterOrder, obj.FilterFc1, obj.FilterFc2, obj.Fs);
            Hd = design(h, 'butter');
            
            % Filter data
            obj.Data = filter(Hd, obj.Data);

            disp(['Finished preprocessing data in ' ...
                num2str(toc(preprocessTimer), '%.2f') ' seconds.'])
            
        end % function preprocess
        
        function downsample(obj, n)
            % downsamples measurement by specified factor           
            % for example, if Fs = 500Hz and n = 2, the measurment will be
            % downsamples to Fs = 250Hz
            % factor has to be an integer
            
            if n < 2 || rem(n,1) ~= 0
                error('Invalid input. Downsampling factor has to an integer bigger or equal 2.')
            end % if
            if rem(obj.Fs, n)
                error('The resulting sampling frequency will not be an integer. Not supported.')
            end 
            
            obj.Data = decimate(obj.Data, n); % decimate filters the data first to avoid anti aliasing
            obj.isNaN = downsample(obj.isNaN, n); % no need for an anti aliasing filter here
            obj.Length = length(obj.Data); % update length
            obj.Fs = obj.Fs/n; % update sampling frequency
            obj.Ts = 1/obj.Fs; % update sampling time
            
            % update seizure start and end
            obj.SeizureStart = floor(obj.SeizureStart/n);
            obj.SeizureEnd = ceil(obj.SeizureEnd/n);
            
            fprintf('Resampled measurement. New sampling rate is %.2fHz.\n', obj.Fs)
            
        end % function resample
        
        function cut(obj, startIndex, endIndex)
            % keeps only the part of the measurement between startIndex and
            % stopIndex (parameters are in samples, not seconds)
            
            if ~isscalar(startIndex) || ~isscalar(endIndex)
                error('Input indices should be scalar values.')
            end % if
             
            % check cut indices
            if startIndex < 1 || startIndex > endIndex
                error('Invalid start index.')
            elseif endIndex > obj.Length
                error('Stop index too large.')
            end % if
            
            % print a warning if some seizures are removed
            if endIndex < obj.SeizureEnd(end) || startIndex > obj.SeizureStart(1)
                warning('Some seizures are getting removed.')
            end
            
            bitCut = obj.SeizureEnd < startIndex | obj.SeizureStart > endIndex;
            
            newStart = obj.SeizureStart;
            newStart(bitCut) = [];
            newEnd = obj.SeizureEnd;
            newEnd(bitCut) = [];
            
            newStart = newStart - startIndex + 1;
            newEnd = newEnd - startIndex + 1;
            
            newLength = endIndex - startIndex + 1;
            
            % do not cut through seizures
            if any(newStart < 0) || any(newEnd > newLength)
                error('Start index or end index within a seizure. Please do not cut seizures... Measurement was not cut.')
            end % if
            
            % all good, do the actual cutting
            obj.SeizureStart = newStart;
            obj.SeizureEnd = newEnd;            
            obj.NumSeizures = length(obj.SeizureStart);
                               
            obj.Data = obj.Data(startIndex:endIndex);
            obj.isNaN = obj.isNaN(startIndex:endIndex);
            obj.Length = newLength;
            obj.Duration = (obj.Length-1)*obj.Ts;
               
        end % function
        
        function plot(obj)
            
            % create figure
            h = figure( 'units', 'pixels', ...
                'Name', [obj.StudyName ', channel ' num2str(obj.StudyChannelNumber)], ...
                'Position', get(0,'ScreenSize'), ... % make it fullscreen
                'Tag', 'MeasurementPlot', ...
                'Visible', 'off'); % plotting while not visible is faster
            
            ax = axes();
            
            tempData = obj.Data;
            time_measurement = 0:obj.Ts:(obj.Length-1)*obj.Ts;
            tempTime = time_measurement;
            
            % replace seizure data with nan so it can be plotted in another
            % color
            for k = 1:length(obj.SeizureStart)
                tempTime(obj.SeizureStart(k):obj.SeizureEnd(k)) = nan;
                tempData(obj.SeizureStart(k):obj.SeizureEnd(k)) = nan;
            end % for k
            l = plot(tempTime, tempData);
            legendstr = {'EEG Data'};
            
            if ~isempty(obj.SeizureStart)
                % plot seizures in a different color
                tempData = nan(size(obj.Data));
                tempTime = nan(size(time_measurement));
                
                for k = 1:length(obj.SeizureStart)
                    tempTime(obj.SeizureStart(k):obj.SeizureEnd(k)) = ...
                        time_measurement(obj.SeizureStart(k):obj.SeizureEnd(k));
                    tempData(obj.SeizureStart(k):obj.SeizureEnd(k)) = ...
                        obj.Data(obj.SeizureStart(k):obj.SeizureEnd(k));
                end % for k
                hold(ax, 'on')
                l(end+1) = plot(tempTime, tempData, 'r');
                hold(ax, 'off')
                legendstr{end+1} = 'Seizure';
            else
                warning('Measurement contains no seizure.')
            end % if ~isempty
            
            legend(ax, l, legendstr)
            xlabel(ax, 'Time [s]')
            ylabel(ax, 'EEG Data')
            grid(ax, 'on')
            h.Visible = 'on';
            
        end % function plot
        
        function details(obj)
            % print out some information about the measurement
            
            days = floor(obj.Duration/(3600*24));
            hours = floor((obj.Duration-days*24*3600)/3600);
            minutes = floor((obj.Duration-days*24*3600-hours*3600)/60);
            seconds = obj.Duration-days*24*3600-hours*3600-minutes*60;
            
            fprintf('** Measurement %s channel %d\n', obj.StudyName, obj.StudyChannelNumber)
            fprintf('Number of seizures: %d\n', obj.NumSeizures)
            fprintf('Duration: %d days %d hours %d minutes %.2f seconds\n', days, hours, minutes, seconds)
            
        end % function details
        
        function out = isSeizure(obj, queryIndices)
            % returns an logical (boolean) array that is true at index i
            % if that index corresponds to a seizure in the measurement and
            % false otherwise
            % input argument queryIndices can be ommited, the function will
            % then query over the whole measurement
            
            
            if nargin == 1
                queryIndices = 1:obj.Length;
            end
            
            out = false(length(queryIndices), 1);
            
            for k = 1:length(obj.SeizureStart)
                    out(obj.SeizureStart(k):obj.SeizureEnd(k)) = true;
            end % for k

            out = out(queryIndices); 
            % yes, it is inefficient to calculate everything and then cut
            % it out but this function is just a quick fix
            
        end % function isSeizure
    end
    
    methods (Access = protected)
        
       function checkData(obj)
            % when setting the Data property, check that the input does not
            % contain any NaN entries. If it does, replace them with their
            % nearest non-NaN neighbor.     
            % NaN values cause issues in certain functions
            if any(isnan(obj.Data))
                disp('Data contains NaN entries. Replacing them...')
                x = transpose(1:obj.Length);
                obj.isNaN = isnan(obj.Data);
                obj.Data = interp1(x(~obj.isNaN), obj.Data(~obj.isNaN), x, 'nearest', 'extrap');
            end % if            
        end % function set.Data
        
    end % methods (Access = protected)
    
end

