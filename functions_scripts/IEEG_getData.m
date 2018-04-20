function measurement = IEEG_getData(studyName, channel, loginName, loginBin)
% get data from IEEG.org and convert it to a measurement struct so it can
% be used by my featureEvaluation script
% please refer to the iEEG matlab toolbox documentation for more details on
% how to download data from the iEEG.org website

session = IEEGSession(studyName, loginName, loginBin);
dataset = session.data;

Fs = dataset.sampleRate;

annotation = session.data.annLayer;
annarray = annotation.getEvents(0);

% event times are in micro-seconds
start_ms = [];
stop_ms = [];

for i = 1:length(annarray)
    if strcmp(annarray(i).type, 'Seizure')
        start_ms(end+1) = annarray(i).start; %#ok<*AGROW>
        stop_ms(end+1) = annarray(i).stop;
    end % if
end % for i

buffer_hours = 0.5;
buffer_samples = floor(buffer_hours*3600*Fs);

start_samples = floor(start_ms*(10^-6)*Fs);
stop_samples = ceil(stop_ms*(10^-6)*Fs);

downloadStart = start_samples(1) - buffer_samples + 1;
downloadEnd = stop_samples(end) + buffer_samples - 1;

shift = -start_samples(1) + buffer_samples + 1;
start_samples = start_samples + shift;
stop_samples = stop_samples + shift;

% download in chunks of 100 millions samples
stepSize = 100e6;
data = [];

for k = downloadStart:stepSize:downloadEnd-stepSize+1
    data = [data; dataset.getvalues(k:k+stepSize-1, channel)];
end % for k
k = k+stepSize;

data = [data; dataset.getvalues(k:downloadEnd, channel)];

targetFs = round(Fs);

measurement = Measurement(data, targetFs, start_samples, stop_samples, 'StudyName', studyName, 'StudyChannelNumber', channel);
try
    measurement.StudyChannelName = dataset.channelLabels{channel, 1};
catch
    warning('Unable to add study channel name.')
end % try