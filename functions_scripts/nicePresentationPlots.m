% generate good looking plots for presentations

lim = [7.728 7.7375]*10^4;
position = [1 1 1440*0.9 900*0.6];

load('/Users/paul/Google Drive/Microchip_Biosignal_Computation/Seizure_Data/Study_005/Study_005_channel1.mat')
measurement.downsample(2)
% measurement_raw = IEEG_getData('Study 005', 1);
measurement_raw = measurement;


%measurement_raw.cut(1, measurement_raw.SeizureEnd(50)+2500*3);

measurement = measurement_raw.copy();

measurement.FilterOrder = 6;
measurement.preprocess();

%% Raw and Processed
measurement_raw.plot();
h1 = gcf();
ax1 = gca();
lines1 = findall(ax1, 'Type', 'line');

measurement.plot();
h2 = gcf();
ax2 = gca();
lines2 = findall(ax2, 'Type', 'line');

hfig = figure('Position', position);
ax_sub1 = subplot(2, 1, 1);
ax_sub2 = subplot(2, 1, 2);
copyobj(lines1, ax_sub1)
copyobj(lines2, ax_sub2)

delete(h1);
delete(h2);

legend(ax_sub1, 'iEEG Data', 'Seizure')
xlabel(ax_sub2, 'Time [s]')
ylabel(ax_sub1, 'Raw iEEG Data')
ylabel(ax_sub2, 'Filtered iEEG Data')
title(ax_sub1, 'Raw vs. Filtered iEEG Signal')
increaseSize(gcf, 'LineWidth', 0.75);
formatplot(hfig, 'LeftMargin', 0.1, 'BottomMargin', 0.12, 'TopMargin', 0.12, 'Gap', 0.04);
ax_sub1.XLim = lim;
ax_sub1.YLim = ax_sub2.YLim;
difference = ax_sub2.XTick(2) - ax_sub2.XTick(1);
ticks = 0:difference:difference*length(ax_sub2.XTickLabel);
ax_sub2.XTickLabel = cellstr(num2str(ticks'))';
hfig.Color = 'white';
hfig.Position = position;
export_fig('processed_EEG', '-jpg')
delete(hfig)

%% Features
lineLength = LineLength(measurement);
thetaBand = SpectralBandPower(measurement, 'ThetaBand', [4 8]);
alphaBand = SpectralBandPower(measurement, 'AlphaBand', [14 32]);
betaBand = SpectralBandPower(measurement, 'BetaBand', [8 12]);
nonlinearEnergy = NonlinearEnergy(measurement);

lineLength.calculate();
thetaBand.calculate();
alphaBand.calculate();
betaBand.calculate();
nonlinearEnergy.calculate();

lineLength.calculateOperatingCharacteristicCurve();
thetaBand.calculateOperatingCharacteristicCurve();
alphaBand.calculateOperatingCharacteristicCurve();
betaBand.calculateOperatingCharacteristicCurve();
nonlinearEnergy.calculateOperatingCharacteristicCurve();

lineLength.findBestThreshold();
thetaBand.findBestThreshold();
alphaBand.findBestThreshold();
betaBand.findBestThreshold();
nonlinearEnergy.findBestThreshold();

lineLength.plot('Measurement', 'Threshold', 'Baseline');
lineLength.AxesHandle(2).XLim = lim;
ax = lineLength.AxesHandle;
h = lineLength.PlotHandle;
ax(2).XLim = lim;
h.Color = 'white';
title(ax(1), 'Line Length Feature')
h.Position = position;
difference = ax(2).XTick(2) - ax(1).XTick(1);
ticks = 0:difference:difference*length(ax(2).XTickLabel);
ax(2).XTickLabel = cellstr(num2str(ticks'))';
increaseSize(gcf, 'LineWidth', 1);
formatplot(h, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
export_fig('lineLength', '-jpg')
close(h);

thetaBand.plot('Measurement', 'Threshold', 'Baseline');
thetaBand.AxesHandle(2).XLim = lim;
ax = thetaBand.AxesHandle;
h = thetaBand.PlotHandle;
ax(2).XLim = lim;
title(ax(1), 'Energy in Theta Band [4-8Hz] Feature')
ax(2).XTickLabel = cellstr(num2str(ticks'))';
h.Color = 'white';
h.Position = position;
increaseSize(gcf, 'LineWidth', 1);
formatplot(h, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
export_fig('thetaBand', '-jpg')
close(h);

alphaBand.plot('Measurement', 'Threshold', 'Baseline');
alphaBand.AxesHandle(2).XLim = lim;
ax = alphaBand.AxesHandle;
h = alphaBand.PlotHandle;
ax(2).XLim = lim;
title(ax(1), 'Energy in Alpha Band [8-14Hz] Feature')
ax(2).XTickLabel = cellstr(num2str(ticks'))';
h.Color = 'white';
h.Position = position;
increaseSize(gcf, 'LineWidth', 1);
formatplot(h, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
export_fig('alphaBand', '-jpg')
close(h);

betaBand.plot('Measurement', 'Threshold', 'Baseline');
betaBand.AxesHandle(2).XLim = lim;
ax = betaBand.AxesHandle;
h = betaBand.PlotHandle;
ax(2).XLim = lim;
title(ax(1), 'Energy in Beta Band [14-32Hz] Feature')
ax(2).XTickLabel = cellstr(num2str(ticks'))';
h.Color = 'white';
h.Position = position;
increaseSize(gcf, 'LineWidth', 1);
formatplot(h, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
export_fig('betaBand', '-jpg')
close(h);

nonlinearEnergy.plot('Measurement', 'Threshold', 'Baseline');
nonlinearEnergy.AxesHandle(2).XLim = lim;
ax = nonlinearEnergy.AxesHandle;
h = nonlinearEnergy.PlotHandle;
ax(2).XLim = lim;
h.Color = 'white';
h.Position = position;
title(ax(1), 'Nonlinear Energy Feature')
ax(2).XTickLabel = cellstr(num2str(ticks'))';
increaseSize(h , 'LineWidth', 1);
formatplot(h, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
export_fig('nonlinearEnergy', '-jpg')
close(h);

% lineLength.ThresholdBaselineFactor = 4;
% lineLength.plot('Measurement', 'Threshold', 'Baseline');
% ax = lineLength.AxesHandle;
% h = lineLength.PlotHandle;
% ax(2).XLim = lim;
% ax(1).YLim = [-500 500];
% ax(1).YTick = [];
% ax(2).YTick = [];
% ax(1).XTick = [];
% ax(2).XTick = [];
% ylabel(ax(1), 'iEEG Data')
% ylabel(ax(2), 'Feature')
% title(ax(1), 'Feature Analysis')
% ax(2).Legend.String{1} = 'Feature';
% h.Color = 'white';
% h.Position = position;
% increaseSize(gcf, 'LineWidth', 2, 'FontSize', 24);
% formatplot(h, 'LeftMargin', 0.08, 'BottomMargin', 0.1, 'TopMargin', 0.1);
% export_fig('exampleFeatureAnalysis', '-jpg')
% close(h);

%% Feature to logical

hfig = figure('Position', position);
ax_sub1 = subplot(2, 1, 1);
ax_sub2 = subplot(2, 1, 2);

time_feature = 0:alphaBand.StepSize*1/alphaBand.MeasurementFs:(alphaBand.Length-1)*alphaBand.StepSize*1/alphaBand.MeasurementFs;

alphaBand.plotFeature(ax_sub1, 'Threshold');
stairs(ax_sub2, time_feature, alphaBand.BitAboveThreshold);

hfig.Color = 'white';
title(ax_sub1, 'Conversion to logical values')
hfig.Position = position;
increaseSize(hfig, 'LineWidth', 1);
formatplot(hfig, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
ax_sub2.XLim = lim;
ax_sub2.YLim = [-0.1 1.1];
difference = ax_sub2.XTick(2) - ax_sub2.XTick(1);
ticks = 0:difference:difference*length(ax_sub2.XTickLabel);
ax_sub2.XTickLabel = cellstr(num2str(ticks'))';
xlabel(ax_sub2, 'Time [s]')
ylabel(ax_sub2, 'Logical Value')
export_fig('thresholding', '-jpg')
close(hfig);

%% Matrix

load('/Users/paul/Google Drive/Microchip_Biosignal_Computation/Seizure_Data/Matrices/Study_005_matrix.mat')
matrix.fillMatrices

channels = [1 4 3];
cols = [1 2 4];

fig = figure('Position', position);
ax = gobjects(length(cols)*length(channels), 1);
for k = 1:length(ax)
    ax(k) = subplot(3, 3, k);
end % for
lim_samples = lim*matrix.Matrix{1}.MeasurementFs/matrix.Matrix{1}.StepSize;
index_samples = lim_samples(1):lim_samples(2);

for k = 1:length(ax)
    [r, c] = ind2sub([length(channels), length(cols)], k);
    bits = squeeze(matrix.FlagMatrix(channels(r), cols(c), index_samples));   
    plot(ax(k), bits)
end % for

increaseSize(fig, 'linewidth', 1.5);

for k = 1:length(ax)
    ax(k).XTickLabel = [];
    [r, c] = ind2sub([length(channels), length(cols)], k);
    if c ~= 1
        ax(k).YTickLabel = [];
    else
        ylabel(ax(k), ['Channel ' num2str(r)])
    end %
    if r == 1
        title(ax(k), ['Feature ' num2str(c)])
    end 
    ax(k).Position(1) = 0.1 + 0.3*(c-1);
    ax(k).Position(2) = 0.95 - 0.3*(r);
    ax(k).Position(3) = 0.29;
    ax(k).Position(4) = 0.29;
end % for k

linkaxes(ax, 'xy')
ax(1).YLim = [-0.1 1.1];
ax(1).XLim = [1 length(index_samples)];

fig.Color = 'w';
export_fig('matrix_flags', '-jpg')
close(fig)

%% Plot weighted detection

mdOLS = MatrixDetection_LinearRegression(matrix);
mdOLS.CalculationFunctionHandle = @(X, y) (X'*X + 0.01*eye(size(X'*X)))\X'*y;
mdOLS.calculateWeights();
mdOLS.calculate();

mdOLS.plotWeightedDetection('Threshold');
fig = gcf();
fig.Position = position;
ax = gca();
title(ax, 'Weighted Seizure Detection')
fig.Color = 'w';
increaseSize(fig);
formatplot(fig, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
ax.XLim = lim;
ax.YLim = [0 100];
difference = ax.XTick(2) - ax.XTick(1);
ticks = 0:difference:difference*length(ax.XTickLabel);
ax.XTickLabel = cellstr(num2str(ticks'))';
export_fig('weighted_detection', '-jpg')
delete(fig)

%% Plot Curve

[h, ax] = mdOLS.CharacteristicCurve.plot('Weighted Matrix Detection');
matrix.bestFeature().CharacteristicCurve.addCurve(h);
ax.Children(1).LineStyle = '--';

title('Sensitivity vs. false alarm rate');
increaseSize(h, 'LineWidth', 1.5);
formatplot(h, 'LeftMargin', 0.08, 'BottomMargin', 0.12, 'TopMargin', 0.12);
ax.XLim = [0 1];
export_fig('sensitivity_plot', '-jpg')
delete(h)