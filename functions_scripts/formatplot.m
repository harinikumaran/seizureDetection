function [hfig, hax] = formatplot(hfig, varargin)
% FORMATPLOT reduces the gap between the elements of the current figure.
% The figure must contain a column of subplots that have the same x-axis
% data since FORMATPLOT removes the x-axis ticks and label of all subplots
% exept the bottom one.
%
% FORMATPLOT applies the function to the current figure
%
% FORMATPLOT(hfig) applies the function to the figure with handle hfig
%
% [hfig, hax] = FORMATPLOT(hfig) additionally returns the figure and axes
% handles
%
% [hfig, hax] = FORMATPLOT(hfig, 'Name', 'Value') accepts Name-Value pairs
% to customize the formatting
%
%
% Possible Name-Value pairs:
% -------------------------------------------------------------------------
% Name              Type      Description (default)
%
% TopMargin         double      normalized margin at the top of the figure
%                               (0.05)
%
% BottomMargin      double      normalized margin at the bottom of the 
%                               figure (0.08)
%
% LeftMargin        double      normalized margin at the left of the figure
%                               (0.08)
%
% RightMargin       double      normalized margin at the right of the 
%                               figure (0.02)
%
% Gap               double      normalized vertical gap between subplots 
%                               (0.02)
%
% LegendOutside     boolean     places the legend on the right next to the
%                               plots (false)
%
% LegendOutsideGap  double      normalized gap between legends and plots if
%                               legends are placed outside of plots (0.02)
%
% AllowTextcut      boolean     allows the axis labels and titles to be
%                               cutoff by the side of the figure (false)
% -------------------------------------------------------------------------
%
% see also SUBPLOT, PLOT
%
%
% Usage example:
%
%   x = linspace(0,4*pi,100);
%   h = figure();
%
%   subplot(3, 1, 1);
%   plot(x, sin(x));
%   legend('sin(x)');
%
%   subplot(3, 1, 2);
%   plot(x, cos(x));
%   legend('cos(x)');
%
%   subplot(3, 1, 3);
%   plot(x, sin(2*x));
%   legend('sin(2*x)');
%   xlabel('x [rad]');
%
%   formatplot(h, 'LegendOutside', true);
%
%
% -------------------------------------------------
% developed using Matlab Version 9.2.0 (R2017a)
% v1.2 - Apr 2018
% Paul Meyer-Rachner - paul@meyer-rachner.email
% -------------------------------------------------

% Version Log:
% -------------------------------------------------------------------------
% v1.2 - 2018-04-11:
% - added 'AllowTextcut' option. False by default. When the plots are
% repositioned, avoid cutting off title or axis labels, which can happen
% when the fontsize is pretty big and the gaps to the edged too small.
%
% v1.1.2 - 2017-10-25:
% - added 'limitrate' argument to drawnow. This speeds up the drawnow
% functions for plots with a high amount of data.
%
% v1.1.1 - 2017-05-05:
% - added example to help
%
% v1.1 - 2017-04-04:
% - the xlabels of all but the bottom subplot are now removed
% -------------------------------------------------------------------------

%% Check version

if verLessThan('matlab', '9.0')
    msg = ['You are using an older version of Matlab. ' ...
        'This function might not work as intended.'];
    warning(msg)
end % if verLessThan

%% Check input arguments
% if there are no input arguments, use current figure
if nargin == 0
    hfig = gcf;
    if isempty(hfig)
        % no figure open
        error('No figure open.')
    end % if isempty(hfig)
end % nargin == 0

% Parse Input
% Initialize inputParser
p = inputParser;

addParameter(p, 'TopMargin', 0.05, ...
    @(x) validateattributes(x, {'numeric', 'nonempty'}, ...
    {'nonempty', 'scalar', '>', 0, '<', 0.5}, mfilename))
addParameter(p, 'BottomMargin', 0.08, ...
    @(x) validateattributes(x, {'numeric', 'nonempty'}, ...
    {'nonempty', 'scalar', '>', 0, '<', 0.5}, mfilename))
addParameter(p, 'LeftMargin', 0.08, ...
    @(x) validateattributes(x, {'numeric', 'nonempty'}, ...
    {'nonempty', 'scalar', '>', 0, '<', 0.5}, mfilename))
addParameter(p, 'RightMargin', 0.02, ...
    @(x) validateattributes(x, {'numeric', 'nonempty'}, ...
    {'scalar', '>', 0, '<', 0.5}, mfilename))
addParameter(p, 'Gap', 0.02, ...
    @(x) validateattributes(x, {'numeric', 'nonempty'}, ...
    {'nonempty', 'scalar', '>', 0, '<', 0.2}, mfilename))
addParameter(p, 'LegendGap', 0.02, ...
    @(x) validateattributes(x, {'numeric', 'nonempty'}, ...
    {'nonempty', 'scalar', '>', 0, '<', 0.1}, mfilename))
addParameter(p, 'LegendOutside', false, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, ...
    {'nonempty', 'scalar', 'binary'}, mfilename))
addParameter(p, 'AllowTextcut', false, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, ...
    {'nonempty', 'scalar', 'binary'}, mfilename))

% Parse input
parse(p, varargin{:});

% Get results
marginTop = p.Results.TopMargin;
marginBottom = p.Results.BottomMargin;
marginLeft = p.Results.LeftMargin;
marginRight = p.Results.RightMargin;
gap = p.Results.Gap;
gapLegend = p.Results.LegendGap;
legendOutside = p.Results.LegendOutside;
allowTextcut = p.Results.AllowTextcut;

%% Check axes
% All subplots have to be in one column for this function to execute
% To check if this requirement is fullfiled, determine the positions of
% each axes

% Get axes handle
hax = findall(hfig, 'Type', 'Axes');

% Stop if no axes are available
if isempty(hax)
    msg = 'No axes found in current figure.';
    error(msg)
end % if isempty

position = zeros(length(hax), 4);
for k = 1:length(hax)
   position(k, :) = hax(k).Position; 
end % for k

% if all x coordinate of the bottom left corner are equal, the requirement
% is fullfiled
if ~all(position(:, 1 == position(1, 1)))
    msg = ['All subplots have to be in one column for ' ... 
        'this function to execute.'];
    error(msg)
end

%% Sort axes from top to bottom

% sort them according to the y coordinate of the bottom left corner
% hax(1) will be highest subplot, hax(end) will be lowest subplot 
[~, index] = sort(position(:, 2), 'descend');
hax = hax(index);

%% Position subplots

% number of subplots
countPlots = length(hax);

% calculate minimum margin required to not cutoff axis labels
requiredMarginBottom = hax(end).Position(2) - hax(end).OuterPosition(2);
requiredMarginLeft = hax(1).Position(1) - hax(1).OuterPosition(1);
requiredMarginTop = hax(1).OuterPosition(2) + hax(1).OuterPosition(4) ...
    - hax(1).Position(2) - hax(1).Position(4);
requiredMarginRight = hax(1).OuterPosition(1) + hax(1).OuterPosition(3) ...
    - hax(1).Position(1) - hax(1).Position(3);

% enlarge margins as to not cutoff text unless specified
if (requiredMarginBottom > marginBottom) && ~allowTextcut
    marginBottom = requiredMarginBottom;
end % if
if (requiredMarginRight > marginRight) && ~allowTextcut
    marginRight = requiredMarginRight;
end % if
if (requiredMarginTop > marginTop) && ~allowTextcut
    marginTop = requiredMarginTop;
end % if
if (requiredMarginLeft > marginLeft) && ~allowTextcut
    marginLeft = requiredMarginLeft;
end % if
    
% calculate width and height of subplots
width = (1 - marginLeft - marginRight); 
height = (1 - marginTop - marginBottom - gap*(countPlots-1))/countPlots; 

% calculate height of first subplot
posy = 1 - marginTop - height;

% position subplot in loop
for k = 1:length(hax)
    hax(k).Position = [marginLeft posy width height];
    posy = posy - height - gap;
end % for k


%% Format axes

% Remove XTick und XLabel from all but last subplot
for k = 1:length(hax) - 1
    hax(k).XTickLabel = [];
    hax(k).XLabel = [];
end % for k

%% Format legend

% only do formatting if requested
if legendOutside
    % get all legends
    hleg = findall(hfig, 'Type', 'Legend');

    % move them outside
    for k = 1:length(hleg)
        hleg(k).Location = 'northeastoutside';
    end % 

    % update
    drawnow

    % position of legend
    positionLegend = zeros(length(hleg), 4);
    for k = 1:length(hleg)
       positionLegend(k, :) = hleg(k).Position;
    end % for k
    maxWidth = max(positionLegend(:, 3));
    width = 1 - marginLeft - marginRight - maxWidth - gapLegend;

    % Resize plots
    for k = 1:length(hax)
        hax(k).Position(3) = width;
    end % for k

    % Reposition legends
    for k = 1:length(hleg)
        hleg(k).Position(1) = marginLeft + width + gapLegend;
    end % for k
end % if legendOutside

%% Link axes and enable grid

% link axes
linkaxes(hax, 'x')

% enable grid
for k = 1:length(hax)
    grid(hax(k), 'on')
end % for k

% update plot
drawnow limitrate

