function increaseSize(hfig, varargin)
% INCREASESIZE increases the font size and the linewidth of the axes and
% lines in the figure whos handle is passed in the input argument
%
%   increaseSize()
%   increaseSize(hfig)
%
% Paul Meyer-Rachner
% paul@meyer-rachner.email
% 2018-02-18

%% parse input
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

defaultLineWidth = 1.5;
defaultFontSize = 20;

addParameter(p, 'LineWidth',defaultLineWidth, @isscalar);
addParameter(p, 'FontSize', defaultFontSize, @isscalar);
parse(p, varargin{:});
lineWidth = p.Results.LineWidth;
fontSize = p.Results.FontSize;

%% Make figure background white
hfig.Color = 'white';

%% increase size of plot elements
children = hfig.Children;
ax = findall(children, 'Type', 'Axes');

% do it for every ax found
for i = 1:length(ax)
    l = findall(ax(i).Children, 'Type', 'Line');
    s = findall(ax(i).Children, 'Type', 'Stair');

    % set Width of Plot and Stair Lines
    set(l, 'LineWidth', lineWidth)
    set(s, 'LineWidth', lineWidth)

    % increase text size of axes
    set(ax(i), 'FontSize', fontSize)
    set(ax(i), 'LineWidth', 0.75)
end



