function [htext] = label(h,textString,varargin)
%LABEL places a label next to your data.  
% 
% This function provides an option between the legend and text or annotation commands
% for labeling data that you plot.  Edward Tufte
% says that data shouldn't stray far from its label, because
% the viewer of a graph should not need to repeatedly move his or her eyes
% back and forth between plotted data and the legend to connect the dots of
% which data are which.  In this spirit, label can be used to place a
% label directly on a plot close to the data it describes.  
%
%% Syntax 
% 
%  label(h,'string')
%  label(...,'location',LocationString)
%  label(...,'TextProperty',PropertyValue)
%  label(...,'slope')
%  h = label(...)
%
%% Description 
% 
% label(h,'string') places 'string' near the leftmost data described by
% handle h. 
%
% label(...,'location',LocationString) specifies location of the string.
% LocationString can be any of the following:
% 
% * 'left' or 'west' (default) 
% * 'right' or 'east' 
% * 'top' or 'north' 
% * 'bottom' or 'south' 
% * 'center' or 'middle' 
% 
% label(...,'TextProperty',PropertyValue) specifies text properties as
% name-value pairs. 
%
% label(...,'slope') attempts to angle text following the local slope of
% the data. 
%
% htext = label(...) returns the handle htext of the newly-created text
% object. 
% 
%% Author Info
% Written by Chad A. Greene of the University of Texas at Austin and its 
% Institute for Geophysics, July 2014. 
% Fixed for R2014b in January 2015. 
% 
% See also annotation, text, and legend. 
%% Initial input error checks: 

assert(ishandle(h)==1,'Unrecognized object handle.')
assert(isempty(get(0,'children'))==0,'No current axes are open.') 
assert(isnumeric(textString)==0,'Label given by textString must be a string.') 
assert(nargin>=2,'Must input an object handle and corresponding label string.') 

%% Set defaults: 

location = 'left'; 
followSlope = false; 

%% Parse inputs

tmp = strncmpi(varargin,'loc',3); 
if any(tmp)
    location = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
end

tmp = strcmpi(varargin,'slope'); 
if any(tmp) 
    followSlope = true; 
    varargin = varargin(~tmp); 
end


%% 

color = get(h,'color'); 
xdata = get(h,'XData'); 
assert(isvector(xdata)==1,'Plotted data must be vector or scalar.') 
ydata = get(h,'YData'); 

gcax = get(gca,'xlim'); 
gcay = get(gca,'ylim'); 

if followSlope
    pbp = kearneyplotboxpos(gca); % A modified version of Kelly Kearney's plotboxpos function is included as a subfunction below.  

    % slope is scaled because of axes and plot box may not be equal and square:
    gcaf = pbp(4)/pbp(3); 
    apparentTheta = atand(gcaf*gradient(ydata,xdata).*(gcax(2)-gcax(1))/(gcay(2)-gcay(1)));

end

% Find indices of data within figure window: 
ind = find(xdata>=gcax(1)&xdata<=gcax(2)&ydata>=gcay(1)&ydata<=gcay(2)); 

switch lower(location)
    case {'left','west','leftmost','westmost'}
        horizontalAlignment = 'left'; 
        verticalAlignment = 'bottom'; 
        x = min(xdata(ind));
        y = ydata(xdata==x);
        textString = [' ',textString]; 
        xi = xdata==x; 
        
    case {'right','east','rightmost','eastmost'}
        horizontalAlignment = 'right'; 
        verticalAlignment = 'bottom'; 
        x = max(xdata(ind)); 
        y = ydata(xdata==x);
        textString = [textString,' ']; 
        xi = xdata==x(1); 
        
    case {'top','north','topmost','northmost'}
        horizontalAlignment = 'left'; 
        verticalAlignment = 'top'; 
        y = max(ydata(ind));
        x = xdata(ydata==y);
        xi = xdata==x(1); 
        
    case {'bottom','south','southmost','bottommost'} 
        horizontalAlignment = 'left'; 
        verticalAlignment = 'bottom'; 
        y = min(ydata(ind));
        x = xdata(ydata==y);
        xi = xdata==x(1); 
        
    case {'center','central','middle','centered','middlemost','centermost'}
        horizontalAlignment = 'center'; 
        verticalAlignment = 'bottom'; 
        xi = round(mean(ind)); 
        if ~ismember(xi,ind)
            xi = find(ind<xi,1,'last'); 
        end
        x = xdata(xi); 
        y = ydata(xi);
        
        
    otherwise
        error('Unrecognized location string.') 
end
 
% Set rotation preferences: 
if followSlope
    theta = apparentTheta(xi); 
else
    theta = 0; 
end


% Create the label: 
htext = text(x(1),y(1),textString,'color',color,'horizontalalignment',horizontalAlignment,...
    'verticalalignment',verticalAlignment,'rotation',theta); 

% Add any user-defined preferences: 
if length(varargin)>1 
    set(htext,varargin{:});
end


% Clean up: 
if nargout == 0
    clear htext
end

end


function pos = kearneyplotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region. THIS IS A
%SLIGHTLY MODIFIED VERSION OF KELLY KEARNEY'S PLOTBOXPOS FUNCTION. 
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    dx = diff(get(h, 'XLim'));
    dy = diff(get(h, 'YLim'));
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis

temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', get(h, 'parent'));
% set(temp, 'Units', currunit); % <-This line commented-out by Chad Greene, specifically for label function.  
pos = get(temp, 'position');
delete(temp);
end
