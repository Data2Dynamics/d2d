function cmap = redwhiteblue(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = floor(m/2);
         
r = [linspace(239,247,n) linspace(247,103,n)]';
g = [linspace(138,247,n) linspace(247,169,n)]';
b = [linspace(98,247,n) linspace(247,207,n)]';

cmap = flipud([r g b]./247);