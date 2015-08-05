function cmap = redwhiteblue(m)
% Interpolates Matlab's redbluecmap for a smooth color gradient from red
% through white to blue.

if nargin < 1, m = size(get(gcf,'colormap'),1); end

tmp = redbluecmap(11);
cmap(:,1) = interp1(linspace(0,1,size(tmp,1)), tmp(:,1), linspace(0,1,m));
cmap(:,2) = interp1(linspace(0,1,size(tmp,1)), tmp(:,2), linspace(0,1,m));
cmap(:,3) = interp1(linspace(0,1,size(tmp,1)), tmp(:,3), linspace(0,1,m));
