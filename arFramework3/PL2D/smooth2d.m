function smooth2d(xnodes)

% smooth2d(xnodes)
%
% Generates a 2D-Profile on a regular grid.
%
% See also: gen2d, vplfrom2d

global ar

if ~exist('xnodes','var') 
    xnodes = ar.ple2d.config.plot.xnodes;
end

%Unprocessed data points from ar.ple2d:
try
    x = ar.ple2d.raw.plpar;
    y = ar.ple2d.raw.predsteps;
    z = ar.ple2d.raw.chi2;
catch
    sprintf(['\n ERROR smooth2d: No 2D-Profile  found. \n, '...
        'Calculate it first by calling Init2DPL and gen2d \n'])
    return
end

% Impose grid data structure by interpolating each parameter profile
% separately
xs = linspace(min(x,[],'all'),max(x,[],'all'),xnodes);
ys = y;
[xq,yq] = meshgrid(xs,ys);
zq = NaN(length(ys),xnodes);
for ii = 1:length(ys)
    q_tmp = ~isnan(x(:,ii));
    x_tmp = x(q_tmp,ii);
    z_tmp = z(q_tmp,ii);
    if length(x_tmp) ~= length(unique(x_tmp))
        [x_tmp,ind_unis] = unique(x_tmp);
        z_tmp = z_tmp(ind_unis);
    end
    zq(ii,:) = interp1(x_tmp,z_tmp,xs');
end

offset = min(zq,[],'all');
zq = zq-offset; %shift global minimum to chi2 = 0;

ar.ple2d.smooth.xq = xq;
ar.ple2d.smooth.yq = yq;
ar.ple2d.smooth.zq = zq;
ar.ple2d.smooth.offset = offset;

end


