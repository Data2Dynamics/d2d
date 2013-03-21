function [nout, xbins] = hist_norm(x, nbins,silent)

if(nargin<2)
    nbins = max([10 length(x)/10]);
end
if(nargin<3)
    silent = false;
end

if(length(nbins) == 1)
    xmax = max(x);
    xmin = min(x);
    xbinsdiff = (xmax-xmin) * 0.05;
    xbins = linspace(xmin-xbinsdiff, xmax+xbinsdiff, nbins);
elseif(length(nbins) > 1)
    xbins = nbins;
end

nb = histc(x,xbins);
nnorm = sum(nb(1:end))*(xbins(2)-xbins(1));
nout = nb(1:end)/nnorm;

if(~silent)
    h = bar(xbins,nout,'histc');
    set(h, 'FaceColor','none');
end