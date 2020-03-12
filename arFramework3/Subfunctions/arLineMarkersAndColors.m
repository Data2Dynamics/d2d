% C = arLineMarkersAndColors([j,n,color,markerstyle,linestyle])
% create line styles and color

function C = arLineMarkersAndColors(j,n,color,markerstyle,linestyle)
global ar

% test function
if(nargin==0 || isempty(j))
    color = [];
    markerstyle = '.';
    linestyle = '-';

    figure(1); clf;
    N = 4;
    h = [];
    for j=1:N
        C = arLineMarkersAndColors(j,N,color,markerstyle,linestyle);
        h(end+1) = plot(randn(1,10), C{:}); %#ok<AGROW>
        hold on
    end
    hold off
    legend(h)
    return
end

Nmax = 9;
if(n>Nmax)
    n=Nmax;
end

if ~isfield(ar.config, 'plotColorSet')
    ar.config.plotColorSet = 'matlab_default';
end

% Legacy colors
if strcmp(ar.config.plotColorSet, 'd2d_legacy')
if(~exist('color','var') || isempty(color))
    if(n==1)
        colors = [0 0 0];
    elseif(n==2)
        colors = [0 0 0; 1 0 0];
    elseif(n==3)
        colors = [0 0 0; 1 0 0; 0 0 1];
    elseif(n==4)
        colors = [0 0 0; 1 0 0; 1 0 1; 0 0 1];
        colors(3,:) = bsxfun(@rdivide, colors(3,:), ...
            sqrt(sum(colors(3,:).^2,2)));
%     elseif(n==5)
%         colors = [0,0,.17; 1,0,0; 0,0,1; 1,0.1,0.72; 1.*.8,0.82*.8,0];
%     elseif(n==6)
%         colors = [0,0,.17; 1,0,0; 0,0,1; 1,0.1,0.72; 1.*.8,0.82*.8,0; 0,.344,0];
%     elseif(n==7)
%         colors = [0,0,.17; 1,0,0; 0,0,1; 1,0.1,0.72; 1.*.8,0.82*.8,0; 0,.344,0; 0, .7, 0];
    else
        colors = jet(n-1);
        colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2,2)));
        colors = [0 0 0; colors];
    end
    jc = mymod(j,length(colors));
    color = colors(jc,:);
end

% dMod colors
elseif strcmp(ar.config.plotColorSet, 'dMod')
if(~exist('color','var') || isempty(color))
    if n < 11
        colorvec = [0 0 0; 
            197 0 11; 
            0 132 209; 
            87,157,28; 
            255,149,14; 
            75,31,111; 
            204,121,167; 
            0,100,0; 
            240,228,66; 
            139,69,19]/255;
        colors = colorvec(1:n,:);
    else
        colors = jet(n-1);
        colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2,2)));
        colors = [0 0 0; colors];
    end
    jc = mymod(j,length(colors));
    color = colors(jc,:);
end

% MATLAB default colors
elseif strcmp(ar.config.plotColorSet, 'matlab_default')
if(~exist('color','var') || isempty(color))
    if(n==1)
        colors = [0 0 0];
    elseif n < 7
        colorvec = [0, 0.4470, 0.7410;
            0.8500, 0.3250, 0.0980;
            0.9290, 0.6940, 0.1250;
            0.4940, 0.1840, 0.5560;
            0.4660, 0.6740, 0.1880;
            0.3010, 0.7450, 0.9330;
            0.6350, 0.0780, 0.1840];
            colors = colorvec(1:n,:);
    else
    colors = jet(n-1);
    colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2,2)));
    colors = [0 0 0; colors];
    end
    jc = mymod(j,length(colors));
    color = colors(jc,:);
end
else
    error('Invalid value in ar.config.plotColorSet. Use ''matlab_default'', ''dMod'' or ''d2d_legacy''.')
end

if(~exist('markerstyle','var') || isempty(markerstyle))
    markers = {'.','o','x','+','*','s','d','v','^','<','>'};
    jm = mymod(j,length(markers));
    markerstyle = markers{jm};
end

if(~exist('linestyle','var') || isempty(linestyle))
    linestyles = {'-',':','-.','--'};   
    jl = mymod(j,length(linestyles));
    linestyle = linestyles{jl};
end

C = {'Color',color,'LineStyle',linestyle,'Marker',markerstyle};


function j = mymod(k,N)
j = mod(k-1,N)+1;