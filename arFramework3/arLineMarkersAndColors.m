% create line styles and color
%
% C = arLineMarkersAndColors(j,n,color,markerstyle,linestyle)

function C = arLineMarkersAndColors(j,n,color,markerstyle,linestyle)

% test function
if(nargin==0 || isempty(j))
    if(~exist('color','var') || isempty(color))
        color = [];
    end
    if(~exist('markerstyle','var') || isempty(markerstyle))
        markerstyle = [];
    end
    if(~exist('linestyle','var') || isempty(linestyle))
        linestyle = [];
    end
    figure(1); clf;
    N = 20;
    h = [];
    for j=1:N
        C = arLineMarkersAndColors(j,color,markerstyle,linestyle);
        h(end+1) = plot(randn(1,10), C{:});
        hold on
    end
    hold off
    legend(h)
    return
end

if(n>7)
    n=7;
end

if(~exist('color','var') || isempty(color))
    colors = lines(n);
    colors = jet(n);
    colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2,2)));
    jc = mymod(j,length(colors));
    color = colors(jc,:);
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

C = {'LineStyle',linestyle,'Marker',markerstyle,'Color',color};


function j = mymod(k,N)
j = mod(k-1,N)+1;