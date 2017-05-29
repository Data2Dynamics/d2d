% create line styles and color
%
% C = arLineMarkersAndColors(j,n,color,markerstyle,linestyle)

function C = arLineMarkersAndColors(j,n,color,markerstyle,linestyle)

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
    elseif(n==5)
        colors = [0,0,.17; 1,0,0; 0,0,1; 1,0.1,0.72; 1.*.8,0.82*.8,0];
    elseif(n==6)
        colors = [0,0,.17; 1,0,0; 0,0,1; 1,0.1,0.72; 1.*.8,0.82*.8,0; 0,.344,0];
    elseif(n==7)
        colors = [0,0,.17; 1,0,0; 0,0,1; 1,0.1,0.72; 1.*.8,0.82*.8,0; 0,.344,0; 0, .7, 0];
    else
        colors = jet(n-1);
        colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2,2)));
        colors = [0 0 0; colors];
    end
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

C = {'Color',color,'LineStyle',linestyle,'Marker',markerstyle};


function j = mymod(k,N)
j = mod(k-1,N)+1;