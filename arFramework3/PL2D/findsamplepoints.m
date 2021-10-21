function [x_sample,y_sample] = findsamplepoints(x,y,levels,sigma,weakness)
% [x_sample,y_sample] = findsamplepoints(x,y,levels,sigma,weakness)
%
% Find intersections of a curve (x,y) with the points specified in levels.
% Additional arguments sigma and weakness add a weak quadratic prior to the function.

if(exist('weakness', 'var')&&exist('sigma', 'var'))
    [y,x] = addweakprior(x,y,sigma,weakness);
end

x_sample = [];
y_sample = [];

for ii = 1:length(levels)
    inters = InterX([x;y],[x;repmat(levels(ii),1,length(x))]);
    x_sample = [x_sample,inters(1,:)];
    y_sample = [y_sample,inters(2,:)];
end

[x_sample,perm_ind] = sort(x_sample);
y_sample = y_sample(perm_ind);
%Sort the sample again in increasing x order
end

function [y,x] = addweakprior(x,y,sigma,weakness)
% Adds a weak quadratic prior to the curve.

[~,ind_min] = min(y); 
y = y + ((x-x(ind_min))/(weakness*sigma)).^2; %weak prior

end