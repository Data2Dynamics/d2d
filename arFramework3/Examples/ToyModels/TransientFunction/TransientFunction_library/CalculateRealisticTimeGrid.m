% tgrid = CalculateRealisticTimeGrid(t,tlim)
% 
% Calculation of a realistic time grid (corresponding to exp. data)
% At least 20 time points are defined and max 100 time points.
% 
% The following rules are used:
%   - more than 100 times available: leave out every second until <100
%   - between 20 and 100 (exp.) times available: leave it
%   - between 4 and 20: Place times between consecutive times until >=20
%       using a quadratic polynomial based on three consecuitve times
%   - 3 times: quadratic polynomial fit (see code for details)
%   - 2 times, i.e. interval: selection time points at log time scale

function tgrid = CalculateRealisticTimeGrid(t,tlim)
if size(t,1)>1
    t = t';
end
if length(t)>5 % an experimental time grid is availble
    t = union(t,tlim);
end
t = t(~isnan(t));
t = t(~isinf(t));

if length(t)>100
    while length(t)>100% leave out every second
        tend = t(end);
        t = t(1:2:end); 
        t = union(t,tend);
    end
    tgrid = t;
elseif length(t)>=20 % leave it
    tgrid = t;
elseif length(t)>=4 % place points in between
    tgrid = t;
    while length(tgrid)<20
        tnew = tgrid;

        p = polyfit(1:3,tgrid(1:3),2); % fit a quadratic polynomial to first three times
        tnew = union(polyval(p,1.5),tnew); % place one time point between the first two
        
        for i=1:(length(tgrid)-3)
            p = polyfit(1:4,tgrid(i:(i+3)),2); % fit a quadratic polynomial to ths time spacing
            tnew = union(polyval(p,2.5),tnew);
        end
        
        p = polyfit(1:3,tgrid((end-2):end),2); % fit a quadratic polynomial to first three times
        tnew = union(polyval(p,2.5),tnew); % place one time point between the last two

        tgrid = tnew;
    end
elseif length(t)==3
    p = polyfit(1:3,t,2); % fit a quadratic polynomial to ths time spacing
    tgrid = polyval(p,linspace(1,3,20)); 
elseif length(t)==2
    if t(1)==0
        tgrid = logspace(log10(t(1)+1),log10(t(2)+1),20)-1;
    elseif min(t)>0
        tgrid = logspace(log10(t(1)),log10(t(2)),20);
    else
        dt = min(t)+1;
        tgrid = logspace(log10(t(1)-dt),log10(t(2)-dt),20)+dt;        
    end
else
    if tlim(1)==0
        tgrid = logspace(log10(tlim(1)+1),log10(tlim(2)+1),20)-1;
    elseif min(tlim)>0
        tgrid = logspace(log10(tlim(1)),log10(tlim(2)),20);
    else
        dt = min(tlim)+1;
        tgrid = logspace(log10(tlim(1)-dt),log10(tlim(2)-dt),20)+dt;        
    end
end    

if size(tgrid,2)>1
    tgrid = tgrid'; % make a column
end
