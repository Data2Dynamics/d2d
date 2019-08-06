% Pay attention: Defintion of bounds depends on the data.
% 
%   This function uses time-offsets, i.e. is for the "retarded" transient
%   fucntion, which can be exemplified by the following:
% amp_sust = 1;
% amp_trans = 4;
% timescale_sust = 1;
% timescale_trans = 1;
% offset_TF = 0;
% toffset = 0;
% tt = log10(10.^t+10.^toffset)-log10(1+10.^toffset);
% y1 = amp_trans.*(1-exp(-tt./timescale_sust)).*exp(-tt./(timescale_sust+timescale_trans));
% y2 = amp_sust.*(1-exp(-tt./timescale_sust));
% y = y1+y2+offset_TF;

function [bounds,boundsNegative] = DefaultLbUbTransient2(qPositive)
global ar
if ~exist('qPositive','var') || isempty(qPositive)
    qPositive = false;  % Is the truth known to be positive?
end

bounds.pLabel = ar.pLabel;
bounds.lb = ar.lb;
bounds.ub = ar.ub;

%% The transient fuctions looks like this:
%     y1 = p4*(1-exp(-t/timescale_sust)).*exp(-p2*t);
%     y2 = p3*(1-exp(-t/timescale_sust));
%     y = y1+y2+p(5);
% 
%     bound.lb = [.5/(t(end)-t(1)),.5/(t(end)-t(1)),1e-10,1e-10,  nanmin(data)-0.5*D];
%     bound.ub = [2/(t(2)-t(1)),2/(t(2)-t(1)),D*2,D*2,  nanmax(data)];

%% signum:
ind = ar.fit_transient.indp.signum;
for i=1:length(ind)
    bounds.lb(ind(i)) = -1;
    bounds.ub(ind(i)) = 1;  
end

%% time-offset: [because of rescaling of the time axis always between 0 and 10
ind = ar.fit_transient.indp.toffset;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    bounds.lb(ind(i)) = nanmin(t)-1;
    bounds.ub(ind(i)) = nanmin(t)+5; % half of the time course, (after rescaling the new time axis has length=10  
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));        
    end
end

%% amplitudes
ind = ar.fit_transient.indp.amp_sust;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    if D==0
        D = 1e-6;
    end
    bounds.lb(ind(i)) = D*1e-6;
    bounds.ub(ind(i)) = D*2;
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

ind = ar.fit_transient.indp.amp_trans;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    if D==0
        D = 1e-6;
    end
    bounds.lb(ind(i)) = D*1e-6;
%     bounds.ub(ind(i)) = D*4; % this choice might be better for more dense data
    bounds.ub(ind(i)) = D*2; % this choice is more conservative and might be better for roughly sampled data.
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

%% offsets
indOffset = ar.fit_transient.indp.offset;
for i=1:length(indOffset)
    [t,y,D] = getData(ar.pLabel{indOffset(i)});
    bounds.lb(indOffset(i)) = nanmin(y)-0.5*D;
%          bounds.lb(indOffset(i)) = max(0, bounds.lb(indOffset(i)));
    
    bounds.ub(indOffset(i)) = nanmax(y);
    if(ar.qLog10(indOffset(i))==1)
        bounds.lb(indOffset(i)) = log10(bounds.lb(indOffset(i)));
        bounds.ub(indOffset(i)) = log10(bounds.ub(indOffset(i)));
    elseif qPositive
        bounds.lb(indOffset(i)) = max(0,bounds.lb(indOffset(i)));
        bounds.ub(indOffset(i)) = max(0,bounds.ub(indOffset(i)));
    end
end
if(sum(ar.qLog10(indOffset(i)))>0)
    warning('Negative offsets are often required. Please set ar.qLog10 for offset parameters =0.');
end

%% timescales
ind = ar.fit_transient.indp.timescale_sust;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    if length(t)<10
        bounds.lb(ind(i)) = (t(2)-t(1))/2;
    else
        bounds.lb(ind(i)) = (t(2)-t(1));
    end
    bounds.lb(ind(i)) = max(bounds.lb(ind),range(t)/100);
    
    bounds.ub(ind(i)) = 2*range(t);
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

ind = ar.fit_transient.indp.timescale_trans;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    bounds.lb(ind(i)) = (t(2)-t(1))/2;
    bounds.lb(ind(i)) = max(bounds.lb(ind),range(t)/100);

    % divide by two because
    % tau_trans=tau_trans+taus_sus in the transient function
    bounds.lb(ind(i)) = bounds.lb(ind(i))/2;
    
    bounds.ub(ind(i)) = 2*range(t);
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

%% SDs
ind = ar.fit_transient.indp.sd;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});

    % lower bound:
    R = range(y);
    if R>0 && ~isnan(R) && ~isinf(R)
        bounds.lb(ind(i)) = R / 1e4;
        if(ar.qLog10(ind(i))==1)
            bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        end
    else % keep bounds from ar.lb
        warning('Constant data: SD bounds are not adjusted.');
    end
    
    % upper bound:
    SD = nanstd(y);
    % handle case SD=0 or SD too small:
    sdlb = bounds.lb(ind(i));
    if(ar.qLog10(ind(i))==1)
        sdlb = 10.^sdlb;
    end
    if isnan(SD) || SD<sdlb
        SD = sdlb*10;
    end
        
    bounds.ub(ind(i)) = SD;
    
    if(ar.qLog10(ind(i))==1)
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

%% for fits in negative direction:
boundsNegative = bounds;

if(~isempty(indOffset))
    ub_unlog = NaN(size(indOffset));
    for i=1:length(indOffset)
        if(ar.qLog10(indOffset(i))==1)
            boundsNegative.lb(indOffset(i)) = log10(nanmin(y));
            boundsNegative.ub(indOffset(i)) = log10(nanmax(y)+0.5*D);
            ub_unlog(i) = 10.^boundsNegative.ub(indOffset(i));
        else
            boundsNegative.lb(indOffset(i)) = nanmin(y);
            boundsNegative.ub(indOffset(i)) = nanmax(y)+0.5*D;
            if qPositive
                boundsNegative.lb(indOffset(i)) = max(0,boundsNegative.lb(indOffset(i)));
                boundsNegative.ub(indOffset(i)) = max(0,boundsNegative.ub(indOffset(i)));
            end
            ub_unlog(i) = boundsNegative.ub(indOffset(i));            
        end
    end
    max_ub_offset = max(ub_unlog);
else
    max_ub_offset = Inf;
end

ind = [ar.fit_transient.indp.amp_sust,ar.fit_transient.indp.amp_trans];
if(~isempty(ind))
    for i=1:length(ind)
        if(ar.qLog10(ind(i))==0) && qPositive
            boundsNegative.lb(ind(i)) = min(max_ub_offset,boundsNegative.lb(ind(i)));
            boundsNegative.ub(ind(i)) = min(max_ub_offset,boundsNegative.ub(ind(i)));
        elseif qPositive  % log
            boundsNegative.lb(ind(i)) = max(-10,log10(min(max_ub_offset,10^boundsNegative.lb(ind(i)))));
            boundsNegative.ub(ind(i)) = max(-10,log10(min(max_ub_offset,10^boundsNegative.ub(ind(i)))));
        end
    end
end

%% check:
if sum(abs(imag(bounds.lb))>0)>0
    error('bounds.lb are complex');
elseif sum(abs(imag(bounds.lb))>0)>0
    error('bounds.ub are complex');
elseif sum(abs(imag(boundsNegative.lb))>0)>0
    error('boundsNegative.lb are complex');
elseif sum(abs(imag(boundsNegative.lb))>0)>0
    error('boundsNegative.ub are complex');
end



function [t,y,D] = getData(pname)
global ar
y = ar.fit_transient.(pname).yexp;
t = ar.fit_transient.(pname).t;

t = unique(t); 
t = t*10/range(t); % the factor for rescaling time 
D = range(y);
