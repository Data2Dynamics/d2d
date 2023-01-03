% Sometimes wider bounds are reasonable, i.e. less conservative
% restriction, e.g. if many data points are available. Then boundfactor can
% be set to >1.

function Initialize_FitTransient2(boundfactor,ms,qPositive)
global ar
if ~exist('qPositive','var') || isempty(qPositive)
    qPositive = false;  % Is the truth known to be positive?
end
if ~exist('boundfactor','var') || isempty(boundfactor)
    boundfactor = 1;  % standard bounds, designed for experimetnal data    
end
if ~exist('ms','var') || isempty(ms)
    ms = 1:length(ar.model);    
end

ar.fit_transient = struct;

ar.fit_transient.indp.toffset    = strmatch('toffset_TF',ar.pLabel);
ar.fit_transient.indp.signum    = strmatch('signum_TF',ar.pLabel);
ar.fit_transient.indp.offset    = strmatch('offset_TF',ar.pLabel);
ar.fit_transient.indp.amp_sust  = strmatch('amp_sust',ar.pLabel);
ar.fit_transient.indp.amp_trans = strmatch('amp_trans',ar.pLabel);
ar.fit_transient.indp.timescale_sust = strmatch('timescale_sust',ar.pLabel);
ar.fit_transient.indp.timescale_trans = strmatch('timescale_trans',ar.pLabel);
ar.fit_transient.indp.sd = strmatch('sd_TF',ar.pLabel);
ar.fit_transient.indp.dummy = strmatch('init____dummy___',ar.pLabel);


ar.fit_transient.indp.all = [ar.fit_transient.indp.toffset,...
    ar.fit_transient.indp.amp_sust,...
    ar.fit_transient.indp.amp_trans,...
    ar.fit_transient.indp.offset,...
    ar.fit_transient.indp.timescale_sust,...
    ar.fit_transient.indp.timescale_trans,...
    ar.fit_transient.indp.sd ,...
    ar.fit_transient.indp.signum];
ar.fit_transient.indp.all = ar.fit_transient.indp.all(:)';

if length(ar.fit_transient.indp.signum) ~= length(ar.fit_transient.indp.toffset) || ...
   length(ar.fit_transient.indp.signum) ~= length(ar.fit_transient.indp.offset) || ...
   length(ar.fit_transient.indp.signum) ~= length(ar.fit_transient.indp.amp_sust) || ...
   length(ar.fit_transient.indp.signum) ~= length(ar.fit_transient.indp.amp_trans) || ...
   length(ar.fit_transient.indp.signum) ~= length(ar.fit_transient.indp.timescale_sust) || ...
   length(ar.fit_transient.indp.signum) ~= length(ar.fit_transient.indp.timescale_trans)
   warning('The number of different parameter-types of the transient function should be equal. \nMay be detection by name does not work. Check parameter names!');
end 

% determine the parameter -> data relationship
% t, yexp, ystd is required for specifying reasonable bounds
pnames = ar.pLabel;
psym = arSym(pnames);
ar.fit_transient.mLink = {};
ar.fit_transient.dLink = {};
for i=1:length(pnames)
    ar.fit_transient.(pnames{i}) = struct;
    ar.fit_transient.(pnames{i}).mLink = [];
    ar.fit_transient.(pnames{i}).dLink = [];
    ar.fit_transient.(pnames{i}).t = [];
    ar.fit_transient.(pnames{i}).yexp = [];
    ar.fit_transient.(pnames{i}).ystd = [];
end
for m=ms
    for d=1:length(ar.model(m).data)
        vars = symvar(arSym(ar.model(m).data(d).fp));
        if length(ar.model(m).data(d).y)>1
            error('Model%i, data%i has multiple observations which is not yet implemented.',m,d)
        end
        
        [~,found] = intersect(psym,vars);
        for i=1:length(found)
            ar.fit_transient.(pnames{found(i)}).mLink = [ar.fit_transient.(pnames{found(i)}).mLink,m];
            ar.fit_transient.(pnames{found(i)}).dLink = [ar.fit_transient.(pnames{found(i)}).dLink,d];
            ar.fit_transient.(pnames{found(i)}).t    = [ar.fit_transient.(pnames{found(i)}).t;ar.model(m).data(d).tExp];
            ar.fit_transient.(pnames{found(i)}).yexp = [ar.fit_transient.(pnames{found(i)}).yexp;ar.model(m).data(d).yExp];
            ar.fit_transient.(pnames{found(i)}).ystd = [ar.fit_transient.(pnames{found(i)}).ystd;ar.model(m).data(d).yExpStd];
            ar.fit_transient.(pnames{found(i)}).ind_arp = found;
        end
    end
end

% dummy
ar.qFit(ar.fit_transient.indp.dummy) = 0;


% signums are constants and non-log:
ar.qFit(ar.fit_transient.indp.signum) = 2;
ar.qLog10(ar.fit_transient.indp.signum) = 0;

ar.qLog10(ar.fit_transient.indp.toffset) = 0;
ar.qLog10(ar.fit_transient.indp.amp_sust) = 1;
ar.qLog10(ar.fit_transient.indp.amp_trans) = 1;
ar.qLog10(ar.fit_transient.indp.offset) = 0;  % often, negative offsets required
ar.qLog10(ar.fit_transient.indp.sd) = 1;


for m=ms
    for d=1:length(ar.model(m).data)
        if sum(ar.model(m).data(d).logfitting~=0)>0
            error('The transient function is design for non-log fitting. If log-trans required, do it before/do it by hand.')       
        elseif sum(ar.model(m).data(d).logplotting~=0)>0
            error('The transient function is design for non-log fitting. If log-trans required, do it before/do it by hand.')       
        end
    end
end


[ar.fit_transient.bounds,ar.fit_transient.boundsNeg] = DefaultLbUbTransient2(qPositive);

ind = setdiff(ar.fit_transient.indp.all,ar.fit_transient.indp.signum);
for i=1:length(ind)
    if ar.qLog10(ind(i))==1
        ar.fit_transient.bounds.lb(ind(i)) = ar.fit_transient.bounds.lb(ind(i))-log10(boundfactor);
        ar.fit_transient.bounds.ub(ind(i)) = ar.fit_transient.bounds.ub(ind(i))+log10(boundfactor);
        ar.fit_transient.boundsNeg.lb(ind(i)) = ar.fit_transient.boundsNeg.lb(ind(i))-log10(boundfactor);
        ar.fit_transient.boundsNeg.ub(ind(i)) = ar.fit_transient.boundsNeg.ub(ind(i))+log10(boundfactor);
    
    else % nonlog
        if ar.fit_transient.bounds.lb(ind(i))>0
            ar.fit_transient.bounds.lb(ind(i)) = ar.fit_transient.bounds.lb(ind(i))/boundfactor;
        else
            ar.fit_transient.bounds.lb(ind(i)) = ar.fit_transient.bounds.lb(ind(i))*boundfactor;
        end
        if ar.fit_transient.bounds.ub(ind(i))>0
            ar.fit_transient.bounds.ub(ind(i)) = ar.fit_transient.bounds.ub(ind(i))*boundfactor;
        else
            ar.fit_transient.bounds.ub(ind(i)) = ar.fit_transient.bounds.ub(ind(i))/boundfactor;
        end
        
        if ar.fit_transient.boundsNeg.lb(ind(i))>0
            ar.fit_transient.boundsNeg.lb(ind(i)) = ar.fit_transient.boundsNeg.lb(ind(i))/boundfactor;
        else
            ar.fit_transient.boundsNeg.lb(ind(i)) = ar.fit_transient.boundsNeg.lb(ind(i))*boundfactor;
        end
        if ar.fit_transient.boundsNeg.ub(ind(i))>0
            ar.fit_transient.boundsNeg.ub(ind(i)) = ar.fit_transient.boundsNeg.ub(ind(i))*boundfactor;
        else
            ar.fit_transient.boundsNeg.ub(ind(i)) = ar.fit_transient.boundsNeg.ub(ind(i))/boundfactor;
        end            
    end
    
    ar.lb(ind(i)) = ar.fit_transient.bounds.lb(ind(i));  % use bounds for positive sign as default
    ar.ub(ind(i)) = ar.fit_transient.bounds.ub(ind(i));
end

ind = [ar.fit_transient.indp.offset;ar.fit_transient.indp.timescale_sust;ar.fit_transient.indp.timescale_trans];
for i=1:length(ind)
    ar.p(ind(i)) = mean([ar.fit_transient.bounds.lb(ind(i)),ar.fit_transient.bounds.ub(ind(i))]);
end
ind = [ar.fit_transient.indp.sd;ar.fit_transient.indp.amp_sust;ar.fit_transient.indp.amp_trans];
for i=1:length(ind)
    ar.p(ind(i)) = 0.1*ar.fit_transient.bounds.lb(ind(i))+0.9*ar.fit_transient.bounds.ub(ind(i));  % more towards upper bound
end

ar.p(ar.fit_transient.indp.toffset) = max(0,ar.lb(ar.fit_transient.indp.toffset));
