% Sometimes broader bounds are reasonable, i.e. less conservative
% restriction, e.g. if many data points are available. Then boundfactor can
% be set to >1.

function Initialize_FitTransient(boundfactor)
if ~exist('boundfactor','var') || isempty(boundfactor)
    boundfactor = 1;  % standard bounds, designed for experimetnal data    
end
global ar

ar.fit_transient = struct;

name_patterns = struct;
name_patterns.signum = 'signum_TF';
name_patterns.toffset = 'toffset_TF';
name_patterns.offset = 'offset_TF';
name_patterns.amp_sust = 'amp_sust';
name_patterns.amp_trans = 'amp_trans';
name_patterns.timescale_sust = 'timescale_sust';
name_patterns.timescale_trans = 'timescale_trans';
name_patterns.sd = 'sd_TF';

% Find parameters of the TF in ar.pLabel
pnamesTF = fieldnames(name_patterns);
ar.fit_transient.pnameTF = pnamesTF;
ar.fit_transient.indp.all = [];
for ip=1:length(pnamesTF)
    ar.fit_transient.indp.(pnamesTF{ip})    = strmatch(name_patterns.(pnamesTF{ip}),ar.pLabel);
    if isempty(ar.fit_transient.indp.(pnamesTF{ip}))
        warning('Initialize_FitTransient.m: ''%s'' not found.',pnamesTF{ip})
    else
        ar.fit_transient.indp.all = [ar.fit_transient.indp.(pnamesTF{ip})];
    end
    if length(ar.fit_transient.indp.(pnamesTF{ip})) ~= length(ar.fit_transient.indp.(pnamesTF{1}))        
        warning('%s: The number of different parameter-types of the transient function should be equal. \nMay be detection by name does not work. Check parameter names!',pnamesTF{ip});
    end
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
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        vars = symvar(arSym(ar.model(m).data(d).fp));
        if length(ar.model(m).data(d).y)>1
            error('Not yet implemented.')
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



% signums are constants and non-log:
ar.qFit(ar.fit_transient.indp.signum) = 2;
ar.qLog10(ar.fit_transient.indp.signum) = 0;

arQlogParameters(ar.fit_transient.indp.amp_sust, 0);
arQlogParameters(ar.fit_transient.indp.amp_trans, 0);
arQlogParameters(ar.fit_transient.indp.offset, 0);  % often, negative offsets required
arQlogParameters(ar.fit_transient.indp.sd, 1);

% ar.qLog10(ar.fit_transient.indp.amp_sust) = 0;
% ar.qLog10(ar.fit_transient.indp.amp_trans) = 0;
% ar.qLog10(ar.fit_transient.indp.offset) = 0;  % often, negative offsets required
% ar.qLog10(ar.fit_transient.indp.sd) = 1;


for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        if sum(ar.model(m).data(d).logfitting~=0)>0
            error('The transient function is design for non-log fitting. If log-trans required, do it before/do it by hand.')       
        elseif sum(ar.model(m).data(d).logplotting~=0)>0
            error('The transient function is design for non-log fitting. If log-trans required, do it before/do it by hand.')       
        end
    end
end


[ar.fit_transient.bounds,ar.fit_transient.boundsNeg] = DefaultLbUbTransient;

ind = ar.fit_transient.indp.all;
for i=1:length(ind)
    if ar.qLog10(ind(i))==1
        ar.fit_transient.bounds.lb(ind(i)) = ar.fit_transient.bounds.lb(ind(i))-log10(boundfactor);
        ar.fit_transient.bounds.ub(ind(i)) = ar.fit_transient.bounds.ub(ind(i))+log10(boundfactor);
        ar.fit_transient.boundsNeg.lb(ind(i)) = ar.fit_transient.boundsNeg.lb(ind(i))-log10(boundfactor);
        ar.fit_transient.boundsNeg.ub(ind(i)) = ar.fit_transient.boundsNeg.ub(ind(i))+log10(boundfactor);
    else
        ar.fit_transient.bounds.lb(ind(i)) = ar.fit_transient.bounds.lb(ind(i))/boundfactor;
        ar.fit_transient.bounds.ub(ind(i)) = ar.fit_transient.bounds.ub(ind(i))*boundfactor;
        ar.fit_transient.boundsNeg.lb(ind(i)) = ar.fit_transient.boundsNeg.lb(ind(i))/boundfactor;
        ar.fit_transient.boundsNeg.ub(ind(i)) = ar.fit_transient.boundsNeg.ub(ind(i))*boundfactor;
    end
end

