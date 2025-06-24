function arDisableHierarchical()
%arDisableHierarchical Restore the original user settings overriden by
%   arInitHierarchical and set ar.config.useHierarchical to false.
%
%   One may say that arInitHierarchical initializes the "hierarchical mode"
%   of optimization and this function switches back to the "normal mode".
%
%   This function in particular attempts to restore the lower and upper
%   bounds for the scale parameters subjected to hierarchical optimization.
%   If the scales computed in the run of hierarchical optimization do not
%   fit the original bounds, then these bounds are adjusted to cover the
%   computed value.
%
%   See also ARINITHIERARCHICAL

global ar

assert(isfield(ar.config,'useHierarchical') && ar.config.useHierarchical, ...
    ['Hierarchical optimization is not enabled. You have not called arInitHierarchical ', ...
     'or there are no scale parameters suitable for hierarchical optimization in your observables.'])

for is = 1:length(ar.scales)
    pLink = ar.scales(is).pLink;
    ar.qFit(pLink) = ar.scales(is).backup_user.qFit;
    ar.qLog10(pLink) = ar.scales(is).backup_user.qLog10;
    ar.lb(pLink) = ar.scales(is).backup_user.lb;
    ar.ub(pLink) = ar.scales(is).backup_user.ub;
    if isnan(ar.p(pLink)) % This happens when this function is called immediately after arInitHierarchical
        ar.p(pLink) = ar.scales(is).backup_user.p;
    elseif ar.qLog10(pLink) % Hierarchical mode always uses linear scale for the scale parameters, so when disabling this mode we apply the log10 scale if relevant
        ar.p(pLink) = log10(ar.p(pLink));
    end
    if ar.lb(pLink) > ar.p(pLink)
        ar.lb(pLink) = ar.lb(pLink) - 0.2*(ar.ub(pLink) - ar.lb(pLink));
    elseif ar.ub(pLink) < ar.p(pLink)
        ar.ub(pLink) = ar.ub(pLink) + 0.2*(ar.ub(pLink) - ar.lb(pLink));
    end
end

ar.config.useHierarchical = false;
