function arCalcHierarchical(sensi)
%arCalcHierarchical([sensi]) Do the following things:
%   - For each data series with hierarchical scale parameter, assign
%     respective x or z values from the xExpSimu or zExpSimu field in the
%     corresponding condition structure.
%   - Compute hierarchical scale parameters, overriding any custom values
%     concerning these parameters.
%   - For each data series with hierarchical scale parameter, use the newly
%     computed scale value to recompute values of observables stored in the
%     yExpSimu field, overriding the results of arSimu.
%   - If sensi is true, do analogous things for the scale gradients and
%     sensitivities of observables. Namely, compute scale gradients based
%     on sxExpSimu or szExpSimu fields in the condition structures and
%     subsequently recompute syExpSimu fields in the data structures.
%
%   sensi       logical indicating whether to calculate scale gradients   [true]
%               and respectively update sensitivities of observables

if ~exist('sensi','var') || isempty(sensi)
    sensi = 1;
end

global ar

% Check settings which are required to carry out hierarchical optimization
assert(isfield(ar.config,'useHierarchical') && ar.config.useHierarchical, ...
    ['Hierarchical optimization is not enabled. You have not called arInitHierarchical ', ...
     'or there are no scale parameters suitable for hierarchical optimization in your observables.'])
for is = 1:length(ar.scales)
    assert(~ar.qLog10(ar.scales(is).pLink),sprintf('Please disable log10 scale for parameter %s when using hierarchical optimization.',ar.pLabel{ar.scales(is).pLink}))
    assert(ar.qFit(ar.scales(is).pLink)~=1,sprintf('Please disable fitting for parameter %s when using hierarchical optimization.',ar.pLabel{ar.scales(is).pLink}))
end

% Check for features which are yet not supported along with hierarchical optimization (but possibly may be supported in the future)
if ar.config.fiterrors==0
haveNanExpStd = false;
for im = 1:length(ar.model)
    for id = 1:length(ar.model(im).data)
        for iy = 1:length(ar.model(im).data(id).fy)
            if ar.model(im).data(id).useHierarchical(iy) && any(isnan(ar.model(im).data(id).yExpStd(:,iy)))
                haveNanExpStd = true;
            end
        end
    end
end
end
useErrorModel = (ar.config.fiterrors==1 || (ar.config.fiterrors==0 && haveNanExpStd));
assert(~useErrorModel,'Hierarchical optimization in combination with model-based errors is not supported yet. Please use experimental errors for the observables which have parameters for hierarchical optimization.')
useCustomResidual = isfield(ar.config,'user_residual_fun') && ~isempty(ar.config.user_residual_fun);
assert(~useCustomResidual,'Hierarchical optimization in combination with custom residuals is not supported yet.')
assert(~isfield(ar,'conditionconstraints'),'Hierarchical optimization in combination with condition constraints is not supported yet.')
assert(sum(ar.type~=0)==0,'Hierarchical optimization in combination with priors other than flat box is not supported yet.')
assert(~isfield(ar,'random'),'Hierarchical optimization in combination with random effects is not supported yet.')
for im = 1:length(ar.model)
    for ic = 1:length(ar.model(im).condition)
        qssEnabled = isfield(ar.model(im).condition(ic), 'qSteadyState') && sum(ar.model(im).condition(ic).qSteadyState==1)>0;
        assert(~qssEnabled,'Hierarchical optimization in combination with qSteadyState is not supported yet.')
    end
end
for im = 1:length(ar.model)
    for id = 1:length(ar.model(im).data)
        for iy = 1:length(ar.model(im).data(id).fy)
            if ar.model(im).data(id).useHierarchical(iy)
                assert(ar.model(im).data(id).logfitting(iy)==0,'Hierarchical optimization in combination with fitting of observables in log10 scale is not supported yet. Please switch off fitting in log10 scale for the observables which have parameters for hierarchical optimization.')
            end
        end
        if any(ar.model(im).data(id).useHierarchical)
            useDetectionLimit = isfield( ar.model(im).data(id), 'resfunction' ) && isstruct( ar.model(im).data(id).resfunction ) && ar.model(im).data(id).resfunction.active;
            assert(~useDetectionLimit,'Hierarchical optimization in combination with detection limits is not supported yet. Please switch off detection limits for the data sets which have parameters for hierarchical optimization.')
        end
    end
end

% For each data, extract corresponding xExpSimu/zExpSimu values and their sensitivities
for im = 1:length(ar.model)
    for id = 1:length(ar.model(im).data)

        % This check is only for performance - continue early if no reason to stuck in this iteration
        if ~any(ar.model(im).data(id).useHierarchical)
            continue
        end

        ip = ar.model(im).data(id).pCondLink;
        ic = ar.model(im).data(id).cLink;
        td = ar.model(im).data(id).tExp;
        tc = ar.model(im).condition(ic).tExp;

        ar.model(im).data(id).xzExpSimu = zeros(length(td),length(ar.model(im).data(id).fy));
        if sensi
            ar.model(im).data(id).sxzExpSimu = zeros(length(td),length(ar.model(im).data(id).fy),length(ar.model(im).data(id).p));
        end
        for iy = 1:length(ar.model(im).data(id).fy)
            if ar.model(im).data(id).useHierarchical(iy)
                ixz = ar.model(im).data(id).xzLink(iy);
                xzType = ar.model(im).data(id).xzType{iy};
                xzExpSimu = ar.model(im).condition(ic).(sprintf('%sExpSimu',xzType))(:,ixz); % TODO: Test the performance impact of sprintf
                ar.model(im).data(id).xzExpSimu(:,iy) = interp1(tc,xzExpSimu,td,'nearest'); % NOTE: td should always be a subset of tc, hence the 'nearest' option
                if sensi
                    sxzExpSimu = ar.model(im).condition(ic).(sprintf('s%sExpSimu',xzType))(:,ixz,:);
                    ar.model(im).data(id).sxzExpSimu(:,iy,ip) = interp1(tc,sxzExpSimu,td,'nearest');
                end
            end
        end

    end
end

% Compute each scale parameter and their gradients
for is = 1:length(ar.scales)
    num = 0;
    den = 0;
    numGrad = 0;
    denGrad = 0;
    for il = 1:length(ar.scales(is).links)

        im = ar.scales(is).links(il).m;
        id = ar.scales(is).links(il).d;
        iy = ar.scales(is).links(il).fy;
        ystd = ar.model(im).data(id).yExpStd(:,iy);
        yExp = ar.model(im).data(id).yExp(:,iy);
        xzExpSimu = ar.model(im).data(id).xzExpSimu(:,iy);
        num = num + sum(yExp.*xzExpSimu./(ystd.^2));
        den = den + sum(xzExpSimu.^2./(ystd.^2));
        if sensi
            sxzExpSimu = squeeze(ar.model(im).data(id).sxzExpSimu(:,iy,:));
            numGrad = numGrad + sum(yExp.*sxzExpSimu./(ystd.^2),1);
            denGrad = denGrad + 2*sum(xzExpSimu.*sxzExpSimu./(ystd.^2),1);
        end
    end
    assert(den>0,sprintf('The solution of your model corresponding to the scale %s is zero at all measurement times. Hierarchical optimization is not feasible.',ar.scales(is).scaleLabel))
    ar.scales(is).scale = num/den;
    if sensi
        ar.scales(is).scaleGrad = (numGrad*den - denGrad*num)/(den^2);
    end
    ar.p(ar.scales(is).pLink) = ar.scales(is).scale;
end

% Update observables and their sensitivities
for is = 1:length(ar.scales)
    scale = ar.scales(is).scale;
    if sensi
        scaleGrad =  ar.scales(is).scaleGrad;
    end
    for il = 1:length(ar.scales(is).links)
        im = ar.scales(is).links(il).m;
        id = ar.scales(is).links(il).d;
        iy = ar.scales(is).links(il).fy;
        xzExpSimu = ar.model(im).data(id).xzExpSimu(:,iy);
        ar.model(im).data(id).yExpSimu(:,iy) = scale.*xzExpSimu;
        if sensi
            sxzExpSimu = squeeze(ar.model(im).data(id).sxzExpSimu(:,iy,:));
            ar.model(im).data(id).syExpSimu(:,iy,:) = scale.*sxzExpSimu + scaleGrad.*xzExpSimu;
        end
    end
end
% NOTE: Alternatively, we could iterate over the data instead of over the scales in the latter loop.
