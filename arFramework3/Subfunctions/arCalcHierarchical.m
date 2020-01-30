function arCalcHierarchical(sensi)
%arCalcHierarchical Do the following things:
%   - For each data series with hierarchical scale parameter, assign
%     respective x or z values from the xExpSimu or zExpSimu field in the
%     corresponding condition structure.
%   - Compute hierarchical scale parameters, overriding any custom values
%     concerning these parameters.
%   - For each data series with hierarchical scale parameter, recompute
%     values of observables using the newly computed scale values,
%     overriding the results of arSimu.
%   - Do analogous things for the scales gradients and observables of
%     sensitivities, based in sxExpSimu or szExpSimu fields in the
%     condition structures.
%
%   The function takes sensi as an argument to expose interfece coherent
%   with the convention of other d2d functions. However, this argument
%   is not used currently. The function always computes scale gradients
%   and sensitivities of observables, regardless of the value of sensi.
%   Put simply, if sxExpSimu or szExpSimu values are not up to date, then
%   the resulting scale gradients and sensitivities of observables are not
%   up to date. This is actually analogous to what happens with sxExpSimu
%   and szExpSimu when calling arSimu with sensi set to false.

if ~exist('sensi','var') || isempty(sensi)
    sensi = 1;
end

global ar

assert(isfield(ar.config,'useHierarchical') && ar.config.useHierarchical, ...
    ['Hierarchical optimization is not enabled. You have not called arInitHierarchical ', ...
     'or there are no scale parameters suitable for hierarchical optimization in your observables.'])
errorFitting = ( ar.config.fiterrors == 1) || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)==1)>0 );
assert(~errorFitting,'Hierarchical optimization in combination with with error fitting is not supported yet.')
assert(sum(ar.type~=0)==0,'Hierarchical optimization in combination with priors other than flat box is not supported yet.')
for im = 1:length(ar.model)
    for id = 1:length(ar.model(im).data)
        useCustomResidual = isfield( ar.model(im).data(id), 'resfunction' ) && isstruct( ar.model(im).data(id).resfunction ) && ar.model(im).data(id).resfunction.active;
        assert(~useCustomResidual,'Please choose between hierarchical optimization and custom residuals.')
    end
end
for is = 1:length(ar.scales)
    assert(~ar.qLog10(ar.scales(is).pLink),sprintf('Please disable log10 scale for parameter %s when using hierarchical optimization.',ar.pLabel{ar.scales(is).pLink}))
    assert(ar.qFit(ar.scales(is).pLink)~=1,sprintf('Please disable fitting for parameter %s when using hierarchical optimization.',ar.pLabel{ar.scales(is).pLink}))
end

% For each data, extract corresponding xExpSimu/zExpSimu values and their sensitivities
for im = 1:length(ar.model)
    for id = 1:length(ar.model(im).data)

        % This check is only for performance - continue early if no reason to stuck in this iteration
        if all(isnan(ar.model(im).data(id).useHierarchical))
            continue
        end

        ip = ar.model(im).data(id).pCondLink;
        ic = ar.model(im).data(id).cLink;
        td = ar.model(im).data(id).tExp;
        tc = ar.model(im).condition(ic).tExp;
        it = ismember(tc,td); % TODO: Test the performance of this function and possibly consider another solution

        ar.model(im).data(id).xzExpSimu = zeros(sum(it),length(ar.model(im).data(id).fy));
        ar.model(im).data(id).sxzExpSimu = zeros(sum(it),length(ar.model(im).data(id).fy),length(ar.model(im).data(id).p));
        for iy = 1:length(ar.model(im).data(id).fy)
            if ar.model(im).data(id).useHierarchical(iy)
                ixz = ar.model(im).data(id).xzLink(iy);
                xzType = ar.model(im).data(id).xzType{iy};
                ar.model(im).data(id).xzExpSimu(:,iy) = ar.model(im).condition(ic).(sprintf('%sExpSimu',xzType))(it,ixz); % TODO: Test the performance impact of sprintf
                ar.model(im).data(id).sxzExpSimu(:,iy,ip) = ar.model(im).condition(ic).(sprintf('s%sExpSimu',xzType))(it,ixz,:);
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
        if ar.config.fiterrors == -1
            ystd = ar.model(im).data(id).yExpStd(:,iy);
        elseif ar.config.fiterrors == 0
            ystd = ar.model(im).data(id).yExpStd(:,iy);
            noSD = isnan(ystd);
            ystd(noSD) = ar.model(im).data(id).ystdExpSimu(noSD,iy);
        end
        yExp = ar.model(im).data(id).yExp(:,iy);
        xzExpSimu = ar.model(im).data(id).xzExpSimu(:,iy);
        sxzExpSimu = squeeze(ar.model(im).data(id).sxzExpSimu(:,iy,:));

        num = num + sum(yExp.*xzExpSimu./(ystd.^2));
        den = den + sum(xzExpSimu.^2./(ystd.^2));

        numGrad = numGrad + sum(yExp.*sxzExpSimu./(ystd.^2),1);
        denGrad = denGrad + 2*sum(xzExpSimu.*sxzExpSimu./(ystd.^2),1);
    end
    ar.scales(is).scale = num/den;
    ar.scales(is).scaleGrad = (numGrad*den - denGrad*num)/(den^2);
    ar.p(ar.scales(is).pLink) = ar.scales(is).scale;
end

% Update observables and their sensitivities
for is = 1:length(ar.scales)
    scale = ar.scales(is).scale;
    scaleGrad =  ar.scales(is).scaleGrad;
    for il = 1:length(ar.scales(is).links)
        im = ar.scales(is).links(il).m;
        id = ar.scales(is).links(il).d;
        iy = ar.scales(is).links(il).fy;
        xzExpSimu = ar.model(im).data(id).xzExpSimu(:,iy);
        sxzExpSimu = squeeze(ar.model(im).data(id).sxzExpSimu(:,iy,:));

        ar.model(im).data(id).yExpSimu(:,iy) = scale.*xzExpSimu;
        ar.model(im).data(id).syExpSimu(:,iy,:) = scale.*sxzExpSimu + scaleGrad.*xzExpSimu;
    end
end
% NOTE: Alternatively, we could iterate over the data instead of over the scales in the latter loop.

% TODO: This function should update also yFineSimu.
