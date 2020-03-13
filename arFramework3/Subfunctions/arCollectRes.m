% arCollectRes(sensi, [debugres])
%
% Collects all residuals, sres and chi2 of the individual data sets and 
% calculates the number of data points. Additional the priors, constr,
% random effects and user defined residuals are collected.
% 
%   sensi          boolean, collect sensitivities
%   debugres  [0]  boolean, fill ar.resinfo with 
%                  additional information 
%
% Function collects residuals and chi2 calculated by arCalcRes(true)
% from the data structs to the top level of the ar struct  
%   - ar.model.data.res         -> ar.res , ar.type=1
%   - ar.model.data.reserr      -> ar.res , ar.type=2
%   - prior                     -> ar.res , ar.type=3
%   - constr                    -> ar.constr
%   - random                    -> ar.res , ar.type=4
%   - ar.model.data.sres        -> ar.sres 
%   - ar.model.data.chi2        -> ar.chi2


function arCollectRes(sensi, debugres)

global ar 

if ( nargin < 2 )
    debugres = 0;
end
if ~isfield(ar,'ndata_res')
    arCalcRes(true)
end
ar.ndata = ar.ndata_res;
ar.nprior = 0;
ar.nrandom = 0;
ar.nconstr = 0;

ar.chi2 = 0;
ar.chi2err = 0;
ar.chi2prior = 0;
ar.chi2random = 0;
ar.chi2constr = 0;
ar.chi2fit = 0;
ar.chi2err_logdataCorrection = 0;

ar.firstorderopt = nan;

if(~isfield(ar,'res'))
    ar.res = [];
end
if(~isfield(ar,'res_type'))  
    ar.res_type = [];
end
if(~isfield(ar,'sres'))
    ar.sres = [];
end
if(~isfield(ar,'constr'))
    ar.constr = [];
end
if(~isfield(ar,'sconstr'))
    ar.sconstr = [];
end
if(~isfield(ar,'L1subtype'))
    ar.L1subtype = ones(size(ar.p));
end

np = length(ar.p);

useMSextension = false;

resindex = 1;
sresindex = 1;

% fit errors?
fiterrors = ( ar.config.fiterrors == 1  || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0 ) );

ar.sres_size = []; % for NEB merger

% data
for jm = 1:length(ar.model)
    if(isfield(ar.model(jm), 'data'))
        nd = length(ar.model(jm).data);
        for jd = 1:nd
            if(ar.model(jm).data(jd).has_yExp)
                ar.chi2 = ar.chi2 + sum(ar.model(jm).data(jd).chi2(ar.model(jm).data(jd).qFit==1));
                
                if(useMSextension && isfield(ar, 'ms_count_snips') && ~isempty(ar.model(jm).data(jd).ms_index))
                    
                else
                    % collect residuals for fitting
                    tmpres = ar.model(jm).data(jd).res(:,ar.model(jm).data(jd).qFit==1);
                    ar.res(resindex:(resindex+length(tmpres(:))-1)) = tmpres;
                    ar.res_type(resindex:(resindex+length(tmpres(:))-1)) = 1;
                    
                    if ( debugres )
                        for jr = 1 : numel( tmpres )
                            ar.resinfo(resindex + jr - 1).type = 'data';
                            ar.resinfo(resindex + jr - 1).m = jm;
                            ar.resinfo(resindex + jr - 1).d = jd;
                        end
                    end
                    resindex = resindex+length(tmpres(:));
                    
                    if( fiterrors )                        
                        ar.chi2err = ar.chi2err + sum(ar.model(jm).data(jd).chi2err(ar.model(jm).data(jd).qFit==1));
                        tmpreserr = ar.model(jm).data(jd).reserr(:,ar.model(jm).data(jd).qFit==1);
                        ar.res(resindex:(resindex+length(tmpreserr(:))-1)) = tmpreserr;
                        ar.res_type(resindex:(resindex+length(tmpreserr(:))-1)) = 2;
                        if ( debugres )
                            for jr = 1 : numel( tmpreserr )
                                ar.resinfo(resindex + jr - 1).type = 'error model';
                                ar.resinfo(resindex + jr - 1).m = jm;
                                ar.resinfo(resindex + jr - 1).d = jd;
                            end
                        end
                        resindex = resindex+length(tmpreserr(:));
                        
                        % add correction term for fitting on log-axis due
                        % to normalization term of log-normal distribution
                        islogfitted = logical(ar.model(jm).data(jd).logfitting);
                        if any(islogfitted)
                            data = ar.model(jm).data(jd).yExp(islogfitted);
                            data = data(~isnan(data));
                            ar.chi2err_logdataCorrection = ar.chi2err_logdataCorrection+ sum(log(log(10)^2 * (10.^data).^2));
                        end

                    end
                    
                    % add correction term for logarithmic fitting
                    if any(ar.model(jm).data(jd).logfitting)
                    end
                    
                    % collect sensitivities for fitting
                    if(ar.config.useSensis && sensi)
                        tmptmpsres = ar.model(jm).data(jd).sres(:,ar.model(jm).data(jd).qFit==1,:);
                        tmpsres = zeros(length(tmpres(:)), np);
                        tmpsres(:,ar.model(jm).data(jd).pLink) = reshape(tmptmpsres, ...
                            length(tmpres(:)), sum(ar.model(jm).data(jd).pLink));
                        ar.sres(sresindex:(sresindex+length(tmpres(:))-1),:) = tmpsres;
                        sresindex = sresindex+length(tmpres(:));
                        
                        if ( fiterrors )
                            tmpsreserr = zeros(length(tmpreserr(:)), np);
                            tmpsreserr(:,ar.model(jm).data(jd).pLink) = reshape(ar.model(jm).data(jd).sreserr(:,ar.model(jm).data(jd).qFit==1,:), ...
                                length(tmpreserr(:)), sum(ar.model(jm).data(jd).pLink));
                            ar.sres(sresindex:(sresindex+length(tmpres(:))-1),:) = tmpsreserr;
                            sresindex = sresindex+length(tmpres(:));
                        end
                    end
                end
            end
        end
    end
    
    %  NEB merger
    if(isfield(ar, 'merger'))
        if(isfield(ar.merger, 'neb'))
            if(isfield(ar.merger.neb, 'state') && strcmp(ar.merger.neb.state,'on'))

                ar.sres_size(jm) = sresindex;

                if jm > 1 && jm < length(ar.model)
                    r = ar.p(ar.merger.neb.ps_index(jm+1,:)) - ar.p(ar.merger.neb.ps_index(jm-1,:));
                    grad = ar.sres([ar.sres_size(jm-1):ar.sres_size(jm)-1],ar.merger.neb.ps_index(jm,:));

                    proj = nan(size(grad));
                    for ig = 1: size(grad,1)
                        proj(ig,:) = (sum(grad(ig,:) .* r)./norm(r) ).* r ./norm(r);
                    end

                    ar.sres([ar.sres_size(jm-1):ar.sres_size(jm)-1] ,ar.merger.neb.ps_index(jm,:)) ...
                        = grad-proj;

                end
            end
        end
    end

end

% constraints
constrindex = 1;
sconstrindex = 1;

% Intercondition constraints
% To do: Further generalize this to arbitrary constraints
if ( isfield( ar, 'conditionconstraints' ) )
    for jC = 1 : length( ar.conditionconstraints )
        m1                  = ar.conditionconstraints(jC).m1;
        m2                  = ar.conditionconstraints(jC).m2;
        relative            = ar.conditionconstraints(jC).relative;
        c1                  = ar.conditionconstraints(jC).c1;
        c2                  = ar.conditionconstraints(jC).c2;
        tLink1              = ar.conditionconstraints(jC).tLink1;
        tLink2              = ar.conditionconstraints(jC).tLink2;
        sd                  = ar.conditionconstraints(jC).sd;
        states1             = ar.conditionconstraints(jC).states1;
        states2             = ar.conditionconstraints(jC).states2;
        pLink1              = ar.model(m1).condition(c1).pLink;
        pLink2              = ar.model(m2).condition(c2).pLink;

        nstates             = length(states1);
        npts                = length(tLink1);
        tmpsres             = zeros(npts*nstates, np);

        % Fetch simulations
        dynamic1            = ar.model(m1).condition(c1).xExpSimu(tLink1,states1);
        dynamic2            = ar.model(m2).condition(c2).xExpSimu(tLink2,states2);

        % Fetch sensitivities w.r.t. p
        sens1               = ar.model(m1).condition(c1).sxExpSimu(tLink1,states1,:);
        sens2               = ar.model(m2).condition(c2).sxExpSimu(tLink2,states2,:);

        % Relative?
        if ( relative )            
            sens1           = sens1 ./ ( log(10) * repmat(dynamic1,1,1,size(sens1,3)) );
            sens2           = sens2 ./ ( log(10) * repmat(dynamic2,1,1,size(sens1,3)) );

            dynamic1        = log10( dynamic1 );
            dynamic2        = log10( dynamic2 );            
        end

        sens1               = arTrafoParameters(sens1,m1,c1,false);
        sens2               = arTrafoParameters(sens2,m2,c2,false);

        % Reshape to fit format desired for sres
        sens1               = reshape( sens1, nstates*npts, sum(pLink1));
        sens2               = reshape( sens2, nstates*npts, sum(pLink2));
        tmpsres(:, pLink1)  = tmpsres(:, pLink1) + sens1 / sd;
        tmpsres(:, pLink2)  = tmpsres(:, pLink2) - sens2 / sd;
        tmpres              = (dynamic1 - dynamic2)./sd;

        % Store
        %ar.res(resindex:resindex+npts*nstates-1) = tmpres;
        %ar.sres(sresindex:(sresindex+npts*nstates-1),:) = tmpsres;
        %sresindex = sresindex+npts*nstates;

        ar.constr(constrindex:constrindex+npts*nstates-1) = tmpres;
        ar.sconstr(sconstrindex:(sconstrindex+npts*nstates-1),:) = tmpsres;

        constrindex  = constrindex+npts*nstates;
        sconstrindex = sconstrindex+npts*nstates;

        ar.nconstr      = ar.nconstr + length(tmpres);
        ar.chi2constr   = ar.chi2constr + sum(sum(tmpres.^2));        
    end
end

% steady state conditions
qRelativeToInitialValue = true;
for jm = 1:length(ar.model)
    nc = length(ar.model(jm).condition);
    for jc = 1:nc
        if(isfield(ar.model(jm).condition(jc), 'qSteadyState'))
            qss = ar.model(jm).condition(jc).qSteadyState==1;
            if(sum(qss)>0)
                if(qRelativeToInitialValue)
                    if(isfield(ar.model(jm).condition(jc), 'xExpSimu') && ...
                            ~isempty(ar.model(jm).condition(jc).xExpSimu))
                        x = ar.model(jm).condition(jc).xExpSimu(1,qss);
                    else
                        x = ar.model(jm).condition(jc).xFineSimu(1,qss);
                    end
                    dxdt = ar.model(jm).condition(jc).dxdt(qss);
                    tmpStdSteadystate = ar.model(jm).condition(jc).stdSteadyState(qss);
                    tmpconstr = (dxdt ./ x) ./ tmpStdSteadystate;
                    validconstr = ~(isnan(tmpconstr) | isinf(tmpconstr));
                    tmpconstr = tmpconstr(validconstr);
                    ar.constr(constrindex:(constrindex+length(tmpconstr(:))-1)) = tmpconstr;
                    constrindex = constrindex+length(tmpconstr(:));
                    ar.nconstr = ar.nconstr + sum(qss);
                    ar.chi2constr = ar.chi2constr + sum(tmpconstr.^2);
                    
                    if(ar.config.useSensis && sensi)
                        tmpsconstr = zeros(length(tmpconstr(:)), np);
                        if(isfield(ar.model(jm).condition(jc), 'sxExpSimu') && ...
                            ~isempty(ar.model(jm).condition(jc).sxExpSimu))
                            dxdp = squeeze(ar.model(jm).condition(jc).sxExpSimu(1,qss,:));
                            if(iscolumn(dxdp)) % transpose dxdp if squeeze returns column vector (sum(qss)==1)
                                dxdp = dxdp';
                            end
                            dxdp = arTrafoParameters(dxdp,jm,jc,false);
                        else
                            dxdp = squeeze(ar.model(jm).condition(jc).sxFineSimu(1,qss,:));
                            if(iscolumn(dxdp)) % transpose dxdp if squeeze returns column vector (sum(qss)==1)
                                dxdp = dxdp';
                            end
                        end
                        ddxdtdp = ar.model(jm).condition(jc).ddxdtdp(qss,:);
                        ddxdtdp = arTrafoParameters(ddxdtdp,jm,jc,false);
                        
                        tmptmpsconstr = bsxfun(@rdivide,ddxdtdp,x') - bsxfun(@times,bsxfun(@rdivide,dxdt,x.^2)', dxdp);
                        tmpsconstr(:,ar.model(jm).condition(jc).pLink) = tmptmpsconstr(validconstr,:);
                        tmpsconstr = tmpsconstr ./ (tmpStdSteadystate(validconstr)'*ones(1,np));
                        tmpsconstr(:,ar.qFit~=1) = 0;
                        
                        ar.sconstr(sconstrindex:(sconstrindex+length(tmpconstr(:))-1),:) = tmpsconstr;
                        sconstrindex = sconstrindex+length(tmpconstr(:));
                    end
                else
                    tmpconstr = ar.model(jm).condition(jc).dxdt(qss) ./ ar.model(jm).condition(jc).stdSteadyState(qss);
                    ar.constr(constrindex:(constrindex+length(tmpconstr(:))-1)) = tmpconstr;
                    constrindex = constrindex+length(tmpconstr(:));
                    ar.nconstr = ar.nconstr + sum(qss);
                    ar.chi2constr = ar.chi2constr + sum(tmpconstr.^2);
                    
                    if(ar.config.useSensis && sensi)
                        tmpsconstr = zeros(length(tmpconstr(:)), np);
                        ddxdtdp = ar.model(jm).condition(jc).ddxdtdp(qss,:);
                        ddxdtdp = arTrafoParameters(ddxdtdp,jm,jc,false);
                        ddxdtdp = bsxfun(@rdivide, ddxdtdp, ar.model(jm).condition(jc).stdSteadyState(qss)');
                        
                        tmpsconstr(:,ar.model(jm).condition(jc).pLink) = ddxdtdp;
                        tmpsconstr(:,ar.qFit~=1) = 0;
                       
                        ar.sconstr(sconstrindex:(sconstrindex+length(tmpconstr(:))-1),:) = tmpsconstr;
                        sconstrindex = sconstrindex+length(tmpconstr(:));
                    end
                end
            end
        end
    end
end

% priors
notpriors = 1:sresindex-1;
for jp=1:np
    if(ar.type(jp) == 0) % flat prior with hard bounds
    elseif(ar.type(jp) == 1) % normal prior
        tmpres = (ar.mean(jp)-ar.p(jp))./ar.std(jp);
        ar.res(resindex) = tmpres;
        ar.res_type(resindex) = 3; 
        if ( debugres )
            ar.resinfo(resindex).type = 'normal prior';
            ar.resinfo(resindex).m = jp;
        end
        resindex = resindex+1;
        if(ar.config.useSensis && sensi)
            tmpsres = zeros(size(ar.p));
            tmpsres(jp) = -1 ./ ar.std(jp);
            ar.sres(sresindex,:) = tmpsres;
            sresindex = sresindex+1;
        end
        ar.ndata = ar.ndata + 1;
        ar.nprior = ar.nprior + 1;
        ar.chi2 = ar.chi2 + tmpres^2;
        ar.chi2prior = ar.chi2prior + tmpres^2;
    elseif(ar.type(jp) == 2) % uniform with normal bounds
        tmpres = 0;
        w = 0.1;
        if(ar.p(jp) < ar.lb(jp))
            tmpres = (ar.p(jp) - ar.lb(jp))*w;
        elseif(ar.p(jp) > ar.ub(jp))
            tmpres = (ar.p(jp) - ar.ub(jp))*w;
        end
        ar.res(resindex) = tmpres;
        ar.res_type(resindex) = 3; 
        if ( debugres )
            ar.resinfo(resindex).type = 'uniform with normal bounds';
            ar.resinfo(resindex).m = jp;
        end
        resindex = resindex+1;
        if(ar.config.useSensis && sensi)
            tmpsres = zeros(size(ar.p));
            tmpsres(jp) = w;
            ar.sres(sresindex,:) = tmpsres;
            sresindex = sresindex+1;
        end
        ar.ndata = ar.ndata + 1;
        ar.nprior = ar.nprior + 1;
        ar.chi2 = ar.chi2 + tmpres^2;
        ar.chi2prior = ar.chi2prior + tmpres^2;
    elseif(ar.type(jp) == 3 && ~isinf(ar.std(jp)) ) % Lasso and Adaptions
        
        threshTo0 = 1e-10;
        switch ar.L1subtype(jp)
            case 1 % L1
                tmpres = sqrt(abs(ar.mean(jp)-ar.p(jp))./abs(ar.std(jp)));
                if(ar.config.useSensis && sensi)
                    tmpsres = zeros(size(ar.p));
                    if abs(ar.mean(jp) - ar.p(jp)) > threshTo0
                        tmpsres(jp) = sign(ar.p(jp) - ar.mean(jp)) ...
                            .* 0.5./sqrt(abs(ar.mean(jp) - ar.p(jp))...
                            .*ar.std(jp));
                    elseif abs(2*ar.res(notpriors)*ar.sres(notpriors,jp)) < ...
                            (1./ar.std(jp))
                        ar.sres(:,jp) = 0;
                    end
                end
            case 2 % Adaptive Lasso
                tmpres = sqrt(ar.lnuweights(jp).*abs(ar.mean(jp)-ar.p(jp))...
                    ./abs(ar.std(jp)));
                if(ar.config.useSensis && sensi)
                    tmpsres = zeros(size(ar.p));
                    if abs(ar.mean(jp) - ar.p(jp)) > threshTo0
                        tmpsres(jp) = ar.lnuweights(jp).*sign(ar.p(jp) - ar.mean(jp)) ...
                            .* 0.5./sqrt(abs(ar.mean(jp) - ar.p(jp))...
                            .*ar.std(jp));
                    elseif abs(2*ar.res(notpriors)*ar.sres(notpriors,jp)) < ...
                            (1./ar.std(jp)* ar.lnuweights(jp))
                        ar.sres(:,jp) = 0;
                    end
                end
            case 3 % Lq
                tmpres = sqrt(abs(ar.mean(jp)-ar.p(jp))...
                    .^(ar.expo(jp))./abs(ar.std(jp)));
                if(ar.config.useSensis && sensi)
                    tmpsres = zeros(size(ar.p));
                    if abs(ar.mean(jp) - ar.p(jp)) > threshTo0
                        tmpsres(jp) = sign(ar.p(jp) - ar.mean(jp)) ...
                            .* 0.5.*(ar.expo(jp)) .* sqrt(abs(ar.mean(jp) - ar.p(jp))...
                            .^(ar.expo(jp)-2) ./ar.std(jp));
                    elseif abs(2*ar.res(notpriors)*ar.sres(notpriors,jp)) < ...
                            (1./ar.std(jp)* ar.expo(jp) ...
                            * threshTo0 .^(ar.expo(jp)-1))
                        ar.sres(:,jp) = 0;
                    end
                end
            case 4 % Elastic Net
                tmpres = sqrt(ar.alpha(jp) ./ ar.std(jp)) ...
                    .* (ar.mean(jp)-ar.p(jp));
                % quadratic contribution
                ar.res(resindex) = tmpres;
                resindex = resindex + 1;
                tmpres = sqrt((1-ar.alpha(jp)) ./ ar.std(jp) ...
                    .* abs(ar.mean(jp)-ar.p(jp)));
                % L1 contribution
                if(ar.config.useSensis && sensi)
                    tmpsres = zeros(size(ar.p));
                    tmpsres(jp) = sqrt(ar.alpha(jp)./ar.std(jp));
                    ar.sres(sresindex,:) = tmpsres;
                    sresindex = sresindex + 1;
                    % quadratic sensitivities are just constants
                    
                    tmpsres = zeros(size(ar.p));
                    if abs(ar.mean(jp) - ar.p(jp)) > threshTo0
                        tmpsres(jp) = sign(ar.p(jp) - ar.mean(jp)) .* 0.5 ...
                            .* sqrt((1-ar.alpha(jp)) ./ ar.std(jp) ...
                            ./ abs(ar.p(jp) - ar.mean(jp)));                        
                    elseif abs(2*ar.res(notpriors)*ar.sres(notpriors,jp)) < ...
                            (1-ar.alpha(jp)) ./ ar.std(jp)
                        ar.sres(:,jp) = 0;
                    end
                end
                    
        end
        
        ar.res(resindex) = tmpres;
        ar.res_type(resindex) = 3; 
        resindex = resindex+1;
        if(ar.config.useSensis && sensi)
            ar.sres(sresindex,:) = tmpsres;
            sresindex = sresindex+1;
        end
        ar.ndata = ar.ndata + 1;
        ar.nprior = ar.nprior + 1;
        ar.chi2 = ar.chi2 + tmpres^2;
        ar.chi2prior = ar.chi2prior + tmpres^2;
    end
end

% grouped priors

if any(ar.type == 5)
    
    indsp = 1:np;
    
    for g = ar.grplas.groups
        
%         w = ar.grplas.weights(g);
        gind = indsp((ar.grplas.grouping == g) & (ar.type == 5) );
        % indices of parameters grouped as g and marked as group lasso
        
%         tmpsum = sqrt(sum((ar.mean(gind) - ar.p(gind)).^2 ...
%            ./(ar.std(gind).^2), 2));
        
        tmpsum = sqrt(...
              (ar.mean(gind) - ar.p(gind))./ar.std(gind) ...
            * ar.grplas.A(gind,gind)...
            * ((ar.mean(gind) - ar.p(gind))./ar.std(gind))');
        
        
        
        tmpres = sqrt(tmpsum);
        ar.res(resindex) = tmpres;
        ar.res_type(resindex) = 3;
        resindex = resindex + 1;
        
        if(ar.config.useSensis && sensi)
            tmpsres = zeros(size(ar.p));
            bA = (ar.mean(gind) - ar.p(gind)) ./ (ar.std(gind).^2) ...
                * (ar.grplas.A(gind,gind) + ar.grplas.A(gind,gind)');
            if tmpsum > 1e-10
                tmpsres(gind) = - 1/4 * bA * (tmpsum .^ (-3/2));
            else
                Aii = diag(ar.grplas.A(gind,gind));
                exc = abs(2*ar.res(notpriors)...
                    * ar.sres(notpriors,gind)) < sqrt(Aii)' ./ ar.std(gind);
                ar.sres(:,gind(exc)) = 0;
            end
            ar.sres(sresindex,:) = tmpsres;
            sresindex = sresindex+1;
        end
        ar.ndata = ar.ndata + 1;
        ar.nprior = ar.nprior + 1;
        ar.chi2 = ar.chi2 + tmpres^2;
        ar.chi2prior = ar.chi2prior + tmpres^2;
    end
    
end

% random effects
if(isfield(ar, 'random'))
    for j=1:length(ar.random)
        [tmpres, tmpsres] = arRandomEffect(ar.p(ar.random{j}));
        ar.res(resindex) = tmpres;
        ar.res_type(resindex) = 5; 
        resindex = resindex+1;
        if(ar.config.useSensis && sensi)
            tmpsres2 = zeros(size(ar.p));
            tmpsres2(ar.random{j}) = tmpsres;
            ar.sres(sresindex,:) = tmpsres2;
            sresindex = sresindex+1;
        end
        ar.ndata = ar.ndata + 1;
        ar.nrandom = ar.nrandom + 1;
        ar.chi2 = ar.chi2 + tmpres^2;
        ar.chi2random = ar.chi2random + tmpres^2;        
    end
end

% % multiple shooting (difference penality)
% if(isfield(ar,'ms_count_snips'))
%     ar.ms_violation = 0;
%     for jm = 1:length(ar.model)
%         for jms = 1:size(ar.model(jm).ms_link,1)
%             tmpres1 = ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).xExpSimu(ar.model(jm).ms_link(jms,4),:);
%             tmpres2 = ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).xExpSimu(ar.model(jm).ms_link(jms,5),:);
%             ar.ms_violation = [ar.ms_violation (log10(tmpres1) - log10(tmpres2)).^2];
%             
%             if(ar.ms_strength>0)
%                 tmpres = sqrt(ar.ms_strength) * (tmpres1 - tmpres2);
%                 ar.res(resindex:(resindex+length(tmpres(:))-1)) = tmpres;
%                 resindex = resindex+length(tmpres(:));
%                 
%                 if(ar.config.useSensis && sensi)
%                     tmpsres1 = zeros(length(tmpres(:)), np);
%                     tmpsres1(:,ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).pLink) = ...
%                         reshape(ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).sxExpSimu(...
%                         ar.model(jm).ms_link(jms,4),:,:), length(tmpres(:)), ...
%                         sum(ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).pLink));
%                     
%                     tmpsres2 = zeros(length(tmpres(:)), np);
%                     tmpsres2(:,ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).pLink) = ...
%                         reshape(ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).sxExpSimu(...
%                         ar.model(jm).ms_link(jms,5),:,:), length(tmpres(:)), ...
%                         sum(ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).pLink));
%                     
%                     tmpsres = sqrt(ar.ms_strength) * (tmpsres1 - tmpsres2);
%                     
%                     for j=find(ar.qLog10==1)
%                         tmpsres(:,j) = tmpsres(:,j) * 10.^ar.p(j) * log(10);
%                     end
%                     
%                     ar.sres(sresindex:(sresindex+length(tmpres(:))-1),:) = tmpsres;
%                     sresindex = sresindex+length(tmpres(:));
%                 end
%             end
%         end
%     end
% end

% % multiple shooting (log ratio penalty)
% if(isfield(ar,'ms_count_snips'))
%     ar.ms_violation = [];
%     for jm = 1:length(ar.model)
%         for jms = 1:size(ar.model(jm).ms_link,1)
%             tmpres1 = ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).xExpSimu(ar.model(jm).ms_link(jms,4),:);
%             tmpres2 = ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).xExpSimu(ar.model(jm).ms_link(jms,5),:);
%             ar.ms_violation = [ar.ms_violation (log10(tmpres1) - log10(tmpres2)).^2];
%             
%             if(ar.ms_strength>0)
%                 tmpres = sqrt(ar.ms_strength) * (log10(tmpres1) - log10(tmpres2));
%                 ar.res(resindex:(resindex+length(tmpres(:))-1)) = tmpres;
%                 resindex = resindex+length(tmpres(:));
%                 
%                 if(ar.config.useSensis && sensi)
%                     tmpsres1 = zeros(length(tmpres(:)), np);
%                     tmpsres1(:,ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).pLink) = ...
%                         reshape(ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).sxExpSimu(...
%                         ar.model(jm).ms_link(jms,4),:,:), length(tmpres(:)), ...
%                         sum(ar.model(jm).condition(ar.model(jm).ms_link(jms,1)).pLink));
%                     
%                     tmpsres2 = zeros(length(tmpres(:)), np);
%                     tmpsres2(:,ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).pLink) = ...
%                         reshape(ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).sxExpSimu(...
%                         ar.model(jm).ms_link(jms,5),:,:), length(tmpres(:)), ...
%                         sum(ar.model(jm).condition(ar.model(jm).ms_link(jms,2)).pLink));
%                     
%                     tmpsres = sqrt(ar.ms_strength) * (tmpsres1./(tmpres1'*ones(1,np)) - (tmpsres2./(tmpres2'*ones(1,np)))) / log(10);
%                     
%                     for j=find(ar.qLog10==1)
%                         tmpsres(:,j) = tmpsres(:,j) * 10.^ar.p(j) * log(10);
%                     end
%                     
%                     ar.sres(sresindex:(sresindex+length(tmpres(:))-1),:) = tmpsres;
%                     sresindex = sresindex+length(tmpres(:));
%                 end
%             end
%         end
%     end
% end

%% user-defined residuals (calculated by ar.config.user_residual_fun)
if isfield(ar.res_user,'res')
    for i=1:length(ar.res_user.res)
        ar.res(resindex) = ar.res_user.res(i);
        ar.res_type(resindex) = ar.res_user.type(i);
        if ( debugres )
            ar.resinfo(resindex).type = 'user residual';
        end
        resindex = resindex + 1;
    end
    ar.chi2 = ar.chi2 + sum(ar.res_user.res.^2);
    ar.ndata = ar.ndata + length(ar.res_user.res);

    if ( sensi )
        for i=1:size(ar.res_user.sres,1)
            ar.sres(sresindex,:) = ar.res_user.sres(i,:);
            sresindex = sresindex + 1;
        end
    end
end


% cut off too long arrays
if(isfield(ar.model, 'data'))
    if(~isempty(ar.res))
        if(length(ar.res)>=resindex)
            ar.res(resindex:end) = [];
        end
        if(length(ar.res_type)>=resindex)
            ar.res_type(resindex:end) = [];
        end
        if(ar.config.useSensis && sensi && size(ar.sres,1)>=sresindex)
            ar.sres(sresindex:end,:) = [];
        end
    end
    if(~isempty(ar.constr))
        if(length(ar.constr)>=constrindex)
            ar.constr(constrindex:end) = [];
        end
        if(ar.config.useSensis && sensi && size(ar.sconstr,1)>=sconstrindex)
            ar.sconstr(sconstrindex:end,:) = [];
        end
    end
end

if fiterrors == 1
    ar.chi2fit = ar.chi2 + ar.chi2err + ar.chi2err_logdataCorrection;
else
    ar.chi2fit = ar.chi2;
end

if(isfield(ar.model, 'data') && ~isempty(ar.res))
    ar.res_NaN = find(isnan(ar.res));
    if(sum(ar.res_NaN)>0)
        arDebugResidual;
        error('%i NaNs in residuals (check ar.res_NaN)', sum(isnan(ar.res)));
    else
        ar = rmfield(ar,'res_NaN');
    end

    if(sensi && ~isempty(ar.sres) && sum(sum(isnan(ar.sres(:,ar.qFit==1))))>0)
        for jm = 1:length(ar.model)
            if(isfield(ar.model(jm), 'data'))
                nd = length(ar.model(jm).data);
                for jd = 1:nd
                    if(ar.model(jm).data(jd).has_yExp)
                        ar.model(jm).data(jd).syExpSimu (:) = 0;
                        ar.model(jm).data(jd).systdExpSimu (:) = 0;
                        ar.model(jm).data(jd).sres(:) = 0;
                        ar.model(jm).data(jd).sreserr(:) = 0;
                    end
                end
            end
        end
        
        if ar.config.fiterrors == -1 % only exp errors
            warning('ar.config.fiterrors = -1 enforces usage of exp. errors. NaN in res or sres occur if no exp. Errors are in the data. Please check.')
        end
        arDebugResidual;
        error('NaN in derivative of residuals: %i', sum(sum(isnan(ar.sres(:,ar.qFit==1)))));
    end
end


if length(ar.res_type) ~= length(ar.res)
    error('arCollectRes.m: length(ar.res_type) ~= length(ar.res)')
end
