% Collects all residuals in 
%   - ar.model.data.res         -> ar.res , ar.type=1
%   - ar.model.data.reserr      -> ar.res , ar.type=2
%   - prior                     -> ar.res , ar.type=3
%   - constr                    -> ar.constr
%   - random                    -> ar.res , ar.type=4
% and calculates chi2 values as well as the number of data points

function arCollectRes(sensi)

global ar 

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

np = length(ar.p);

useMSextension = false;

resindex = 1;
sresindex = 1;

% fit errors?
fiterrors = ( ar.config.fiterrors == 1  || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0 ) );

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
                    resindex = resindex+length(tmpres(:));
                    
                    if( fiterrors )
                        ar.chi2err = ar.chi2err + sum(ar.model(jm).data(jd).chi2err(ar.model(jm).data(jd).qFit==1));
                        tmpreserr = ar.model(jm).data(jd).reserr(:,ar.model(jm).data(jd).qFit==1);
                        ar.res(resindex:(resindex+length(tmpreserr(:))-1)) = tmpreserr;
                        ar.res_type(resindex:(resindex+length(tmpreserr(:))-1)) = 2;
                        resindex = resindex+length(tmpreserr(:));
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

        % Determine sensitivities w.r.t. log10(p) for the logtransformed ones
        trafo1              = ar.qLog10( pLink1 ) .* log(10) .* 10.^ar.p( pLink1 );
        trafo2              = ar.qLog10( pLink2 ) .* log(10) .* 10.^ar.p( pLink2 );
        for a = 1 : sum( pLink1 )
            sens1(:,:,a)    = sens1(:,:,a) .* trafo1(a);
        end
        for a = 1 : sum( pLink2 )
            sens2(:,:,a)    = sens2(:,:,a) .* trafo2(a);
        end

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
                            dxdp(:,ar.model(jm).condition(jc).qLog10 == 1) = bsxfun(@times, ... % log trafo
                                dxdp(:,ar.model(jm).condition(jc).qLog10 == 1), ...
                                ar.model(jm).condition(jc).pNum(ar.model(jm).condition(jc).qLog10 == 1) * log(10));
                        else
                            dxdp = squeeze(ar.model(jm).condition(jc).sxFineSimu(1,qss,:));
                            if(iscolumn(dxdp)) % transpose dxdp if squeeze returns column vector (sum(qss)==1)
                                dxdp = dxdp';
                            end
                        end
                        ddxdtdp = ar.model(jm).condition(jc).ddxdtdp(qss,:);
                        ddxdtdp(:,ar.model(jm).condition(jc).qLog10 == 1) = bsxfun(@times, ... % log trafo
                            ddxdtdp(:,ar.model(jm).condition(jc).qLog10 == 1), ...
                            ar.model(jm).condition(jc).pNum(ar.model(jm).condition(jc).qLog10 == 1) * log(10));
                        
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
                        ddxdtdp(:,ar.model(jm).condition(jc).qLog10 == 1) = bsxfun(@times, ... % log trafo
                            ddxdtdp(:,ar.model(jm).condition(jc).qLog10 == 1), ...
                            ar.model(jm).condition(jc).pNum(ar.model(jm).condition(jc).qLog10 == 1) * log(10));
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
    elseif(ar.type(jp) == 3 && ~isinf(ar.std(jp))) % L1 prior
        tmpres = sqrt(abs((ar.mean(jp)-ar.p(jp))./ar.std(jp)));
        ar.res(resindex) = tmpres;
        ar.res_type(resindex) = 3; 
        resindex = resindex+1;
        if(ar.config.useSensis && sensi)
            tmpsres = zeros(size(ar.p));
            if abs(ar.mean(jp) - ar.p(jp)) > 1e-10
                tmpsres(jp) = sign(ar.p(jp)-ar.mean(jp)) ./ (2*ar.std(jp).*sqrt(abs((ar.mean(jp)-ar.p(jp))./ar.std(jp))));
            elseif abs(2*ar.res(notpriors)*ar.sres(notpriors,jp)) < 1/ar.std(jp)
                ar.sres(:,jp) = 0;
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
        ar.res_type(resindex) = 4; 
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
        resindex = resindex + 1;
    end
    ar.chi2 = ar.chi2 + sum(ar.res_user.res.^2);
    ar.ndata = ar.ndata + length(ar.res_user.res);
end

for i=1:size(ar.res_user.sres,1)
    ar.sres(sresindex,:) = ar.res_user.sres(i,:);
    sresindex = sresindex + 1;
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
    ar.chi2fit = ar.chi2 + ar.chi2err;
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
        error('NaN in derivative of residuals: %i', sum(sum(isnan(ar.sres(:,ar.qFit==1)))));
    end
end


if length(ar.res_type) ~= length(ar.res)
    error('arCollectRes.m: length(ar.res_type) ~= length(ar.res)')
end
