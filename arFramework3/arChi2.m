% Calculate chi^2 value and number of data points
%
% arChi2(sensi, pTrial, dynamics, doSimu)
%   sensi:          propagate sensitivities         [false]
%                   this argument is passed to arSimu
%   pTrial:         trial parameter of fitting
%   dynamics:       evaluate dynamics               [true]
%   doSimu          should arSimu be called         [true]
% 
% or
%
% ar = arChi2(ar, sensi, pTrial, dynamics)
%   ar:             d2d model/data structure
% 
% Possible calls:
% arChi2        then:
%               qglobalar = true
%               sensi = true
% 
% arChi2(ar)    here, the global ar is overwritten by the argument
%               qglobalar = false
%               sensi = true
% 
% arChi2(ar,sensi)
% arChi2(sensi)
%               like arChi2, but sensi can be set 
% 
% arChi2(ar,sensi,ptrial)
% arChi2(sensi,ptrial)
%               like arChi2(ar,sensi), but 
%               silent = true
%               ar.p is set to ptrial
%               this is the only possiblity to set 'silent' to true
% 
% arChi2(sensi,ptrial,dynamics)
% arChi2(ar,sensi,ptrial,dynamics)
%               dynamics is passed to arSimu, otherwise arSimu is called
%               without this argument, i.e. with its default 
%               
% arChi2(sensi,ptrial,dynamics,doSimu)
% arChi2(ar,sensi,ptrial,dynamics,doSimu)
%           doSimu  can be set to false, then the residuals are calculated
%           without updating the model trajectories (e.g. if ar.qFit or ar.
%           model.data.qFit) has been changed.

function varargout = arChi2(varargin)

global ar

nargs = nargin;
% The possiblity providing ar as an argument and to use of qglobalar==0 is
% obsolete because the gloal "ar" is overwritten anyway in arSimu 
% Implementation due to backwards compability:
if nargs>0 && isstruct(varargin{1})
    ar = varargin{1};
    varargin = varargin(2:end);
    nargs = nargs-1;    
    qglobalar = false;
else
    qglobalar = true;
end

if nargs>=1 && ~isempty(varargin{1})
    sensi = varargin{1};
else 
    sensi = true;
end

if nargs>=2 && ~isempty(varargin{2})
    pTrial = varargin{2};
	ar.p(ar.qFit==1) = pTrial;
    silent = true;
else
    silent = false;
end

if nargs>=3 && ~isempty(varargin{3})
    dynamics = varargin{3};
else
    dynamics = [];
end

if nargs>=4 && ~isempty(varargin{4})
    doSimu = varargin{4};
else
    doSimu = true;
end

if(~isfield(ar, 'fevals'))
    ar.fevals = 0; 
end
ar.fevals = ar.fevals + 1;

ar.ndata = 0;
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
if(~isfield(ar,'sres'))
    ar.sres = [];
end
if(~isfield(ar,'constr'))
    ar.constr = [];
end
if(~isfield(ar,'sconstr'))
    ar.sconstr = [];
end

nm = length(ar.model);
np = length(ar.p);

if(~isfield(ar.config,'useFitErrorCorrection'))
    ar.config.useFitErrorCorrection = true;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end
if ar.config.useFitErrorMatrix==1
    ar.ndata_err = 0;
end


% correction for error fitting
for jm = 1:nm
    if(isfield(ar.model(jm), 'data'))
        nd = length(ar.model(jm).data);
        for jd = 1:nd
            if(ar.model(jm).data(jd).has_yExp)
                ar.ndata = ar.ndata + sum(ar.model(jm).data(jd).ndata(ar.model(jm).data(jd).qFit==1));
                if(ar.config.useFitErrorMatrix == 1 && ar.config.fiterrors_matrix(jm,jd) == 1)
                    ar.ndata_err = ar.ndata_err + sum(ar.model(jm).data(jd).ndata(ar.model(jm).data(jd).qFit==1));
                end
            end
        end
    end
end
if( (ar.config.useFitErrorMatrix==0 && ar.ndata>0 && ar.config.fiterrors==1 && ar.config.useFitErrorCorrection) ...
        || (ar.config.useFitErrorMatrix==1 && ar.ndata>0 && sum(sum(ar.config.fiterrors_matrix==1))>0 && ar.config.useFitErrorCorrection) )
    if(ar.ndata-sum(ar.qError~=1 & ar.qFit==1) < sum(ar.qError~=1 & ar.qFit==1))
        ar.config.fiterrors_correction = 1;
        if(~ar.config.fiterrors_correction_warning)
            warning('ar.config.fiterrors_correction_warning : turning off bias correction, not enough data'); %#ok<WNTAG>
            ar.config.fiterrors_correction_warning = true;
        end
    else
        ar.config.fiterrors_correction = ar.ndata/(ar.ndata-sum(ar.qError~=1 & ar.qFit==1));
        ar.config.fiterrors_correction_warning = false;
    end
else
    ar.config.fiterrors_correction = 1;
end

atol = ar.config.atol;
rtol = ar.config.rtol;
maxsteps = ar.config.maxsteps;
qPositiveX = cell(1,length(ar.model));

if(~isfield(ar.config, 'nCVRestart'))
    nCVRestart = 1;
else
    nCVRestart = ar.config.nCVRestart;
end

for i = 1:nCVRestart
    try
        if doSimu
            if(qglobalar)  % since ar is overwritten anyway in arSimu, the possiblity to use of qglobalar obsolete
                arSimu(sensi, ~isfield(ar.model(jm), 'data'), dynamics);
            else
                ar = arSimu(ar, sensi, ~isfield(ar.model(jm), 'data'), dynamics);
            end
        end
        has_error = false;
        break
    catch error_id
        has_error = true;
        if(~silent)
            disp(error_id.message);
        end
        if nCVRestart > 1
            if strcmp(error_id.identifier,'MATLAB:UndefinedFunction')
                break
            else
                error_printed = 0;
                for m=1:length(ar.model)
                    for c=1:length(ar.model(m).condition)
                        if(ar.model(m).condition(c).status==-1)
                            % CV_TOO_MUCH_WORK
                            if(isempty(qPositiveX{m}))
                                qPositiveX{m} = ar.model(m).qPositiveX;
                                ar.model(m).qPositiveX(:) = 0;
                            else
                                ar.config.maxsteps = (1+.2*(i-1))*maxsteps;
                                if(~error_printed)
                                    arFprintf(1, 'Integration error, restarting %d / %d with 20%% increased maxsteps.\n',i-1,nCVRestart)
                                    error_printed = 1;
                                end
                            end
                        elseif(ar.model(m).condition(c).status<-1)
                            ar.config.atol = (1+.05*i)*atol;
                            ar.config.rtol = (1+.05*i)*rtol;
                            if(~error_printed)
                                arFprintf(1, 'Integration error, restarting %d / %d with 5%% increased precision.\n',i,nCVRestart)
                                error_printed = 1;
                            end
                        end
                    end
                end

            end
        end 
    end
end
ar.config.atol = atol;
ar.config.rtol = rtol;
ar.config.maxsteps = maxsteps;
for m=1:length(ar.model)
    if(~isempty(qPositiveX{m}))
        ar.model(m).qPositiveX = qPositiveX{m};
    end
end
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(isfield(ar.model(m).condition(c), 'xExpSimu'))
            if(sum((min(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),[],1) ./ range(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),1) < -ar.config.rtol) & (min(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),[],1) < -ar.config.atol)) > 0)
                arFprintf(1, 'Negative state in model %d condition %d detected that is defined as positive! Double-check model definition!\nPlot negative states by calling ar.model(%d).qPositiveX(:) = 0; with subsequent arPlot call.\n',m,c,m)
            end
        end
    end
end

useMSextension = false;

resindex = 1;
sresindex = 1;

% data
for jm = 1:nm
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
                    resindex = resindex+length(tmpres(:));
                    
                    if( (ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1) || ...
                            (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==1) )
                        ar.chi2err = ar.chi2err + sum(ar.model(jm).data(jd).chi2err(ar.model(jm).data(jd).qFit==1));
                        tmpreserr = ar.model(jm).data(jd).reserr(:,ar.model(jm).data(jd).qFit==1);
                        ar.res(resindex:(resindex+length(tmpreserr(:))-1)) = tmpreserr;
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
                        
                        if( (ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1) || ...
                                (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==1) )
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

% steady state conditions
qRelativeToInitialValue = true;
for jm = 1:nm
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
for jp=1:np
    if(ar.type(jp) == 0) % flat prior with hard bounds
    elseif(ar.type(jp) == 1) % normal prior
        tmpres = (ar.mean(jp)-ar.p(jp))./ar.std(jp);
        ar.res(resindex) = tmpres;
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
    elseif(ar.type(jp) == 3) % L1 prior
        tmpres = sqrt(abs((ar.mean(jp)-ar.p(jp))./ar.std(jp)));
        ar.res(resindex) = tmpres;
        resindex = resindex+1;
        if(ar.config.useSensis && sensi)
            tmpsres = zeros(size(ar.p));
            if ar.mean(jp) ~= ar.p(jp)
                tmpsres(jp) = sign(ar.p(jp)-ar.mean(jp)) ./ (2*ar.std(jp).*sqrt(abs((ar.mean(jp)-ar.p(jp))./ar.std(jp))));
            else
                tmpsres(jp) = 0;
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
%     for jm = 1:nm
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
%     for jm = 1:nm
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

% cut off too long arrays
if(isfield(ar.model, 'data'))
    if(~isempty(ar.res))
        if(length(ar.res)>=resindex)
            ar.res(resindex:end) = [];
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

if( (ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1) || ...
        (ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0) )
    ar.chi2fit = ar.chi2 + ar.chi2err;
else
    ar.chi2fit = ar.chi2;
end

% set Inf for errors
if(has_error)
    ar.res(:) = Inf;
    ar.chi2 = Inf;
    ar.chi2err = Inf;
    ar.chi2fit = Inf;
end

% calculate first order optimality criterion
if(sensi)
    res = [ar.res ar.constr];
    sres = [];
    if(~isempty(ar.sres))
        sres = ar.sres(:, ar.qFit==1);
    end
    if(~isempty(ar.sconstr))
        sres = [sres; ar.sconstr(:, ar.qFit==1)];
    end
    g = -2*res*sres; % gradient
    if(~isempty(g))
        onbound = [ar.p(ar.qFit==1)==ar.ub(ar.qFit==1); ar.p(ar.qFit==1)==ar.lb(ar.qFit==1)];
        exbounds = [g>0; g<0];
        qred = sum(onbound & exbounds,1)>0;
        ar.firstorderopt = norm(g(~qred));
%         fprintf('first order optimality criterion %f (%i)\n', ar.firstorderopt, -sum(qred));
    else
        qred = nan;
        ar.firstorderopt = nan;
    end
end

if(~silent)
    if(ar.ndata>0)
        if(ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==1)
            arFprintf(1, '-2*log(L) = %g, %i data points, ', ...
                2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit, ar.ndata);
        elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
            arFprintf(1, '-2*log(L) = %g, %i data points, ', ...
                2*ar.ndata_err*log(sqrt(2*pi)) + ar.chi2fit, ar.ndata);
        else
            arFprintf(1, 'global chi^2 = %g, %i data points, ', ar.chi2fit, ar.ndata);
        end
    end
    arFprintf(1, '%i free parameters', sum(ar.qFit==1));
    if(ar.chi2constr ~=0)
        arFprintf(1, ', %g violation of %i constraints', ar.chi2constr, ar.nconstr);
    end
    if(ar.chi2prior ~=0)
        arFprintf(1, ', %g violation of %i prior assumptions', ar.chi2prior, ar.nprior);
    end
    if(sensi)
        arFprintf(1, ', first order optimality criterion %g (%i)', ar.firstorderopt, -sum(qred));
    end
    arFprintf(1, '\n');
end

if(has_error)
    rethrow(error_id)
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
        for jm = 1:nm
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

if(nargout>0 && ~qglobalar)
    varargout{1} = ar;
else
    varargout = cell(0);
end

