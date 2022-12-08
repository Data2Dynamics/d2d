%function [] = arCollectResCov(sensi)
% This function creates a "covariance version" ar.resCov of the residuals 
% array ar.res. This is done by replacing the ordinary residuals of type 
% ar.res_type==1 and ar.res_type==2 (see arCollectRes) by the transformed 
% ones (see arCalcResCov). This function is called subsequent to
% arCalcResCov in certain merit functions in arFit (optimizer 20 and 21) as
% well as in arGetMerit (if ar.chi2cov is used).
%
% The following fields in the ar struct are filled by this function:
%      ar.resCov
%      ar.sresCov
%      ar.chi2cov
%      ar.chi2err_cov
%      ar.chi2cov_fit
% They are used in the "covariance optimizers", i.e. ar.config.optimizer=20
% or 21.
%
% see also arCalcResCov, arFit, arGetMerit, arCollectRes

function [] = arCollectResCov(sensi)
if ~exist('sensi','var') || isempty(sensi)
    sensi = 1;
end

global ar

if(~isfield(ar,'resCov'))
    ar.resCov = [];
end
if(~isfield(ar,'sresCov'))
    ar.sresCov = [];
end
ar.chi2cov = 0;
ar.chi2err_cov = 0;
ar.chi2cov_fit = 0;

np = length(ar.p);

useMSextension = false;

resindex = 1;
sresindex = 1;

% fit errors?
fiterrors = ( ar.config.fiterrors == 1  || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0 ) );

ar.sres_size = []; % for NEB merger

for jm = 1:length(ar.model)
    if(isfield(ar.model(jm), 'data'))
        for jd = 1:length(ar.model(jm).data)
            if(ar.model(jm).data(jd).has_yExp)
                ar.chi2cov = ar.chi2cov + sum(ar.model(jm).data(jd).chi2cov(ar.model(jm).data(jd).qFit==1));
                
                if(useMSextension && isfield(ar, 'ms_count_snips') && ~isempty(ar.model(jm).data(jd).ms_index))
                    
                else
                    % collect residuals for fitting
                    tmpres = ar.model(jm).data(jd).resCov(:,ar.model(jm).data(jd).qFit==1);
                    ar.resCov(resindex:(resindex+length(tmpres(:))-1)) = tmpres;
                    %{
                    ar.res_type(resindex:(resindex+length(tmpres(:))-1)) = 1;
                    
                    if ( debugres )
                        for jr = 1 : numel( tmpres )
                            ar.resinfo(resindex + jr - 1).type = 'data';
                            ar.resinfo(resindex + jr - 1).m = jm;
                            ar.resinfo(resindex + jr - 1).d = jd;
                        end
                    end
                    %}
                    resindex = resindex+length(tmpres(:));
                    
                    if( fiterrors )                        
                        ar.chi2err_cov = ar.chi2err_cov + sum(ar.model(jm).data(jd).chi2err_cov(ar.model(jm).data(jd).qFit==1));
                        tmpreserr = ar.model(jm).data(jd).reserrCov(:,ar.model(jm).data(jd).qFit==1);
                        ar.resCov(resindex:(resindex+length(tmpreserr(:))-1)) = tmpreserr;
                        %{
                        ar.res_type(resindex:(resindex+length(tmpreserr(:))-1)) = 2;
                        if ( debugres )
                            for jr = 1 : numel( tmpreserr )
                                ar.resinfo(resindex + jr - 1).type = 'error model';
                                ar.resinfo(resindex + jr - 1).m = jm;
                                ar.resinfo(resindex + jr - 1).d = jd;
                            end
                        end
                        %}
                        resindex = resindex+length(tmpreserr(:));
                        
                        %{
                        % add correction term for fitting on log-axis due
                        % to normalization term of log-normal distribution
                        islogfitted = logical(ar.model(jm).data(jd).logfitting);
                        if any(islogfitted)
                            data = ar.model(jm).data(jd).yExp(islogfitted);
                            data = data(~isnan(data));
                            ar.chi2err_logdataCorrection = ar.chi2err_logdataCorrection+ sum(log(log(10)^2 * (10.^data).^2));
                        end
                        %}
                    end
                    
                    % add correction term for logarithmic fitting
                    if any(ar.model(jm).data(jd).logfitting)
                    end
                    
                    % collect sensitivities for fitting
                    if(ar.config.useSensis && sensi)
                        tmptmpsres = ar.model(jm).data(jd).sresCov(:,ar.model(jm).data(jd).qFit==1,:);
                        tmpsres = zeros(length(tmpres(:)), np);
                        tmpsres(:,ar.model(jm).data(jd).pLink) = reshape(tmptmpsres, ...
                            length(tmpres(:)), sum(ar.model(jm).data(jd).pLink));
                        ar.sresCov(sresindex:(sresindex+length(tmpres(:))-1),:) = tmpsres;
                        sresindex = sresindex+length(tmpres(:));
                        
                        if ( fiterrors )
                            tmpsreserr = zeros(length(tmpreserr(:)), np);
                            tmpsreserr(:,ar.model(jm).data(jd).pLink) = reshape(ar.model(jm).data(jd).sreserrCov(:,ar.model(jm).data(jd).qFit==1,:), ...
                                length(tmpreserr(:)), sum(ar.model(jm).data(jd).pLink));
                            ar.sresCov(sresindex:(sresindex+length(tmpres(:))-1),:) = tmpsreserr;
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
                    grad = ar.sresCov([ar.sres_size(jm-1):ar.sres_size(jm)-1],ar.merger.neb.ps_index(jm,:));

                    proj = nan(size(grad));
                    for ig = 1: size(grad,1)
                        proj(ig,:) = (sum(grad(ig,:) .* r)./norm(r) ).* r ./norm(r);
                    end

                    ar.sresCov([ar.sres_size(jm-1):ar.sres_size(jm)-1] ,ar.merger.neb.ps_index(jm,:)) ...
                        = grad-proj;
                end
                
            end
        end
    end
end

ar.chi2cov_fit = ar.chi2cov + ar.chi2err_cov;

end

