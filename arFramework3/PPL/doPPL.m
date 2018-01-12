% calculate prediction bands via profile likelihood
% 
% m: model to be used
%
% c: condition or data, respectively. Check ar.model(m).data(c).condition
% and ar.model(m).condition(c).dLink for explicit condition
%
% ix: Vector of states used for CI profile. See ar.model(m).xNames or
% ar.model(m).data(c).yNames for the assignment
%
% t: Vector of times where full Profile is computed (e.g. as comparison to
% integrated profile), whichT flag is the vector index for starting time of
% integration.
%
% takeY: =true ? observation is used for profile integration (data struct), otherwise condition is used (default is true)
%
% options: provide an option struct (type PPL_options for help)
%
% Helge Hass, 2014 (helge.hass@fdm.uni-freiburg.de)

function doPPL(m, c, ix, t, takeY, options) % model, condition, states of interest, 

    global ar
    
    if nargin < 6
      options = [];
    end
    if nargin < 2
        fprintf('Please specify the model, condition/data index, which states should be integrated \n, the time points and if an observation or internal state should be used. \n');
        return;  
    end

    if(~exist('takeY','var') || isempty(takeY))
        takeY = true;
        error('Not specified whether PBs on observation or internal state should be calculated. \n');
    end

    if(~exist('ix','var') || isempty(ix))
        fprintf('No specific state given, thus all are taken!\n');
        if(takeY)
            ix = 1:length(ar.model(m).data(c).yNames);
        else
            ix = 1:length(ar.model(m).xNames);
        end
    end
    confirm_options = PPL_options(options);
    fprintf(confirm_options)
    % Initialize PPL struct, set values to NaN that are newly computed
    t = PPL_init(m,c,t,ix,takeY);
    whichT = ar.ppl.options.whichT;
    dir = ar.ppl.options.dir;
    
    pReset = ar.p;
    chi2start = arGetMerit('chi2');

    tic;
    
    %Set vars distinguishing different setups
    if(takeY)
        data_cond = 'data';
    else
        data_cond = 'condition';
    end

    if(dir==1)
        high_low = '_high';
    elseif(dir==-1)
        high_low = '_low';
    end

    if(~ar.ppl.options.doPPL)
        ppl_vpl = '_vpl';
    else
        ppl_vpl = '';
    end
    
    if(ar.ppl.options.doPPL && ~ar.ppl.options.onlyProfile)
        warning('Integration of confidence bands are not maintained and are not robust yet. Try approximating them by prediction bands with narrowing measurement error. \n')
        return;
    end 
    
    %loop over states
    for jx = 1:length(ix)
        arSimu(false, true, true);
        %try to set standard dev of data point appropriately
        if(ar.ppl.xstd_auto)  
            if(takeY)
                [~,it_first] = min(abs(ar.model(m).data(c).tExp-t(whichT))); 
                 if(~isnan(ar.model(m).data(c).yExpStd(it_first,ix(jx))))
                     ar.ppl.options.xstd = ar.model(m).data(c).yExpStd(it_first,ix(jx));
                 elseif(ar.config.fiterrors~= -1 && takeY && ~isnan(ar.model(m).data(c).ystdFineSimu(1,ix(jx))))
                    ar.ppl.options.xstd = ar.model(m).data(c).ystdFineSimu(it_first,ix(jx));
                 end
            elseif(~takeY)
                ar.ppl.options.xstd = max(ar.model(m).condition(c).xFineSimu(:,ix(jx)))/10;
            end
            warning('Standard deviation of the auxiliary data point was set to %.2f, reset it manually if a different value is desired! \n',ar.ppl.options.xstd)
        end
        if(isnan(ar.ppl.options.xstd) || ar.ppl.options.xstd == 0)
            ar.ppl.options.xstd = 0.1;
            warning('The standard deviation for the prediction profile likelihood could not be set! Check your measurement errors or set the option "xstd" appropriately! \n')
        end
        
        %Run prediction profile likelihood for given time points
        xstart_ppl(m, c, ix(jx), t, ar.ppl.options.doPPL, ar.ppl.options.xstd, pReset, chi2start, whichT, takeY, true, 0, [], ar.ppl.options.onlyProfile);
        
        %Without integration, stop calculation here
        if(ar.ppl.options.onlyProfile)
            if(takeY)
                arLink(true,0.,true,ix(jx), c, m,NaN);
            end
            ar.p = pReset;
            if(jx==length(ix))
                arWaitbar(-1);
                toc;
            end
            continue;
        end
        
        %Start integration function  
        if(ar.model(m).(data_cond)(c).ppl.(['lb_fit' ppl_vpl])(whichT,ix(jx))==-Inf || ar.model(m).(data_cond)(c).ppl.(['ub_fit' ppl_vpl])(whichT,ix(jx))==Inf)
            fprintf('Starting points for Profile at t=%d not defined, check PPL computation!', t(whichT));
            if(takeY)
                arLink(true,0.,true,ix(jx), c, m,NaN);
            end
            ar.p = pReset;
            break;
        else    
            for dir_tmp = [-1 1]
                if(dir_tmp == -1)
                    high_low_tmp = '_low';
                    ub_lb = 'lb';
                else
                    high_low_tmp = '_high';
                    ub_lb = 'ub';
                end
                if(~isnan(ar.model(m).(data_cond)(c).ppl.(['kind' high_low_tmp ppl_vpl])(whichT, ix(jx))))
                    ar.model(m).(data_cond)(c).ppl.(['x' high_low_tmp ppl_vpl])(1, ix(jx))=ar.model(m).(data_cond)(c).ppl.([ub_lb '_fit' ppl_vpl])(whichT,ix(jx));
                    ar.model(m).(data_cond)(c).ppl.(['ps' high_low_tmp])(1, ix(jx),:)=squeeze(ar.model(m).(data_cond)(c).ppl.ps(whichT,ix(jx),ar.model(m).(data_cond)(c).ppl.(['kind' high_low_tmp ppl_vpl])(whichT, ix(jx)),:));
                end
                if(dir==0)
                    ppl_calc(m, c, ix(jx), ar.model(m).(data_cond)(c).ppl.(['x' high_low_tmp ppl_vpl])(1, ix(jx)), ar.model(m).(data_cond)(c).ppl.(['ps' high_low_tmp])(1, ix(jx),:), t(whichT), ar.ppl.options.doPPL, takeY, dir_tmp, ar.ppl.options.stepsize, ar.ppl.options.xstd, ar.ppl.options.ed_steps, pReset, chi2start, ar.ppl.options.backward, ar.ppl.options.fineInt);
                end
            end
            if(dir~=0)
                ppl_calc(m, c, ix(jx), ar.model(m).(data_cond)(c).ppl.(['x' high_low ppl_vpl])(1, ix(jx)), ar.model(m).(data_cond)(c).ppl.(['ps' high_low])(1, ix(jx),:), t(whichT), ar.ppl.options.doPPL, takeY, dir, ar.ppl.options.stepsize, ar.ppl.options.xstd, ar.ppl.options.ed_steps, pReset, chi2start, ar.ppl.options.backward, ar.ppl.options.fineInt);
            end
        end
        %End of integration
        arWaitbar(-1);
        toc;
        if(takeY)
            arLink(true,0.,true,ix(jx), c, m,NaN);
        end
        ar.p = pReset;

    end
    %reset config fields 
    if(isfield(ar.ppl,'fittederrors'))
        ar.config.fiterrors=ar.ppl.fittederrors;
        ar.qFit(ar.qError==1)=ar.ppl.fit_bkp;
    end
    arCalcMerit(); 
    warning('Prediction profile plotting is enabled. To return to normal plots type ar.config.ploterrors = 0;');
    ar.config.ploterrors = -1;
    if(~takeY)
        ar.model(m).qPlotXs(ar.model(m).condition(c).dLink(1))=1;        
    else
        ar.model(m).qPlotYs(c) = 1;        
    end
end

