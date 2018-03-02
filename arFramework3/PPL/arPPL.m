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

function arPPL(m, c, ix, t, takeY, options) % model, condition, states of interest, 

    global ar
    ar.config.optim.Jacobian = 'on';
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
    
    % optimizer settings (set only once)
    fittederrors=ar.config.fiterrors;
    ar.config.fiterrors=0;
    fit_bkp = ar.qFit(ar.qError==1);
    ar.qFit(ar.qError==1)=2;
    
    %fill temporary struct
    ppl_general_struct = struct;
    ppl_general_struct.m = m;
    ppl_general_struct.c = c;
    ppl_general_struct.t = t;
    ppl_general_struct.takeY = takeY;
    ppl_general_struct.x_vector = ix;
    
    confirm_options = PPL_options(options);
    fprintf(confirm_options)
    % Initialize PPL struct, set values to NaN that are newly computed
    t = PPL_init(ppl_general_struct);
    whichT = ar.ppl.options.whichT;
    dir = ar.ppl.options.dir;
    
    ppl_general_struct.pReset = ar.p;
    ppl_general_struct.chi2start = arGetMerit('chi2');

    tic;
    
    %Set vars distinguishing different setups
    if(takeY)
        data_cond = 'data';
        x_y = 'y';
    else
        data_cond = 'condition';
        x_y = 'x';
    end
    ppl_general_struct.data_cond = data_cond;
    ppl_general_struct.x_y = x_y;
    

    if(~ar.ppl.options.doPPL)
        ppl_vpl = 'vpl_';
    else
        ppl_vpl = 'ppl_';
    end
    ppl_general_struct.ppl_vpl = ppl_vpl;
    if(ar.ppl.options.doPPL && ar.ppl.options.integrate)       
        %reset config fields 
        ar.config.fiterrors=fittederrors;
        ar.qFit(ar.qError==1)=fit_bkp;
        error('Integration of confidence bands are not maintained and are not robust yet. Try approximating them by prediction bands with narrowing measurement error. \n')
    end 
    
    %loop over states
    for jx = 1:length(ix)
        arSimu(false, true, true);
        ppl_general_struct.jx = ix(jx);
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
            warning('Standard deviation of the auxiliary data point was set to %.2e, reset it manually if a different value is desired! \n',ar.ppl.options.xstd)
        end
        if(isnan(ar.ppl.options.xstd) || ar.ppl.options.xstd == 0)
            ar.ppl.options.xstd = 0.1;
            warning('The standard deviation for the prediction profile likelihood could not be set! Check your measurement errors or set the option "xstd" appropriately! \n')
        end
        
        %Run prediction profile likelihood for given time points
        arPredictionProfile(t, ppl_general_struct, true, 0, []);
        
        %Without integration, stop calculation here
        if(~ar.ppl.options.integrate)
            if(takeY)
                arLink(true,0.,true,ix(jx), c, m,NaN);
            end
            ar.p = ppl_general_struct.pReset;
            if(jx==length(ix))
                arWaitbar(-1);
                toc;
            end
            continue;
        end
        
        %Start integration function  
        if(ar.model(m).(data_cond)(c).ppl.([ppl_vpl 'lb_threshold_profile'])(whichT,ix(jx))==-Inf || ar.model(m).(data_cond)(c).ppl.([ppl_vpl 'ub_threshold_profile'])(whichT,ix(jx))==Inf)
            fprintf('Starting points for Profile at t=%d not defined, check PPL computation!', t(whichT));
            if(takeY)
                arLink(true,0.,true,ix(jx), c, m,NaN);
            end
            ar.p = ppl_general_struct.pReset;
            break;
        else    
            if(dir==0)
                for dir_tmp = [-1 1]                   
                    arIntegratePredBand(ppl_general_struct, dir_tmp);
                end
            else
                arIntegratePredBand(ppl_general_struct, dir);
            end
        end
        %End of integration
        arWaitbar(-1);
        toc;
        if(takeY)
            arLink(true,0.,true,ix(jx), c, m,NaN);
        end
        ar.p = ppl_general_struct.pReset;

    end
    %reset config fields 
    ar.config.fiterrors=fittederrors;
    ar.qFit(ar.qError==1)=fit_bkp;
    
    arCalcMerit(); 
    warning('Prediction profile plotting is enabled. To return to normal plots type ar.config.ploterrors = 0;');
    ar.config.ploterrors = -1;
    if(~takeY)
        ar.model(m).qPlotXs(ar.model(m).condition(c).dLink(1))=1;        
    else
        ar.model(m).qPlotYs(c) = 1;        
    end
end

