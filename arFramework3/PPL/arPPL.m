% arPPL(m, c, [ix], t, [takeY])
% 
% Major function for calculating prediction bands via profile likelihood
% for a dynamic variable or observable in a specific condition at specified
% time points. Additional time points can be appended by re-using arPPL.
% 
%   m   model index 
%   c	condition or data index, depending on takeY. 
%       Check ar.model(m).data(c).condition
%       and ar.model(m).condition(c).dLink for explicit condition
%   ix  Vector of states used for CI profile. See ar.model(m).xNames or
%       ar.model(m).data(c).yNames for the assignment
%   t   Vector of times where full Profile is computed (e.g. as comparison
%       to integrated profile), whichT flag is the vector index for
%       starting time of integration.  
% takeY [true] 
%       if takeY==true: an observation ar.model.data.y is used for profile
%       integration (data struct) 
%       if takeY==false: a ar.model.condition.x is used 
%
% Initialize calculation by calling options = PPL_options to use default 
% options or specify options by calling PPL_options with name-value pairs. 
%
% Helge Hass, 2014 (helge.hass@fdm.uni-freiburg.de)
% Tried to be documented by Clemens/Tim.
% 
% Example:
% >> PPL_options('Integrate',false,'doPPL',true) 
%       %Set option to calculate prediction profiles without integration
% >> arPPL(1,1,1,10,false)
%       %Run algorithm for time t = 10
% >> arPlotPPL(1,1,1,10,false)
%       %Plot prediction profile for t = 10
%
% See also: PPL_options, arPlotPPL, arPredictionProfile

function arPPL(m, c, ix, t, takeY) % model, condition, states of interest, 

    global ar
    ar.config.optim.Jacobian = 'on';
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
    if (~isfield(ar.ppl,'options')) || (isempty(ar.ppl.options))
        fprintf(['\n ERROR arPPL: Please initialize calculation by calling \n',...
            '       options = PPL_options \n',...
            'or call PPL_options with specified name-value pairs. \n'])
        return;
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
    
    % Integrated prediction bands are not directly supported! But to generate 
    % them anyway:
    % Alternative 1: Integrated validation bands with sufficiently small 
    %       measurement error.
    % Alternative 2: Prediction Profiles for several time points.
    if(ar.ppl.options.doPPL && ar.ppl.options.integrate)       
        %reset config fields 
        ar.config.fiterrors=fittederrors;
        ar.qFit(ar.qError==1)=fit_bkp;
        error(['Integration of confidence bands are not maintained and are not robust yet.',...
            ' Try approximating them by prediction bands with narrowing measurement error.'])
    end 
    
    % Loop over states
    for jx = 1:length(ix)
        arSimu(false, true, true); %is already done in PPL_init?
        ppl_general_struct.jx = ix(jx);
        % Try to set standard deviation of data point appropriately if no
        % value has been specified (xstd_auto = 0).
        if(ar.ppl.xstd_auto)  
            if(takeY)
                [~,it_first] = min(abs(ar.model(m).data(c).tExp-t(whichT))); 
                if(~isnan(ar.model(m).data(c).yExpStd(it_first,ix(jx))))
                    ar.ppl.options.xstd = ar.model(m).data(c).yExpStd(it_first,ix(jx));
                    % Takes experimental error of closest data point
                elseif(ar.config.fiterrors~= -1 && ~isnan(ar.model(m).data(c).ystdFineSimu(1,ix(jx))))
                    [~,it_first_fine] = min(abs(ar.model(m).data(c).tFine-t(whichT)));
                    ar.ppl.options.xstd = ar.model(m).data(c).ystdFineSimu(it_first_fine,ix(jx));
                    % Alternatively use simulated data error
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

