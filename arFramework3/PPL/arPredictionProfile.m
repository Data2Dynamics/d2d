% [xFit, ps] = arPredictionProfile(t, general_struct, save, dir, xFit)
%
% This function computes the validation/prediction profile for various input 
% time points. It is called in arPPL and arIntegratePredBand.
%
%   t:              Vector of times where full Profile is computed
%   general_struct: General quantities from the higher level function, like
%                   model/data index, state index, ppl or vpl, etc.
%   save:           Saves profiles in ar. If save is turned on, both directions
%                   are performed
%   dir:            Direction in which profile is calculated  
%   xFit, ps:       Has something to do with the integration functionality
%
% The available documentation focused on the non-integration aspects.
% Calculation of prediction profiles is performed by mapping the validation
% profiles to prediction profiles by use of an existing relation.
%
% Written by Helge, tried to be documented by Clemens/Tim
%
% See also: arPPL, arIntegratePredBand

function [xFit, ps] = arPredictionProfile(t, general_struct, save, dir, xFit)

global ar;
arWaitbar(0);

if(~exist('dir','var'))
    dir = 0;
end
if(~exist('xFit','var'))
    xFit = 1.;
    save = true;
end

% In arPPL, xFit = []  at this point

%fill temporary variables
pReset = ar.p;
data_cond = general_struct.data_cond;
x_y = general_struct.x_y;
ppl_vpl = general_struct.ppl_vpl;
m=general_struct.m;
c=general_struct.c; 
jx=general_struct.jx; 
takeY=general_struct.takeY;
chi2start = general_struct.chi2start;

integrate = ar.ppl.options.integrate;
xstd = ar.ppl.options.xstd;
doPPL=ar.ppl.options.doPPL;

for ts = 1:length(t)
    t_tmp = t(ts);
    if(isnan(t_tmp))
        continue;
    end
    %if(save)
    ts_tmp = find(ar.model(m).(data_cond)(c).ppl.ts_profile(:,jx) == t_tmp);
    % index of the entries corresponding to time t_tmp and states jx
    
    %Skip calculation if it already exists
    if(save && ((doPPL && ~isnan(ar.model(m).(data_cond)(c).ppl.ppl_ub_threshold_profile(ts_tmp,jx))) ...
            || (~doPPL && ~isnan(ar.model(m).(data_cond)(c).ppl.vpl_ub_threshold_profile(ts_tmp,jx)))) ...
            )
%            && ~(ts == ar.ppl.options.whichT && ar.model(m).(data_cond)(c).ppl.([ppl_vpl 'ub_threshold_profile'])(ts_tmp, jx) ~= ...
%            ar.model(m).(data_cond)(c).ppl.band.(['xs_' ppl_vpl 'upperBand'])(1, jx))) 
        % last line to recalculate profile if starting time for integration is changed
        % But save = false for integration anyway, does the last condition
        % matter?
        
       fprintf(['\n The prediction profile you want to compute for t=%d and state %d already exists. \n',...
           'If you want to overwrite them, delete the value in \n',...
           'ar.model.data/condition.ppl.ppl/vpl_ub_threshold_profile \n'],t_tmp,jx)

       continue; 
    end            
   
    % Add the current optimal prediction as a data point:
    [~,it_first] = min(abs(ar.model(m).(data_cond)(c).tFine-t_tmp));
    arLink(true, t_tmp, takeY, jx, c, m, ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(it_first,jx), xstd);
    % Calling arLink in this fashion adds the corresponding data point to the ar-struct
    % Note that this is only approximately the best prediction due to the
    % finite resolution of tFine.
    arCalcMerit(0,[],1); 
    % Recalculate objective function and simulate trajectory values at 
    % data times (especially for new data point!)
    [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
    if(takeY && length(find(ar.model(m).(data_cond)(c).tExp==t_tmp))>1)
        it = it+1;               
    end
    xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);
    % xSim is now the exact optimal prediction (without adding a data
    % point!, adding a data point affects the position of the optimum)
    arLink(true, t_tmp, takeY, jx, c, m, xSim, xstd);
    % Replace the approximate optimal prediction with the exact version
    if(ar.ppl.qLog10)
        xSim = log10(xSim);
    end
    fprintf('Calculating PPL for t=%d and state number x=%i \n',t(ts),jx);
    
    % When calling this function in arPPL, save = 1. Thus, up and down
    % directions are both performed no matter what was specified in dir. 
    % However, save = 0 if arPredictionProfile is called from 
    % arIntegratePredBand such that dir might have an impact in this
    % scenario.
    
    % go up
    if(save || dir==1)
        % This is where the profile is actually sampled:
        [xtrial_up, xfit_up, ppl_up, vpl_up, ps_up] = ...
            ppl(general_struct, it, t_tmp, xSim, 1);    
               
        ar.p = pReset;  
        arCalcMerit();       
    else
        % Why is this alternative necessary?
        xtrial_up = xSim*1.01;
        xfit_up = xSim*1.01;
        ppl_up = chi2start+0.1;
        vpl_up = chi2start+0.1;
        ps_up = ar.p;
    end
    % go down 
    if(save || dir==-1)
        [xtrial_down, xfit_down, ppl_down, vpl_down, ps_down] = ...
            ppl(general_struct, it, t_tmp, xSim, -1);
        % Reset parameters
        ar.p = pReset;
    else
        xtrial_down = xSim*0.99;
        xfit_down = xSim*0.99;
        ppl_down = chi2start+0.1;
        vpl_down = chi2start+0.1;
        ps_down = ar.p;
        ar.p = pReset;
    end
    
    % Reset data point
    if(takeY)
        arLink(true,0.,true,jx, c, m,NaN);
    end
    % Reset objective function
    arCalcMerit();
    % This should presumably now equal chi2start again
    
    % Join results of both directions:
    ps_tmp = [flipud(ps_down); pReset; ps_up];
    xfit_tmp = [fliplr(xfit_down) xSim xfit_up];
        % data points for the validation profile
    xtrial_tmp = [fliplr(xtrial_down) xSim xtrial_up];
        % predictions for the prediction profile
    ppl_tmp = [fliplr(ppl_down) chi2start ppl_up];
    vpl_tmp = [fliplr(vpl_down) chi2start vpl_up];
    if(doPPL)
        fitted_tmp = [fliplr(xfit_down) xSim xfit_up];
        merit_tmp = [fliplr(ppl_down) chi2start ppl_up];
    else
        fitted_tmp = [fliplr(xtrial_down) xSim xtrial_up];
        merit_tmp = [fliplr(vpl_down) chi2start vpl_up];
    end    
    q_chi2good = merit_tmp <= chi2start+ar.ppl.dchi2;    
    q_nonnan = ~isnan(merit_tmp); 
    % Query:
    % 1. chi2start is the objetive function value without adding a data
    %       point
    % 2. the objective function changes after adding a data point
    % 3. the profile is calculated with the additional data points
    %       objective
    % 4. Thus: Is it considered, that the profile values are not
    %       necessarily comparable to chi2start? (maybe by use of the
    %       custom merit function?)    


    % calculate CI point-wise fit
    lb_tmp = min(fitted_tmp(q_chi2good));
    ub_tmp = max(fitted_tmp(q_chi2good));
    find_tmp = find(merit_tmp > chi2start+ar.ppl.dchi2);
    % Check downward direction:
    if(length(vpl_down)==1 || sum(find_tmp<length(vpl_down)+1)==0 || (~save && dir==1))
        % If downward direction has not been attempted or no downward merit 
        % function value exceeds threshold
        lb_tmp = -Inf;             
        kind_low_tmp = -Inf;                
        fprintf('No -95 lower bound for PPL for t=%d\n',t(ts));
        % ...But alpha level may actually be changed. 
    else
        if(lb_tmp==min(fitted_tmp(q_nonnan)))
            warning(['Multiple likelihood values are assigned to the same model fit. ' ...
                     'check model uncertainty and fits, or set more strict integrator tolerances!'])
            % I don't get it
        end
        kind_low_tmp = find(q_chi2good==1,1,'first');
            % Index of lower bound
        if(length(kind_low_tmp)>1)
            kind_low_tmp = kind_low_tmp(1);
        end
        %Interpolate between just under threshold to just above the
        % threshold to get a good confidence bound:
        lb_tmp = interp1(merit_tmp([kind_low_tmp kind_low_tmp-1]), ...
        fitted_tmp([kind_low_tmp kind_low_tmp-1]), chi2start+ar.ppl.dchi2);
        % This should work because there exists a value exceeding the
        % threshold by the first if statement
    end
    % Check upward direction:
    if(length(vpl_up)==1 || sum(find_tmp>length(vpl_down)+1)==0 || (~save && dir==-1))
        ub_tmp = Inf;
        kind_high_tmp = Inf;        
        fprintf('No -95 lower bound for PPL for t=%d\n',t(ts));
    else
        if(ub_tmp==max(fitted_tmp(q_nonnan)))
            warning(['Multiple likelihood values are assigned to the same model fit. ' ...
                     'check model uncertainty and fits, or set more strict integrator tolerances!'])
        end
        kind_high_tmp = find(q_chi2good==1,1,'last');   
        if(length(kind_high_tmp)>1)
            kind_high_tmp = kind_high_tmp(end);
        end
        ub_tmp = interp1(merit_tmp([kind_high_tmp kind_high_tmp+1]), ...
        fitted_tmp([kind_high_tmp kind_high_tmp+1]), chi2start+ar.ppl.dchi2);
    end      
    
    if integrate && ~save
        if((dir==1 && isinf(ub_tmp)) || (dir==-1 && isinf(lb_tmp)))
            fprintf('ERROR: no bound found at t=%d \n', t_tmp);
            if(takeY)
                arLink(true,0.,true,jx, c, m,NaN);
            end
            break;
        end
        if(dir==1)        
            xFit = ub_tmp;
            ps = ps_tmp(kind_high_tmp,:);
        elseif(dir==-1)        
            xFit = lb_tmp;
            ps = ps_tmp(kind_low_tmp,:);                    
        end
    end
    % Save results in ar-struct:
    if(save)
        ar.model(m).(data_cond)(c).ppl.xtrial_profile(ts_tmp, jx,:) = xtrial_tmp;
        ar.model(m).(data_cond)(c).ppl.xfit_profile(ts_tmp, jx,:) = xfit_tmp;
        ar.model(m).(data_cond)(c).ppl.ppl_likelihood_profile(ts_tmp, jx,:) = ppl_tmp;
        ar.model(m).(data_cond)(c).ppl.vpl_likelihood_profile(ts_tmp, jx,:) = vpl_tmp;
        ar.model(m).(data_cond)(c).ppl.ps_profile(ts_tmp, jx,:,:) = ps_tmp;
        ar.model(m).(data_cond)(c).ppl.([ppl_vpl 'lb_threshold_profile'])(ts_tmp, jx) = lb_tmp;
        ar.model(m).(data_cond)(c).ppl.([ppl_vpl 'ub_threshold_profile'])(ts_tmp, jx) = ub_tmp;       
        if(ts == ar.ppl.options.whichT && integrate)
            ar.model(m).(data_cond)(c).ppl.band.tFine_band(1,jx) = t(ts);
            for dir = [-1 1]
                if(dir == 1)
                    high_low = 'upperBand';
                    xFit = ub_tmp;
                    ps = ps_tmp(kind_high_tmp,:);
                else
                    high_low = 'lowerBand';
                    xFit = lb_tmp;
                    ps = ps_tmp(kind_low_tmp,:);     
                end           
                ar.model(m).(data_cond)(c).ppl.band.(['xs_' ppl_vpl high_low])(1, jx) = xFit;  
                ar.model(m).(data_cond)(c).ppl.band.(['ps_' high_low])(1, jx,:) = ps; 
            end
        end
    end
end
%Extend LB/UB in ar struct to have values over the fine time grid:
    if(save && ~integrate && dir==0 && length(t)>1)        
        struct_vec = {'FineUB','FineLB'};
        low_high_vec = {'ub','lb'};
        for ilh=1:2
            struct_string = struct_vec{ilh};
            ppl_string = [ppl_vpl low_high_vec{ilh} '_threshold_profile'];
            if(takeY)
                struct_string = ['y' struct_string];
            else
                struct_string = ['x' struct_string];
            end            
            %Interpolate boundary values between time points to obtain
            %boundary values on fine time grid:
            ar.model(m).(data_cond)(c).(struct_string)(1:length(ar.model(m).(data_cond)(c).tFine),jx) = ...
                interp1(ar.model(m).(data_cond)(c).ppl.ts_profile(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
                ar.model(m).(data_cond)(c).ppl.(ppl_string)(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
                ar.model(m).(data_cond)(c).tFine,'pchip',NaN);   
            
        end
    end
end


function [xtrial, xfit, ppl, vpl, ps] = ppl(general_struct, it, t, xSim, direction)

% Calculates the validation/prediction profile.
% Prediction Profiles are implictly calculated by use of validation
% profiles.

global ar

m=general_struct.m;
c=general_struct.c; 
ix=general_struct.jx; 
takeY=general_struct.takeY;
chi2start = general_struct.chi2start;
qLog10 = ar.ppl.qLog10;
xstd = ar.ppl.options.xstd;
doPPL=ar.ppl.options.doPPL;
n = ar.ppl.options.n_steps_profile;
tcount = 1;

xtrial = nan(1,n);
xfit = nan(1,n);
ppl = nan(1,n);
vpl = nan(1,n);
ps = nan(n,length(ar.p));

% These are the step sizes in which data points are changed
% This is a fixed stepsize! => Inefficient
% Better: See function VPL
dx = sqrt(ar.ppl.dchi2*ar.ppl.options.rel_increase) * xstd;

if(takeY)
    xLabel = myNameTrafo(ar.model(m).data(c).y{ix});
else
    xLabel = myNameTrafo(ar.model(m).x{ix});    
end

% Start from the optimum:
xExp = xSim;

for j = 1:n
    if(toc>tcount)
        if(direction>0)
            arWaitbar((j), n, sprintf('PPL (up) for %s at t=%g %i/%i', xLabel, t, j, n));
        else
            arWaitbar((j), n, sprintf('PPL (down) for %s at t=%g %i/%i', xLabel, t, j, n));
        end
        tcount = tcount + 0.5; % update every half second
    end
    
    % Make the step:
    xExp = xExp + direction*dx;
    arLink(true,t,takeY,ix, c,m,xExp,xstd);
    
    try
        % Fit to new optimum
        arPPLFit;
    catch exception
        fprintf('ERROR in PPL integration (%s)\n', exception.message);
        if(takeY)
            arLink(true,0.,true,ix, c, m,NaN);
        end
        break;
    end
    
    xtrial(j) = xExp; % data point (unfitted)    
    arCalcMerit(0, ar.p(ar.qFit==1),1) %Simulate the predictions
    if(takeY)
        xSim = ar.model(m).data(c).yExpSimu(it,ix); 
    else
        xSim = ar.model(m).condition(c).xExpSimu(it,ix); 
    end
    xfit(j) = xSim; % corresponding prediction (fitted)
    
    % Get PPL values by their relation to VPL
    if(takeY)
        ppl(j) = arGetMerit('chi2') - (xtrial(j)-xfit(j)).^2/xstd.^2;
        vpl(j) = arGetMerit('chi2');
        
    else
        ppl(j) = arGetMerit('chi2');
        vpl(j) = arGetMerit('chi2') + (xtrial(j)-xfit(j)).^2/xstd.^2;
    end
    ps(j,:) = ar.p;

    % Stop if the confidence threshold is crossed:
    if((doPPL && ppl(j) > chi2start+ar.ppl.dchi2*1.2) || (~doPPL && vpl(j) > chi2start+ar.ppl.dchi2*1.2))
        break
    end
end

    function arPPLFit
        % Fit with custom merit function (why is this necessary?)
        
        [pFit, ~, ~, exitflag] = ...
            lsqnonlin(@ppl_merit_fkt, ar.p(ar.qFit==1), ar.lb(ar.qFit==1), ar.ub(ar.qFit==1), ar.config.optim);
        ar.p(ar.qFit==1) = pFit;
        if(exitflag < 1)
            outputstr = '';
            switch exitflag
                case 1
                    outputstr = 'LSQNONLIN converged to a solution.';
                case 2
                    outputstr = 'Change in X too small.';
                case 3
                    outputstr = 'Change in RESNORM too small.';
                case 4
                    outputstr = 'Computed search direction too small.';
                case 0
                    outputstr = 'Too many function evaluations or iterations.';
                case -1
                    outputstr = 'Stopped by output/plot function.';
                case -2
                    outputstr = 'Bounds are inconsistent.';
                case -3
                    outputstr = 'Regularization parameter too large (Levenberg-Marquardt).';
                case -4
                    outputstr = 'Line search failed.';
            end
            %fprintf('lsqnonlin finished after %i interations: %s\n', ar.fit.output.iterations, outputstr);
        end
    end


    function [res, sres] = ppl_merit_fkt(pTrial)
        
        arCalcMerit(ar.config.useSensis, pTrial, 1)
        
        res = [ar.res ar.constr];
        %res = ar.res;
        if(nargout>1 && ar.config.useSensis)
            sres = ar.sres(:, ar.qFit==1);
            if(~isempty(ar.sconstr))
                sres = [ar.sres(:, ar.qFit==1); ar.sconstr(:, ar.qFit==1)];
            end
        end
        
        
        if(~takeY)% && ~doPPL)
            xSim = ar.model(m).condition(c).xExpSimu(it,ix);
            if(qLog10)
                    xSim = log10(xSim);
            end	
            res(end+1) = (xExp-xSim)/xstd;
        
            if(nargout>1 && ar.config.useSensis)
                sx_tmp = arTrafoParameters(ar.model(m).condition(c).sxExpSimu,m,c,false);
                sxSim = zeros(1,length(ar.p));           
                sxSim(ar.model(m).condition(c).pLink) = ...
                        squeeze(sx_tmp(it,ix,:))';
%                 for j10=find(ar.qLog10==1)
%                     sxSim(j10) = sxSim(j10) * 10.^ar.p(j10) * log(10);
%                 end
                if(qLog10)
                    sxSim = sxSim / 10^xSim / log(10);
                end

                sres(end+1,:) = -sxSim(ar.qFit==1) / xstd;
            end    
        end

    end
end

function str = myNameTrafo(str)
str = strrep(str, '_', '\_');
end