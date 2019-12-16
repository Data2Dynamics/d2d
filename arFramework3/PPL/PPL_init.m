% t = PPL_init(ppl_general_struct)
% 
% Set ppl config and options which is stored in ar.ppl. This function is
% called by arPPL.m
% 
% ppl_general_struct    PPL result struct
% 
% Written by Helge, tried to be documented by Clemens.

    
function t = PPL_init(ppl_general_struct)
    
    global ar;
    
    %Set temporal variables
    m = ppl_general_struct.m;
    c = ppl_general_struct.c;
    t = ppl_general_struct.t;
    takeY = ppl_general_struct.takeY;
    ix = ppl_general_struct.x_vector;
    ar.ppl.qLog10=0;

    if(~isfield(ar.ppl,'xstd_auto') || isempty(ar.ppl.xstd_auto))
        ar.ppl.xstd_auto = 1;
    end
    
    if(~isfield(ar.ppl.options,'tEnd'))
        if(ar.ppl.options.backward)
            if(takeY)
                ar.ppl.options.tEnd = ar.model(m).data(c).tFine(1);
            else
                ar.ppl.options.tEnd = ar.model(m).condition(c).tFine(1); 
            end
        else
            if(takeY)
               ar.ppl.options.tEnd = ar.model(m).data(c).tFine(end);
            else
               ar.ppl.options.tEnd = ar.model(m).condition(c).tFine(end); 
            end
        end
    end

    if(~isfield(ar.ppl.options,'stepsize'))
        if(takeY)
            ar.ppl.options.stepsize = 1/(size(ar.model(m).data(c).tFine,1)/(ar.model(m).data(c).tLim(end)-ar.model(m).data(c).tLim(1)));
        else
            ar.ppl.options.stepsize = 1/(size(ar.model(m).condition(c).tFine,1)/(ar.model(m).condition(c).tFine(end)-ar.model(m).condition(c).tstart));
        end
    end
    ar.ppl.nsteps=abs(floor((ar.ppl.options.tEnd-t(ar.ppl.options.whichT)) / ar.ppl.options.stepsize));
    n = ar.ppl.options.n_steps_profile;
    nsteps = ar.ppl.nsteps;

    %Set correction strength
    if(~isfield(ar.ppl.options,'gamma'))
        ar.ppl.options.gamma = 1./ar.ppl.options.stepsize;
    end
    
    ar.ppl.dchi2 = arChi2inv(1-ar.ppl.options.alpha_level, 1);
    ar.ppl.dchi2;
    arCalcMerit();

    ar.ppl.chi2_threshold = arGetMerit('chi2')+ar.ppl.dchi2;

    
    %Set vars distinguishing different setups
    if(takeY)
        data_cond = 'data';
        x_y = 'y';
    else
        data_cond = 'condition';
        x_y = 'x';
    end
    if((~isfield(ar.model(m).(data_cond)(c),'ppl') || (isfield(ar.model(m).(data_cond)(c),'ppl') && isempty(ar.model(m).(data_cond)(c).ppl))))
    
    %Create fields for prediction profiles
        ar.model(m).(data_cond)(c).ppl.xtrial_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);    
        ar.model(m).(data_cond)(c).ppl.xfit_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);      
        ar.model(m).(data_cond)(c).ppl.vpl_likelihood_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);       
        ar.model(m).(data_cond)(c).ppl.ppl_likelihood_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);       
        ar.model(m).(data_cond)(c).ppl.ps_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1, length(ar.p)); 
        ar.model(m).(data_cond)(c).ppl.ppl_lb_threshold_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.ppl_ub_threshold_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.vpl_lb_threshold_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.vpl_ub_threshold_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));    
        ar.model(m).(data_cond)(c).ppl.ts_profile = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));           
        if(size(t,1)==1)
            ar.model(m).(data_cond)(c).ppl.ts_profile(:,ix) = repmat(t',1,length(ix));
        else
            ar.model(m).(data_cond)(c).ppl.ts_profile(:,ix) = repmat(t,1,length(ix));
        end
       
        %Create fields for bands in ar struct
        ar.model(m).(data_cond)(c).([x_y 'FineLB']) = nan(length(ar.model(m).(data_cond)(c).tFine),size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).([x_y 'FineUB']) = nan(length(ar.model(m).(data_cond)(c).tFine),size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        
        if(isfield(ar.model(m).(data_cond)(c),'zExpSimu'))
            ar.model(m).(data_cond)(c).('zFineLB') = nan(length(ar.model(m).(data_cond)(c).tFine),size(ar.model(m).(data_cond)(c).('zExpSimu'),2));
            ar.model(m).(data_cond)(c).('zFineUB') = nan(length(ar.model(m).(data_cond)(c).tFine),size(ar.model(m).(data_cond)(c).('zExpSimu'),2));
        end
        if(isfield(ar.model(m).(data_cond)(c),'uExpSimu'))
            ar.model(m).(data_cond)(c).('uFineLB') = nan(length(ar.model(m).(data_cond)(c).tFine),size(ar.model(m).(data_cond)(c).('uExpSimu'),2));
            ar.model(m).(data_cond)(c).('uFineUB') = nan(length(ar.model(m).(data_cond)(c).tFine),size(ar.model(m).(data_cond)(c).('uExpSimu'),2));
        end
        
        %Create fields for prediction bands
        ar.model(m).(data_cond)(c).ppl.band.tFine_band = nan(nsteps+1,size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));          
        ar.model(m).(data_cond)(c).ppl.band.xs_ppl_lowerBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));       
        ar.model(m).(data_cond)(c).ppl.band.xs_vpl_lowerBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));   
        ar.model(m).(data_cond)(c).ppl.band.ppl_likelihood_lowerBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));     
        ar.model(m).(data_cond)(c).ppl.band.vpl_likelihood_lowerBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));     
        ar.model(m).(data_cond)(c).ppl.band.xs_ppl_upperBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));      
        ar.model(m).(data_cond)(c).ppl.band.xs_vpl_upperBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));  
        ar.model(m).(data_cond)(c).ppl.band.ppl_likelihood_upperBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));    
        ar.model(m).(data_cond)(c).ppl.band.vpl_likelihood_upperBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));          
        ar.model(m).(data_cond)(c).ppl.band.ps_lowerBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), length(ar.p));    
        ar.model(m).(data_cond)(c).ppl.band.ps_uperBand = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), length(ar.p)); 
    end
    arSimu(false, true, true);
    if(ar.ppl.options.integrate)
        for jx=1:length(ix)
            %Delete fields for new integration
            ar.model(m).(data_cond)(c).ppl.band.tFine_band(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.xs_ppl_lowerBand(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.xs_vpl_lowerBand(:,ix(jx)) = nan;
            if(size(ar.model(m).(data_cond)(c).([x_y 'FineSimu']),2) < ix(jx))
                error('The state you want to calculate does not exist!')
            end
            ar.model(m).(data_cond)(c).ppl.band.xs_ppl_upperBand(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.xs_vpl_upperBand(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.ppl_likelihood_lowerBand(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.vpl_likelihood_lowerBand(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.ppl_likelihood_upperBand(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.vpl_likelihood_upperBand(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.band.ps_lowerBand(:,ix(jx),:) = nan;
            ar.model(m).(data_cond)(c).ppl.band.ps_uperBand(:,ix(jx),:) = nan;    
            
            ar.model(m).(data_cond)(c).([x_y 'FineLB'])(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).([x_y 'FineUB'])(:,ix(jx)) = nan;      
        end
    end
    for jt=1:length(t)
        for jx=1:length(ix)
            cur_t = find(ar.model(m).(data_cond)(c).ppl.ts_profile(:,ix(jx))==t(jt));
            %Append fields for profiles if new time points are added
            
            if(isempty(cur_t))
%                 if(size(ar.model(m).(data_cond)(c).ppl.ts_profile,1)>=jt)
%                     ar.model(m).(data_cond)(c).ppl.ts_profile(jt,ix(jx)) = t(jt);
%                 else                
                    ar.model(m).(data_cond)(c).ppl.ts_profile = [ar.model(m).(data_cond)(c).ppl.ts_profile; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2))];
                    ar.model(m).(data_cond)(c).ppl.ts_profile(end,ix(jx)) = t(jt);
                    ar.model(m).(data_cond)(c).ppl.xtrial_profile = [ar.model(m).(data_cond)(c).ppl.xtrial_profile; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1)];
                    ar.model(m).(data_cond)(c).ppl.xfit_profile = [ar.model(m).(data_cond)(c).ppl.xfit_profile; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1)];
                    ar.model(m).(data_cond)(c).ppl.vpl_likelihood_profile = [ar.model(m).(data_cond)(c).ppl.vpl; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1)];
                    ar.model(m).(data_cond)(c).ppl.ppl_likelihood_profile = [ar.model(m).(data_cond)(c).ppl.ppl; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1)];
                    ar.model(m).(data_cond)(c).ppl.ps_profile = [ar.model(m).(data_cond)(c).ppl.ps; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1, length(ar.p))];
                    ar.model(m).(data_cond)(c).ppl.ppl_lb_threshold_profile = [ar.model(m).(data_cond)(c).ppl.ppl_lb_threshold_profile; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2))];
                    ar.model(m).(data_cond)(c).ppl.ppl_ub_threshold_profile = [ar.model(m).(data_cond)(c).ppl.ppl_ub_threshold_profile; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2))];
                    ar.model(m).(data_cond)(c).ppl.vpl_lb_threshold_profile = [ar.model(m).(data_cond)(c).ppl.vpl_lb_threshold_profile; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2))];
                    ar.model(m).(data_cond)(c).ppl.vpl_ub_threshold_profile = [ar.model(m).(data_cond)(c).ppl.vpl_ub_threshold_profile; nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2))];
%                 end
            end
        end
    end

%     whichT = find(ar.model(m).(data_cond)(c).ppl.ts_profile(:,ix(1)) == t(whichT));
%     t = ar.model(m).(data_cond)(c).ppl.ts_profile(:,ix(1));
end