function [t, whichT] = PPL_init(m,c,t,ix,gammas, onlyProfile, whichT,takeY)
    global ar;
    n = ar.ppl.n;
    nsteps = ar.ppl.nsteps;
    %Set vars distinguishing different setups
    if(takeY)
        data_cond = 'data';
        x_y = 'y';
    else
        data_cond = 'condition';
        x_y = 'x';
    end
    if((~isfield(ar.model(m).(data_cond)(c),'ppl') || (isfield(ar.model(m).(data_cond)(c),'ppl') && isempty(ar.model(m).(data_cond)(c).ppl)) ...
        ))%|| (isfield(ar.model(m).(data_cond)(c),'ppl') && ~isequal(ar.model(m).(data_cond)(c).ppl.tstart(:,ix(1))',t))))
        ar.model(m).(data_cond)(c).ppl.ix = ix(:);   
        ar.model(m).(data_cond)(c).ppl.xtrial = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);
        ar.model(m).(data_cond)(c).ppl.xfit = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);
        ar.model(m).(data_cond)(c).ppl.vpl = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);
        ar.model(m).(data_cond)(c).ppl.ppl = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1);
        ar.model(m).(data_cond)(c).ppl.ps = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1, length(ar.p));
        ar.model(m).(data_cond)(c).ppl.lb_fit = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.ub_fit = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.lb_fit_vpl = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.ub_fit_vpl = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.kind_high = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.kind_low = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.kind_high_vpl = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.kind_low_vpl = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.tstart = nan(length(t), size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        if(size(t,1)==1)
            ar.model(m).(data_cond)(c).ppl.tstart(:,ix) = repmat(t',1,length(ix));
        else
            ar.model(m).(data_cond)(c).ppl.tstart(:,ix) = repmat(t,1,length(ix));
        end
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
        ar.model(m).(data_cond)(c).ppl.t = nan(nsteps+1,size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.x_low = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.x_low_vpl = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.ppl_low = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.vpl_low = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.x_high = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.x_high_vpl = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.ppl_high = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.vpl_high = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.ps_low = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), length(ar.p));
        ar.model(m).(data_cond)(c).ppl.ps_high = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), length(ar.p));
        ar.model(m).(data_cond)(c).ppl.gamma = ones(size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2),1);
        ar.model(m).(data_cond)(c).ppl.corr_high = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));
        ar.model(m).(data_cond)(c).ppl.corr_low = nan(nsteps, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2));    
    end
    arSimu(false, true, true);
    if(~onlyProfile)
        for jx=1:length(ix)
            ar.model(m).(data_cond)(c).ppl.t(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.x_low(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.x_low_vpl(:,ix(jx)) = nan;
            if(size(ar.model(m).(data_cond)(c).([x_y 'FineSimu']),2) < ix(jx))
                error('The state you want to calculate does not exist!')
            end
            ar.model(m).(data_cond)(c).ppl.x_orig(:,ix(jx)) = ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(:,ix(jx));
            ar.model(m).(data_cond)(c).ppl.x_high(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.x_high_vpl(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.ps_low(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.ps_high(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.gamma(ix(jx)) = gammas(jx);
            ar.model(m).(data_cond)(c).ppl.ppl_low(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.vpl_low(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.ppl_high(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.vpl_high(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.corr_high(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).ppl.corr_low(:,ix(jx)) = nan; 
            ar.model(m).(data_cond)(c).([x_y 'FineLB'])(:,ix(jx)) = nan;
            ar.model(m).(data_cond)(c).([x_y 'FineUB'])(:,ix(jx)) = nan;
        end
    end
    for jt=1:length(t)
        for jx=1:length(ix)
            cur_t = find(ar.model(m).(data_cond)(c).ppl.tstart(:,ix(jx))==t(jt));
            
            if(isempty(cur_t))               
                ar.model(m).(data_cond)(c).ppl.tstart = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.tstart];
                ar.model(m).(data_cond)(c).ppl.tstart(1,ix(jx)) = t(jt);
                ar.model(m).(data_cond)(c).ppl.xtrial = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1); ar.model(m).(data_cond)(c).ppl.xtrial];
                ar.model(m).(data_cond)(c).ppl.xfit = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1); ar.model(m).(data_cond)(c).ppl.xfit];
                ar.model(m).(data_cond)(c).ppl.vpl = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1); ar.model(m).(data_cond)(c).ppl.vpl ];
                ar.model(m).(data_cond)(c).ppl.ppl = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1); ar.model(m).(data_cond)(c).ppl.ppl ];
                ar.model(m).(data_cond)(c).ppl.ps = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2), 2*n+1, length(ar.p)); ar.model(m).(data_cond)(c).ppl.ps];
                ar.model(m).(data_cond)(c).ppl.lb_fit = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.lb_fit];
                ar.model(m).(data_cond)(c).ppl.ub_fit = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.ub_fit];
                ar.model(m).(data_cond)(c).ppl.lb_fit_vpl = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.lb_fit_vpl];
                ar.model(m).(data_cond)(c).ppl.ub_fit_vpl = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.ub_fit_vpl];
                ar.model(m).(data_cond)(c).ppl.kind_high = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.kind_high];
                ar.model(m).(data_cond)(c).ppl.kind_low = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.kind_low];
                ar.model(m).(data_cond)(c).ppl.kind_high_vpl = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.kind_high_vpl];
                ar.model(m).(data_cond)(c).ppl.kind_low_vpl = [nan(1, size(ar.model(m).(data_cond)(c).([x_y 'ExpSimu']),2)); ar.model(m).(data_cond)(c).ppl.kind_low_vpl];                         
            end
        end
    end

    whichT = find(ar.model(m).(data_cond)(c).ppl.tstart(:,ix(1)) == t(whichT));
    t = ar.model(m).(data_cond)(c).ppl.tstart(:,ix(1));
end