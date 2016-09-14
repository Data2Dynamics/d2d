function [t, whichT] = PPL_init(m,c,t,ix,gammas, onlyProfile, whichT,takeY)
    global ar;
    n = ar.ppl.n;
    nsteps = ar.ppl.nsteps;
    if(~takeY && (~isfield(ar.model(m).condition(c),'ppl') || (isfield(ar.model(m).condition(c),'ppl') && isempty(ar.model(m).condition(c).ppl)) ...
        ))%|| (isfield(ar.model(m).condition(c),'ppl') && ~isequal(ar.model(m).condition(c).ppl.tstart(:,ix(1))',t))))
        ar.model(m).condition(c).ppl.ix = ix(:);   
        ar.model(m).condition(c).ppl.xtrial = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2), 2*n+1);
        ar.model(m).condition(c).ppl.xfit = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2), 2*n+1);
        ar.model(m).condition(c).ppl.vpl = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2), 2*n+1);
        ar.model(m).condition(c).ppl.ppl = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2), 2*n+1);
        ar.model(m).condition(c).ppl.ps = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2), 2*n+1, length(ar.p));
        ar.model(m).condition(c).ppl.lb_fit = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.ub_fit = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.lb_fit_vpl = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.ub_fit_vpl = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.kind_high = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.kind_low = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.kind_high_vpl = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.kind_low_vpl = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.tstart = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
        if(size(t,1)==1)
            ar.model(m).condition(c).ppl.tstart(:,ix) = repmat(t',1,length(ix));
        else
            ar.model(m).condition(c).ppl.tstart(:,ix) = repmat(t,1,length(ix));
        end
        ar.model(m).condition(c).xFineLB = nan(length(ar.model(m).condition(c).tFine),size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).xFineUB = nan(length(ar.model(m).condition(c).tFine),size(ar.model(m).condition(c).xExpSimu,2));

        ar.model(m).condition(c).ppl.t = nan(nsteps+1,size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.x_low = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.x_low_vpl = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.ppl_low = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.vpl_low = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.x_high = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.x_high_vpl = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.ppl_high = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.vpl_high = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.ps_low = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2), length(ar.p));
        ar.model(m).condition(c).ppl.ps_high = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2), length(ar.p));
        ar.model(m).condition(c).ppl.gamma = ones(size(ar.model(m).condition(c).xExpSimu,2),1);
        ar.model(m).condition(c).ppl.corr_high = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));
        ar.model(m).condition(c).ppl.corr_low = nan(nsteps, size(ar.model(m).condition(c).xExpSimu,2));    
    elseif(takeY && (~isfield(ar.model(m).data(c),'ppl') || (isfield(ar.model(m).data(c),'ppl') && isempty(ar.model(m).data(c).ppl)) ...
            ))%|| (isfield(ar.model(m).data(c),'ppl') && ~isequal(ar.model(m).data(c).ppl.tstart(:,ix(1))',t))))
            ar.model(m).data(c).ppl.ix = ix(:);   
            ar.model(m).data(c).ppl.xtrial = nan(length(t), size(ar.model(m).data(c).yExpSimu,2), 2*n+1);
            ar.model(m).data(c).ppl.xfit = nan(length(t), size(ar.model(m).data(c).yExpSimu,2), 2*n+1);
            ar.model(m).data(c).ppl.ppl = nan(length(t), size(ar.model(m).data(c).yExpSimu,2), 2*n+1);
            ar.model(m).data(c).ppl.vpl = nan(length(t), size(ar.model(m).data(c).yExpSimu,2), 2*n+1);
            ar.model(m).data(c).ppl.ps = nan(length(t), size(ar.model(m).data(c).yExpSimu,2), 2*n+1, length(ar.p));
            ar.model(m).data(c).ppl.lb_fit = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.ub_fit = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.lb_fit_vpl = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.ub_fit_vpl = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.kind_high = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.kind_low = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.kind_high_vpl = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.kind_low_vpl = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.tstart = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
            if(size(t,1)==1)
                ar.model(m).data(c).ppl.tstart(:,ix) = repmat(t',1,length(ix));
            else
                ar.model(m).data(c).ppl.tstart(:,ix) = repmat(t,1,length(ix));
            end
            
            ar.model(m).data(c).yFineLB = nan(length(ar.model(m).data(c).tFine),size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).yFineUB = nan(length(ar.model(m).data(c).tFine),size(ar.model(m).data(c).yExpSimu,2));

            ar.model(m).data(c).ppl.t = nan(nsteps+1,size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.x_low = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.x_low_vpl = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.ppl_low = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.vpl_low = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.x_high = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.x_high_vpl = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.ppl_high = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.vpl_high = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.ps_low = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2), length(ar.p));
            ar.model(m).data(c).ppl.ps_high = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2), length(ar.p));
            ar.model(m).data(c).ppl.gamma = ones(size(ar.model(m).data(c).yExpSimu,2),1);
            ar.model(m).data(c).ppl.corr_high = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
            ar.model(m).data(c).ppl.corr_low = nan(nsteps, size(ar.model(m).data(c).yExpSimu,2));
    end
    arSimu(false, true, true);
    if(~onlyProfile)
        for jx=1:length(ix)
            if(~takeY)

                ar.model(m).condition(c).ppl.t(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.x_low(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.x_low_vpl(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.x_orig(:,ix(jx)) = ar.model(m).condition(c).xFineSimu(:,ix(jx));
                ar.model(m).condition(c).ppl.x_high(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.x_high_vpl(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.ps_low(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.ps_high(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.gamma(ix(jx)) = gammas(jx);
                ar.model(m).condition(c).ppl.ppl_low(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.vpl_low(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.ppl_high(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.vpl_high(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.corr_high(:,ix(jx)) = nan;
                ar.model(m).condition(c).ppl.corr_low(:,ix(jx)) = nan;
            else

                ar.model(m).data(c).ppl.t(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.x_low(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.x_low_vpl(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.x_orig(:,ix(jx)) = ar.model(m).data(c).yFineSimu(:,ix(jx));
                ar.model(m).data(c).ppl.x_high(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.x_high_vpl(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.ps_low(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.ps_high(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.gamma(ix(jx)) = gammas(jx);
                ar.model(m).data(c).ppl.ppl_low(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.vpl_low(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.ppl_high(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.vpl_high(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.corr_high(:,ix(jx)) = nan;
                ar.model(m).data(c).ppl.corr_low(:,ix(jx)) = nan;
            end
        end
    end
    for jt=1:length(t)
        for jx=1:length(ix)
            if(~takeY)
                cur_t = find(ar.model(m).condition(c).ppl.tstart(:,ix(jx))==t(jt));
            elseif(takeY)
                cur_t = find(ar.model(m).data(c).ppl.tstart(:,ix(jx))==t(jt));
            end
            if(isempty(cur_t))
                if(~takeY)
                    ar.model(m).condition(c).ppl.tstart = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.tstart];
                    ar.model(m).condition(c).ppl.tstart(1,ix(jx)) = t(jt);
                    ar.model(m).condition(c).ppl.xtrial = [nan(1, size(ar.model(m).condition(c).xExpSimu,2), 2*n+1); ar.model(m).condition(c).ppl.xtrial];
                    ar.model(m).condition(c).ppl.xfit = [nan(1, size(ar.model(m).condition(c).xExpSimu,2), 2*n+1); ar.model(m).condition(c).ppl.xfit];
                    ar.model(m).condition(c).ppl.vpl = [nan(1, size(ar.model(m).condition(c).xExpSimu,2), 2*n+1); ar.model(m).condition(c).ppl.vpl ];
                    ar.model(m).condition(c).ppl.ppl = [nan(1, size(ar.model(m).condition(c).xExpSimu,2), 2*n+1); ar.model(m).condition(c).ppl.ppl ];
                    ar.model(m).condition(c).ppl.ps = [nan(1, size(ar.model(m).condition(c).xExpSimu,2), 2*n+1, length(ar.p)); ar.model(m).condition(c).ppl.ps];
                    ar.model(m).condition(c).ppl.lb_fit = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.lb_fit];
                    ar.model(m).condition(c).ppl.ub_fit = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.ub_fit];
                    ar.model(m).condition(c).ppl.lb_fit_vpl = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.lb_fit_vpl];
                    ar.model(m).condition(c).ppl.ub_fit_vpl = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.ub_fit_vpl];
                    ar.model(m).condition(c).ppl.kind_high = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.kind_high];
                    ar.model(m).condition(c).ppl.kind_low = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.kind_low];
                    ar.model(m).condition(c).ppl.kind_high_vpl = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.kind_high_vpl];
                    ar.model(m).condition(c).ppl.kind_low_vpl = [nan(1, size(ar.model(m).condition(c).xExpSimu,2)); ar.model(m).condition(c).ppl.kind_low_vpl];          
                else
                    ar.model(m).data(c).ppl.tstart = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.tstart];
                    ar.model(m).data(c).ppl.tstart(1,ix(jx)) = t(jt);
                    ar.model(m).data(c).ppl.xtrial = [nan(1, size(ar.model(m).data(c).yExpSimu,2), 2*n+1); ar.model(m).data(c).ppl.xtrial];
                    ar.model(m).data(c).ppl.xfit = [nan(1, size(ar.model(m).data(c).yExpSimu,2), 2*n+1); ar.model(m).data(c).ppl.xfit];
                    ar.model(m).data(c).ppl.vpl = [nan(1, size(ar.model(m).data(c).yExpSimu,2), 2*n+1); ar.model(m).data(c).ppl.vpl ];
                    ar.model(m).data(c).ppl.ppl = [nan(1, size(ar.model(m).data(c).yExpSimu,2), 2*n+1); ar.model(m).data(c).ppl.ppl ];
                    ar.model(m).data(c).ppl.ps = [nan(1, size(ar.model(m).data(c).yExpSimu,2), 2*n+1, length(ar.p)); ar.model(m).data(c).ppl.ps];
                    ar.model(m).data(c).ppl.lb_fit = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.lb_fit];
                    ar.model(m).data(c).ppl.ub_fit = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.ub_fit];
                    ar.model(m).data(c).ppl.lb_fit_vpl = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.lb_fit_vpl];
                    ar.model(m).data(c).ppl.ub_fit_vpl = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.ub_fit_vpl];
                    ar.model(m).data(c).ppl.kind_high = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.kind_high];
                    ar.model(m).data(c).ppl.kind_low = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.kind_low];
                    ar.model(m).data(c).ppl.kind_high_vpl = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.kind_high_vpl];
                    ar.model(m).data(c).ppl.kind_low_vpl = [nan(1, size(ar.model(m).data(c).yExpSimu,2)); ar.model(m).data(c).ppl.kind_low_vpl];
                end
            end
        end
    end
    if(~takeY)
        whichT = find(ar.model(m).condition(c).ppl.tstart(:,ix(1)) == t(whichT));
        t = ar.model(m).condition(c).ppl.tstart(:,ix(1));
    else          
        whichT = find(ar.model(m).data(c).ppl.tstart(:,ix(1)) == t(whichT));
        t = ar.model(m).data(c).ppl.tstart(:,ix(1));
    end
end