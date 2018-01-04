function [xFit, ps] = xstart_ppl(m, c, jx, t, doPPL, xstd, pReset, chi2start, whichT, takeY, save, dir, xFit, onlyProfile)
global ar;
arWaitbar(0);

if(~exist('dir','var'))
    dir = 0;
end
if(~exist('xFit','var'))
    xFit = 1.;
    save = true;
end
if(~exist('onlyProfile','var'))
    onlyProfile = false;
end
if(takeY)
    data_cond = 'data';
    x_y = 'y';
else
    data_cond = 'condition';
    x_y = 'x';
end
if(~doPPL)
    ppl_vpl = '_vpl';
else
    ppl_vpl = '';
end
ntot = ar.ppl.options.n_steps_profile;
tcount = 1;
for ts = 1:length(t)
    t_tmp = t(ts);
    if(isnan(t_tmp))
        continue;
    end
    %if(save)
    ts_tmp = find(ar.model(m).(data_cond)(c).ppl.tstart(:,jx) == t_tmp);
    if(ts == whichT && save)
        ar.model(m).(data_cond)(c).ppl.t(1,jx) = t(ts);
    end
    %Skip calculation if it already exists
    if(save && ((doPPL && ~isnan(ar.model(m).(data_cond)(c).ppl.ub_fit(ts_tmp,jx))) ...
            || (~doPPL && ~isnan(ar.model(m).(data_cond)(c).ppl.ub_fit_vpl(ts_tmp,jx)))))
       fprintf('The prediction profile you want to compute for t=%d and state %d already exists. If you want to overwrite them, delete the value at ar.model.data/condition.ppl.ub_fit(_vpl) \n',t_tmp,jx)
       continue; 
    end            
   
    %Get value of trajectory at time point
    [~,it_first] = min(abs(ar.model(m).(data_cond)(c).tFine-t_tmp));            
    arLink(true, t_tmp, takeY, jx, c, m, ar.model(m).(data_cond)(c).([x_y 'FineSimu'])(it_first,jx), xstd);
    arCalcMerit(0,[],1);
    [~,it] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));
    if(takeY && length(find(ar.model(m).(data_cond)(c).tExp==t_tmp))>1)
        it = it+1;               
    end
    xSim = ar.model(m).(data_cond)(c).([x_y 'ExpSimu'])(it,jx);
    arLink(true, t_tmp, takeY, jx, c, m, xSim, xstd);
    if(ar.ppl.qLog10)
        xSim = log10(xSim);
    end
    fprintf('Calculating PPL for t=%d and state x=%i \n',t(ts),jx);
    % go up
    if(save || dir==1)
        npre = 0;
        [xtrial_up, xfit_up, ppl_up, vpl_up, ps_up, tcount] = ppl(m, c, jx, it, t_tmp, doPPL, xSim, ...
            chi2start, 1, ar.ppl.qLog10, ntot, xstd, npre, tcount, takeY);    
               
        ar.p = pReset;  
        arCalcMerit();       
    else
        xtrial_up = xSim*1.01;
        xfit_up = xSim*1.01;
        ppl_up = chi2start+0.1;
        vpl_up = chi2start+0.1;
        ps_up = ar.p;
        tcount = 1;
    end
    % go down 
    if(save || dir==-1)
        if(dir==-1)
            npre = 0;
        end
        [xtrial_down, xfit_down, ppl_down, vpl_down, ps_down, tcount] = ppl(m, c, jx, it, t_tmp, doPPL, xSim, ...
            chi2start, -1, ar.ppl.qLog10, ntot, xstd, npre, tcount, takeY);
        % reset parameters
        ar.p = pReset;
    else
        xtrial_down = xSim*0.99;
        xfit_down = xSim*0.99;
        ppl_down = chi2start+0.1;
        vpl_down = chi2start+0.1;
        ps_down = ar.p;
        tcount = 1;
        ar.p = pReset;
    end
    
    %Reset data point
    if(takeY)
        arLink(true,0.,true,jx, c, m,NaN);
    end
    arCalcMerit();
    
    ps_tmp = [flipud(ps_down); pReset; ps_up];
    
    if(doPPL)
        xfit_tmp = [fliplr(xfit_down) xSim xfit_up];
        xtrial_tmp = [fliplr(xtrial_down) xSim xtrial_up];
        ppl_tmp = [fliplr(ppl_down) chi2start ppl_up];
        vpl_tmp = [fliplr(vpl_down) chi2start vpl_up];
        fitted_tmp = [fliplr(xfit_down) xSim xfit_up];
        merit_tmp = [fliplr(ppl_down) chi2start ppl_up];
        q_chi2good = merit_tmp <= chi2start+ar.ppl.dchi2;
        q_nonnan = ~isnan(merit_tmp);         
    else
        xfit_tmp = [fliplr(xfit_down) xSim xfit_up];
        xtrial_tmp = [fliplr(xtrial_down) xSim xtrial_up];
        ppl_tmp = [fliplr(ppl_down) chi2start ppl_up];
        vpl_tmp = [fliplr(vpl_down) chi2start vpl_up];
        fitted_tmp = [fliplr(xtrial_down) xSim xtrial_up];
        merit_tmp = [fliplr(vpl_down) chi2start vpl_up];
        q_chi2good = merit_tmp <= chi2start+ar.ppl.dchi2;
        q_nonnan = ~isnan(merit_tmp); 
    end

    % calculate CI point-wise fit
    lb_tmp = min(fitted_tmp(q_chi2good));
    ub_tmp = max(fitted_tmp(q_chi2good));
    find_tmp = find(merit_tmp > chi2start+ar.ppl.dchi2);

    if(length(vpl_down)==1 || sum(find_tmp<length(vpl_down)+1)==0 || (~save && dir==1))
        lb_tmp = -Inf;             
        kind_low_tmp = -Inf;                
        fprintf('No -95 PPL for t=%d\n',t(ts));
    else
        if(lb_tmp==min(fitted_tmp(q_nonnan)))
            warning(['Multiple likelihood values are assigned to the same model fit. ' ...
                     'check model uncertainty and fits, or set more strict integrator tolerances!'])
        end
        kind_low_tmp = min(find(q_chi2good==1));
        if(length(kind_low_tmp)>1)
            kind_low_tmp = kind_low_tmp(1);
        end
        lb_tmp = interp1(merit_tmp([kind_low_tmp kind_low_tmp-1]), ...
        fitted_tmp([kind_low_tmp kind_low_tmp-1]), chi2start+ar.ppl.dchi2);
    end
    if(length(vpl_up)==1 || sum(find_tmp>length(vpl_down)+1)==0 || (~save && dir==-1))
        ub_tmp = Inf;
        kind_high_tmp = Inf;        
        fprintf('No +95 PPL for t=%d\n',t(ts));
    else
        if(ub_tmp==max(fitted_tmp(q_nonnan)))
            warning(['Multiple likelihood values are assigned to the same model fit. ' ...
                     'check model uncertainty and fits, or set more strict integrator tolerances!'])
        end
        kind_high_tmp = max(find(q_chi2good==1));   
        if(length(kind_high_tmp)>1)
            kind_high_tmp = kind_high_tmp(end);
        end
        ub_tmp = interp1(merit_tmp([kind_high_tmp kind_high_tmp+1]), ...
        fitted_tmp([kind_high_tmp kind_high_tmp+1]), chi2start+ar.ppl.dchi2);
    end      
    
    if ~onlyProfile
        if((dir==1 && isinf(ub_tmp)) || (dir==-1 && isinf(lb_tmp)))
            fprintf('ERROR: no bound found at t=%d \n', t_tmp);
            if(takeY)
                arLink(true,0.,true,ix(jx), c, m,NaN);
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
    if(save)
        ar.model(m).(data_cond)(c).ppl.xtrial(ts_tmp, jx,:) = xtrial_tmp;
        ar.model(m).(data_cond)(c).ppl.xfit(ts_tmp, jx,:) = xfit_tmp;
        ar.model(m).(data_cond)(c).ppl.ppl(ts_tmp, jx,:) = ppl_tmp;
        ar.model(m).(data_cond)(c).ppl.vpl(ts_tmp, jx,:) = vpl_tmp;
        ar.model(m).(data_cond)(c).ppl.ps(ts_tmp, jx,:,:) = ps_tmp;
        ar.model(m).(data_cond)(c).ppl.(['lb_fit' ppl_vpl])(ts_tmp, jx) = lb_tmp;
        ar.model(m).(data_cond)(c).ppl.(['ub_fit' ppl_vpl])(ts_tmp, jx) = ub_tmp;
        ar.model(m).(data_cond)(c).ppl.(['kind_high' ppl_vpl])(ts_tmp, jx) = kind_high_tmp;
        ar.model(m).(data_cond)(c).ppl.(['kind_low' ppl_vpl])(ts_tmp, jx) = kind_low_tmp;                
    end
end
%write LB/UB in ar struct
    if(save && onlyProfile && dir==0 && length(t)>1)        
        struct_vec = {'FineUB','FineLB'};
        low_high_vec = {'ub_fit','lb_fit'};
        for ilh=1:2
            struct_string = struct_vec{ilh};
            ppl_string = [low_high_vec{ilh} ppl_vpl];
            if(takeY)
                struct_string = ['y' struct_string];
            else
                struct_string = ['x' struct_string];
            end            
            ar.model(m).(data_cond)(c).(struct_string)(ar.model(m).(data_cond)(c).tFine<=max(ar.model(m).(data_cond)(c).ppl.tstart(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)))),jx) = ...
                        interp1(ar.model(m).(data_cond)(c).ppl.tstart(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
                        ar.model(m).(data_cond)(c).ppl.(ppl_string)(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
                        ar.model(m).(data_cond)(c).tFine(ar.model(m).(data_cond)(c).tFine<=max(ar.model(m).(data_cond)(c).ppl.tstart(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx))))),...
                        'pchip',NaN);            
        end
    end
end


function [xtrial, xfit, ppl, vpl, ps, tcount] = ppl(m, c, ix, it, t, doPPL, xSim,  chi2start, direction, qLog10, n, xstd, npre, tcount, takeY)

global ar

xtrial = nan(1,n);
xfit = nan(1,n);
ppl = nan(1,n);
vpl = nan(1,n);
ps = nan(n,length(ar.p));

dx = sqrt(ar.ppl.dchi2*ar.ppl.options.rel_increase) * xstd;

if(takeY)
    xLabel = myNameTrafo(ar.model(m).data(c).y{ix});
else
    xLabel = myNameTrafo(ar.model(m).x{ix});    
end
xExp = xSim;
for j = 1:n
    if(toc>tcount)
        if(direction>0)
            arWaitbar((j+npre), n, sprintf('PPL (up) for %s at t=%g %i/%i', xLabel, t, j, n));
        else
            arWaitbar((j+npre), n, sprintf('PPL (down) for %s at t=%g %i/%i', xLabel, t, j, n));
        end
        tcount = tcount + 0.5; % update every half second
    end
    
    xExp = xExp + direction*dx;
    arLink(true,t,takeY,ix, c,m,xExp,xstd);
    
    try
        arPPLFit;
    catch exception
        fprintf('ERROR in PPL integration (%s)\n', exception.message);
        if(takeY)
            arLink(true,0.,true,ix, c, m,NaN);
        end
        break;
    end
    
    xtrial(j) = xExp;    
    arCalcMerit(0, ar.p(ar.qFit==1),1)
    if(takeY)
        xSim = ar.model(m).data(c).yExpSimu(it,ix); 
    else
        xSim = ar.model(m).condition(c).xExpSimu(it,ix); 
    end
    xfit(j) = xSim;
    if(takeY)
        ppl(j) = arGetMerit('chi2') - (xtrial(j)-xfit(j)).^2/xstd.^2;
        vpl(j) = arGetMerit('chi2');
        
    else
        ppl(j) = arGetMerit('chi2');
        vpl(j) = arGetMerit('chi2') + (xtrial(j)-xfit(j)).^2/xstd.^2;
    end
    ps(j,:) = ar.p;
    %ppl(j)
    %xtrial(j)
    if((doPPL && ppl(j) > chi2start+ar.ppl.dchi2*1.2) || (~doPPL && vpl(j) > chi2start+ar.ppl.dchi2*1.2))
        break
    end
end

    function arPPLFit
        
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
            
                sxSim = zeros(1,length(ar.p));           
                sxSim(ar.model(m).condition(c).pLink) = ...
                        squeeze(ar.model(m).condition(c).sxExpSimu(it,ix,:))';
                for j10=find(ar.qLog10==1)
                    sxSim(j10) = sxSim(j10) * 10.^ar.p(j10) * log(10);
                end
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