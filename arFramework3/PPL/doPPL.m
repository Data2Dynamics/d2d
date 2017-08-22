% calculate prediction profile likelihood
% via explicit integration step followed by optimization if needed
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
    ar.config.SimuPPL=1;
    ar.ppl.qLog10=0;

    if nargin < 6
      options = [];
    end
    if nargin < 2
        fprintf('Please specify the model, condition/data index, which states should be integrated \n, the time points and if an observation or internal state should be used. \n');
        return;  
    end

    if(~exist('takeY','var') || isempty(takeY))
        takeY = true;
        fprintf('Not specified whether PBs on observation or internal state should be calculated. \n If not specified, observation is taken!\n');
    end

    if(~exist('ix','var') || isempty(ix))
        fprintf('No specific state given, thus all are taken!\n');
        if(takeY)
            ix = 1:length(ar.model(m).data(c).yNames);
        else
            ix = 1:length(ar.model(m).xNames);
        end
    end

    if(~isfield(ar.ppl,'xstd_auto'))
        ar.ppl.xstd_auto = 0;
    end
    confirm_options = PPL_options(options);
    fprintf(confirm_options)
    ar.ppl.n = ar.ppl.options.n_start;
    if(~isfield(ar.ppl.options,'tEnd'))
        if(takeY)
           ar.ppl.options.tEnd = ar.model(m).data(c).tFine(end);
        else
           ar.ppl.options.tEnd = ar.model(m).condition(c).tFine(end); 
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

    %Set correction strength
    if(~isfield(ar.ppl.options,'gammas'))
        ar.ppl.options.gammas = ones(size(ix))*1./ar.ppl.options.stepsize;
    end
    if(length(ar.ppl.options.gammas) == 1 || length(ar.ppl.options.gammas)~=length(ix))
        ar.ppl.options.gammas = repmat(ar.ppl.options.gammas(1),1,length(ix));
    end
    if(length(ar.ppl.options.gammas) ~= length(ix))
        error('Argument gammas has an incorrect length');
    end

    % optimizer settings (set only once)
    if(ar.config.fiterrors~=-1 && ~ar.ppl.options.onlyProfile)
        ar.ppl.fittederrors=ar.config.fiterrors;
        ar.config.fiterrors=0;
        ar.ppl.fit_bkp = ar.qFit(strncmp(ar.pLabel,'sd_',3));
        ar.qFit(strncmp(ar.pLabel,'sd_',3))=2;
    end
    if(ar.config.useSensis)
        ar.config.optim.Jacobian = 'on';
    else
        ar.config.optim.Jacobian = 'off';
    end
    ar.config.optim.Algorithm = 'trust-region-reflective';

    ar.ppl.dchi2 = chi2inv(1-ar.ppl.options.alpha_level, 1);
    ar.ppl.dchi2;
    dir = ar.ppl.options.dir;
    arCalcMerit();

    ar.ppl.chi2_95 = arGetMerit('chi2') + arGetMerit('chi2err')+ar.ppl.dchi2 + 0.5;

    % Initialize PPL struct, set values to NaN that are newly computed
    [t, whichT] = PPL_init(m,c,t,ix,ar.ppl.options.gammas, ar.ppl.options.onlyProfile, ar.ppl.options.whichT,takeY);

    pReset = ar.p;
    chi2start = arGetMerit('chi2') + arGetMerit('chi2err') ;

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
    
    for jx = 1:length(ix)
        if(ar.ppl.xstd_auto)    
            if(ar.config.fiterrors~= -1 && takeY && ~isnan(ar.model(m).data(c).ystdFineSimu(1,ix(jx))))
                xstd = ar.model(m).data(c).ystdFineSimu(1,ix(jx));
            elseif(ar.config.fiterrors==-1 && takeY)
                [~,it_first] = min(abs(ar.model(m).data(c).tExp-t(whichT))); 
                 if(~isnan(ar.model(m).data(c).yExpStd(it_first,ix(jx))))
                     xstd = ar.model(m).data(c).yExpStd(it_first,ix(jx));
                 end
            elseif(~takeY)
                arSimu(0,1,1);
                xstd = max(ar.model(m).condition(c).xFineSimu(:,ix(jx)))/10;
            end

        end
        xstart_ppl(m, c, ix(jx), t, ar.ppl.options.doPPL, ar.ppl.options.xstd, pReset, chi2start, whichT, takeY, true, 1, [], ar.ppl.options.onlyProfile);
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
        %loop through time  
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
        arWaitbar(-1);
        toc;
        if(takeY)
            arLink(true,0.,true,ix(jx), c, m,NaN);
        end
        ar.p = pReset;

    end
    if(~ar.config.fiterrors && ar.ppl.fittederrors  && ~ar.ppl.options.onlyProfile) 
        ar.config.fiterrors=ar.ppl.fittederrors;
        ar.qFit(strncmp(ar.pLabel,'sd_',3))=ar.ppl.fit_bkp;
    end
    ar.config.SimuPPL=0;
    arCalcMerit(); 
    if(~takeY)
        ar.model(m).qPlotXs(ar.model(m).condition(c).dLink(1))=1;        
    else
        ar.model(m).qPlotYs(c) = 1;        
    end
end

function [xFit, ps] = xstart_ppl(m, c, jx, t, doPPL, xstd, pReset, chi2start, whichT, takeY, save, dir, xFit, onlyProfile)
global ar;

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
arWaitbar(0);
ntot = ar.ppl.options.n_start;
tcount = 1;
for ts = 1:length(t)
    t_tmp = t(ts);
    %if(save)
    ts_tmp = find(ar.model(m).(data_cond)(c).ppl.tstart(:,jx) == t_tmp);
    if(ts == whichT && save)
        ar.model(m).(data_cond)(c).ppl.t(1,jx) = t(ts);
    end
    if(save && ((doPPL && ~isnan(ar.model(m).(data_cond)(c).ppl.ub_fit(ts_tmp,jx))) ...
            || (~doPPL && ~isnan(ar.model(m).(data_cond)(c).ppl.ub_fit_vpl(ts_tmp,jx)))))
       continue; 
    end            
   
    [~,it_first] = min(abs(ar.model(m).(data_cond)(c).tExp-t_tmp));            
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
            chi2start, 1, ar.ppl.qLog10, ar.ppl.options.n_start, xstd, npre, ntot, tcount, takeY);    
        % go down        
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
    
    if(save || dir==-1)
        if(dir==-1)
            npre = 0;
        end
        [xtrial_down, xfit_down, ppl_down, vpl_down, ps_down, tcount] = ppl(m, c, jx, it, t_tmp, doPPL, xSim, ...
            chi2start, -1, ar.ppl.qLog10, ar.ppl.options.n_start, xstd, npre, ntot, tcount, takeY);
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

    if(lb_tmp==min(fitted_tmp(q_nonnan)))
        lb_tmp = -Inf;                
        fprintf('No -95 PPL for t=%d\n',t(ts));
    else
        kind_low_tmp = find(fitted_tmp==lb_tmp);
        if(length(kind_low_tmp)>1)
            kind_low_tmp = kind_low_tmp(1);
        end
        lb_tmp = interp1(merit_tmp([kind_low_tmp kind_low_tmp-1]), ...
        fitted_tmp([kind_low_tmp kind_low_tmp-1]), chi2start+ar.ppl.dchi2);
    end
    if(ub_tmp==max(fitted_tmp(q_nonnan)))
        ub_tmp = Inf;
        fprintf('No +95 PPL for t=%d\n',t(ts));
    else
        kind_high_tmp = find(fitted_tmp==ub_tmp);   
        if(length(kind_high_tmp)>1)
            kind_high_tmp = kind_high_tmp(end);
        end
        ub_tmp = interp1(merit_tmp([kind_high_tmp kind_high_tmp+1]), ...
        fitted_tmp([kind_high_tmp kind_high_tmp+1]), chi2start+ar.ppl.dchi2);
    end        
    
    if ~onlyProfile
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
end


function [xtrial, xfit, ppl, vpl, ps, tcount] = ppl(m, c, ix, it, t, doPPL, xSim,  chi2start, direction, qLog10, n, xstd, npre, ntot, tcount, takeY)

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
            arWaitbar((j+npre), ntot, sprintf('PPL (up) for %s at t=%g %i/%i', xLabel, t, j, n));
        else
            arWaitbar((j+npre), ntot, sprintf('PPL (down) for %s at t=%g %i/%i', xLabel, t, j, n));
        end
        tcount = tcount + 0.5; % update every half second
    end
    
    xExp = xExp + direction*dx;
    try
        arPPLFit;
    catch exception
        fprintf('ERROR PPL: going to lower bound (%s)\n', exception.message);
        break;
    end
    
    xtrial(j) = xExp;    
    arLink(true,t,takeY,ix, c,m,xExp,xstd);
    arCalcMerit(0, ar.p(ar.qFit==1),1)
    if(takeY)
        xSim = ar.model(m).data(c).yExpSimu(it,ix); 
    else
        xSim = ar.model(m).condition(c).xExpSimu(it,ix); 
    end
    xfit(j) = xSim;
    if(takeY)
        ppl(j) = arGetMerit('chi2') + arGetMerit('chi2err') - (xtrial(j)-xfit(j)).^2/xstd.^2;
        vpl(j) = arGetMerit('chi2') + arGetMerit('chi2err');
        
    else
        ppl(j) = arGetMerit('chi2') + arGetMerit('chi2err');
        vpl(j) = arGetMerit('chi2') + arGetMerit('chi2err') + (xtrial(j)-xfit(j)).^2/xstd.^2;
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
        
        arLink(true,t,takeY,ix, c,m,xExp,xstd);
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

function ppl_calc(m, c, jx, xFit, p, t, doPPL, takeY, dir, stepsize, xstd, ed_steps, pReset, chi2start, backward, fineInt)
    global ar
    ar.ppl.xFit_tmp = xFit;
    chi2 = chi2start + ar.ppl.dchi2;
    xSim2 = NaN;
    xSim = NaN;
    xSim3 = NaN;
    npre=0;
    if(takeY)
        data_cond = 'data';
    else
        data_cond = 'condition';
    end
          
    if(dir==1)
        high_low = '_high';
    else
        high_low = '_low';
    end
    if(~doPPL)
        ppl_vpl = '_vpl';
    else
        ppl_vpl = '';
    end
    if(isnan(xFit))
        error('Starting value is NaN, check calculation of profiles! \n')
    end
    t_dir = 1;
    t_tmp=t;
    jt = 0;
    if(backward)
        jt = ar.ppl.nsteps-1;
        t_dir = -1;
    end
    qLog10 = ar.ppl.qLog10;
    if(takeY)
        xLabel = myNameTrafo(ar.model(m).data(c).y{jx});
        gamma_tmp = ar.model(m).data(c).ppl.gamma(jx);
    else
        xLabel = myNameTrafo(ar.model(m).x{jx});
        gamma_tmp = ar.model(m).condition(c).ppl.gamma(jx);
    end    
    nsteps=ar.ppl.nsteps;
    ar.p=squeeze(p);
    if(size(ar.p,2)~=1)
        ar.p=ar.p';
    end
    tcount = 1;
    arWaitbar(0);
    while jt<nsteps
        jt = jt+1;
        corr_tmp = 0;

        [~,it_orig] = min(abs(ar.model(m).(data_cond)(c).tFine-t_tmp-t_dir*stepsize));
        it_getfRHS = it_orig;
        x1_orig=ar.model(m).(data_cond)(c).ppl.x_orig(it_orig,jx);

        if(ar.model(m).(data_cond)(c).tFine(it_orig)-t_tmp-t_dir*stepsize > 0)
            x2_orig=x1_orig;   
            if(it_orig>1);
                x1_orig=ar.model(m).(data_cond)(c).ppl.x_orig(it_orig-1,jx);
            end
        else
            it_orig=it_orig+1;
            if(size(ar.model(m).(data_cond)(c).ppl.x_orig,1)>=it_orig);
                x2_orig=ar.model(m).(data_cond)(c).ppl.x_orig(it_orig,jx);
            else
                x2_orig = x1_orig;
                it_orig=it_orig-1;
            end
        end
        if(it_orig>1)
            x_orig = x1_orig + (x2_orig-x1_orig)/(ar.model(m).(data_cond)(c).tFine(it_orig)-ar.model(m).(data_cond)(c).tFine(it_orig-1))*(t_tmp+t_dir*stepsize-ar.model(m).(data_cond)(c).tFine(it_orig-1));
        else
           x_orig = x1_orig; 
        end

        if(toc>tcount)        
            if(dir==1)
                string_tmp = 'upper';
            else
                string_tmp = 'lower';            
            end
            arWaitbar((jt+npre), nsteps, sprintf(['PPL-integration (' string_tmp ' bound) for %s at t=%g %i/%i'], xLabel, t_tmp, jt, nsteps));
            tcount = tcount + 0.5; % update every half second
        end        

        if(ed_steps==true)% && ~doPPL)               
            %VPL part    
            if(~doPPL)
               [chi2, xSim, xSim2, xSim3, it] = PPL_chi2(t_tmp, true, m, c, jx, takeY, qLog10, doPPL, t_dir*stepsize, ar.ppl.xFit_tmp, xstd);        
               [dps, dx, gamma_tmp] = VPL_Int(t_tmp, m, c, jx, takeY, qLog10, t_dir*stepsize, ar.ppl.xFit_tmp, xstd, gamma_tmp, fineInt); 
            else          
                [dps, dx, gamma_tmp] = PPL_Int(t_tmp, m, c, jx, takeY, qLog10, t_dir*stepsize, ar.ppl.xFit_tmp, xstd, dir, gamma_tmp); 
            end              

            region_fac=1;
            if((ar.ppl.xFit_tmp + dx(1)*stepsize < x_orig && dir==1) || (ar.ppl.xFit_tmp + dx(1)*stepsize > x_orig && dir==-1))
               region_fac =  abs(ar.ppl.xFit_tmp - x_orig) / (2* abs(dx(1)) *stepsize);
            end

            region_fac_p=ones(1,length(dps));

            if(sum(ar.p(ar.qFit==1) + dps'*stepsize > ar.ub(ar.qFit==1))>0 || (sum(ar.p(ar.qFit==1) + dps'*stepsize < ar.lb(ar.qFit==1))>0))
%                 if(trust_radius)
                    region_fac_p = min(abs((ar.p(ar.qFit==1) - ar.ub(ar.qFit==1))./ (2*dps'*stepsize)));
%                 else
%                     region_fac_p = abs( (ar.p(ar.qFit==1) - ar.ub(ar.qFit==1)) ./ (2*dps'*stepsize));
%                     region_fac_p(region_fac_p>1) = 1;
%                 end
                    still_bad = abs((ar.p(ar.qFit==1) - ar.lb(ar.qFit==1)) ./ (2*dps'.*region_fac_p*stepsize)) < region_fac_p;
               if(sum( still_bad > 0))
                    region_low = abs( (ar.p(ar.qFit==1) - ar.lb(ar.qFit==1)) ./ (2*dps'*stepsize));
%                     if(trust_radius)
                        region_fac_p = min(region_low);
%                     else
%                         region_fac_p(still_bad) = region_low(still_bad);
%                     end
               end
            end

            if(max(region_fac_p) < 0.1)
                region_fac_p = zeros(1,length(dps));
            end     
            ar.p(ar.qFit==1)=ar.p(ar.qFit==1) + dps'*stepsize.*region_fac_p;
            ar.ppl.xFit_tmp = ar.ppl.xFit_tmp + dx(1)*stepsize*region_fac;
            ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
            ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);            
            t_tmp = t_tmp + t_dir*stepsize;  
            if(doPPL)
                [chi2, xSim] = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, doPPL, t_dir*stepsize, ar.ppl.xFit_tmp, xstd);
                ar.ppl.xFit_tmp = xSim;
            end
            last_corr = 0;
            if(jt>10)
                last_corr = sum(ar.model(m).(data_cond)(c).ppl.(['corr' high_low])(jt-10:jt,jx));
            end
            if(((fineInt || (max(abs(dps))<1.e-4 || ~isempty(find(~region_fac_p,1))) ...
                    || (ar.ppl.xFit_tmp<x_orig && dir==1) || (ar.ppl.xFit_tmp>x_orig && dir==-1)) ...
                    && last_corr==0) || mod(jt/(floor(nsteps/10)),1)==0)%  || max(abs(dps))==0 || (doPPL && ((curve_xSim>0 && dir ==1) || (curve_xSim<0 && dir==-1))))
                [chi2, xSim] = PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL,chi2);
                corr_tmp = 1;
                if(doPPL)
                    ar.ppl.xFit_tmp = xSim;
                    if(abs(chi2-ar.ppl.chi2_95 + 0.5)>0.5 )
                       chi2 = PPL_doSim_calc(chi2, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL);              
                    end
                end
            elseif(~doPPL)
                [chi2, xSim] = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, doPPL, stepsize, ar.ppl.xFit_tmp, xstd);
            end

            %[chi2, xSim] = PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, 0);
            if(~doPPL && abs(chi2-ar.ppl.chi2_95 + 0.5)>1 )
                fprintf('trying parallel step at t=%i \n',t_tmp-t_dir*stepsize);

                ar.p(ar.qFit==1)=ar.p(ar.qFit==1) - dps'*stepsize.*region_fac_p;
                t_tmp = t_tmp - t_dir*stepsize;
                if(doPPL)
                   [~, ar.ppl.xFit_tmp] = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, doPPL, stepsize, ar.ppl.xFit_tmp, xstd); 
                else
                    ar.ppl.xFit_tmp = ar.ppl.xFit_tmp - dx(1)*stepsize*region_fac;
                end
                ar.ppl.xFit_tmp = getxFit(ar.ppl.xFit_tmp);
                %ar.ppl.xFit_tmp = xSim2;
                t_tmp = t_tmp + t_dir*stepsize;
                [chi2, ar.ppl.xFit_tmp] = PPL_corr(4, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL, chi2);    
            end            
            ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
            ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);        
            if(sum(isnan(dx))>0)               
               ar.p=pReset;
               if(takeY)
                    arLink(true,ar.model(m).data(c).tExp(1),true,jx, c, m,NaN);
               end
               fprintf('ERROR IN STEP AT T=%d \n', t_tmp);
               return;
            end            
        else

            xSim = getxFit(ar.ppl.xFit_tmp);
            t_tmp = t_tmp + t_dir*stepsize;
            ar.ppl.xFit_tmp = xSim;
            chi2 = PPL_doSim_calc(chi2, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL); 

        end

        if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 || (ar.ppl.xFit_tmp<x_orig && dir==1) || (ar.ppl.xFit_tmp>x_orig && dir==-1))

            i_count=0;
            while((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )  
                fprintf('Correction at t=%d \n',t_tmp);
                corr_tmp = 1;
                if(i_count==1 && doPPL)
                    PPL_corr(3, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL,0);
                elseif(i_count==1 && ~doPPL)
                    PPL_corr(2, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL,0);
                end
                chi2 = PPL_doSim_calc(chi2, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL);

                i_count = i_count + 1 ;
                if((~fineInt && (doPPL && i_count>2) || (~doPPL && i_count>4)) || ...
                   (fineInt && (doPPL && i_count>5) || (~doPPL && i_count>6)))
                    break;
                end
            end

            if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
                fprintf('Try to restart at t=%d diff in chi2 is %d \n',t_tmp,abs(chi2 - ar.ppl.chi2_95 + 0.5));
                ar.p = pReset;
                arLink(true,0.,takeY,jx, c, m,NaN);
                [ar.ppl.xFit_tmp, ar.p] = xstart_ppl(m, c, jx, t_tmp, doPPL, xstd, pReset, chi2start, 10, takeY, false, dir, ar.ppl.xFit_tmp);                
                chi2 = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, doPPL, stepsize, ar.ppl.xFit_tmp, xstd);               
                if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
                    fprintf('Check the Profile at t=%d for inconsistencies, diff in chi2 is %d \n',t_tmp,abs(chi2 - ar.ppl.chi2_95 + 0.5));
                end
            end
        end       

        if(doPPL)
            ppl_vpl = 'ppl';
            ar.model(m).(data_cond)(c).ppl.(['x' high_low])(jt+t_dir*1, jx)=xSim;
        else
            ppl_vpl = 'vpl';
            ar.model(m).(data_cond)(c).ppl.(['x' high_low '_' ppl_vpl])(jt+t_dir*1, jx)=ar.ppl.xFit_tmp;
        end

        ar.ppl.chi2_tmp = chi2;
        ar.model(m).(data_cond)(c).ppl.(['corr' high_low])(jt+t_dir*1, jx)=corr_tmp;
        ar.model(m).(data_cond)(c).ppl.([ppl_vpl high_low])(jt+t_dir*1, jx)=chi2;
        ar.model(m).(data_cond)(c).ppl.t(jt+t_dir*1,jx)=t_tmp;
        ar.model(m).(data_cond)(c).ppl.(['ps' high_low])(jt+t_dir*1, jx,:)=ar.p;  

        if(t_dir==-1 && jt>1)
            jt=jt-2;
        elseif(t_dir==-1 && jt<=1)
            fprintf('Backward integration stopped because lower time bound hit, proceeding with normal integration \n');
            break
        end

    end
    %write LB/UB in ar struct
    if(dir==1)
        struct_string = 'FineUB';
        ppl_string = 'x_high';
    else
        struct_string = 'FineLB';
        ppl_string = 'x_low';
    end
    if(~doPPL)
        ppl_string = [ppl_string '_vpl'];
    end
    if(takeY)
        data_cond = 'data';
        struct_string = ['y' struct_string];
    else
        data_cond = 'condition';
        struct_string = ['x' struct_string];
    end
    ar.model(m).(data_cond)(c).(struct_string)(ar.model(m).(data_cond)(c).tFine<=max(ar.model(m).(data_cond)(c).ppl.t(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx)),jx) = ...
            interp1(ar.model(m).(data_cond)(c).ppl.t(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
            ar.model(m).(data_cond)(c).ppl.(ppl_string)(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx),...
            ar.model(m).(data_cond)(c).tFine(ar.model(m).(data_cond)(c).tFine<=max(ar.model(m).(data_cond)(c).ppl.t(~isnan(ar.model(m).(data_cond)(c).ppl.(ppl_string)(:,jx)),jx))),...
            'pchip',NaN);    
  
    function [xFit_par] = getxFit(xFit)
        xFit_par = xFit;
        x1=ar.model(m).(data_cond)(c).ppl.x_orig(it_orig,jx);            
        xFit_par = xFit + (x1-xFit)/(ar.model(m).(data_cond)(c).tFine(it_orig)-t_tmp)*t_dir*stepsize;     
        if(isnan(xFit_par))
           ar.p=pReset;
           if(takeY)
                arLink(true,ar.model(m).(data_cond)(c).tExp(1),true,jx, c, m,NaN);
           end
            fprintf('ERROR IN STEP AT T=%d \n', t_tmp);
            return;
        end
    end   
end

    

function str = myNameTrafo(str)
str = strrep(str, '_', '\_');
end

function inv = chi2inv (x, n)
if (nargin ~= 2)
    error ('chi2inv: you must give two arguments');
end

if (~isscalar (n))
    [retval, x, n] = common_size(x, n);
    if (retval > 0)
        error ('chi2inv: x and n must be of common size or scalar');
    end
end

inv = gaminv(x, n / 2, 2);
end

function inv = gaminv (x, a, b)
if (nargin ~= 3)
    error ('gaminv: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gaminv: x, a and b must be of common size or scalars');
    end
end

sz = size (x);
inv = zeros (sz);

k = find ((x < 0) | (x > 1) | isnan (x) | ~(a > 0) | ~(b > 0));
if (any (k))
    inv (k) = NaN;
end

k = find ((x == 1) & (a > 0) & (b > 0));
if (any (k))
    inv (k) = Inf;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
    if (~isscalar(a) || ~isscalar(b))
        a = a (k);
        b = b (k);
        y = a .* b;
    else
        y = a * b * ones (size (k));
    end
    x = x (k);
    l = find (x < eps);
    if (any (l))
        y(l) = sqrt (eps) * ones (length (l), 1);
    end
    
    y_old = y;
    for i = 1 : 100
        
        h     = (gamcdf (y_old, a, b) - x) ./ gampdf (y_old, a, b);
        y_new = y_old - h;
        ind   = find (y_new <= eps);
        if (any (ind))
            y_new (ind) = y_old (ind) / 10;
            h = y_old - y_new;
        end
        if (max (abs (h)) < sqrt (eps))
            break;
        end
        y_old = y_new;
    end
    
    inv (k) = y_new;
end
end

function pdf = gampdf (x, a, b)
if (nargin ~= 3)
    error ('gampdf: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gampdf: x, a and b must be of common size or scalars');
    end
end

sz = size(x);
pdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    pdf (k) = NaN;
end

k = find ((x > 0) & (a > 0) & (a <= 1) & (b > 0));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = (x(k) .^ (a - 1)) ...
            .* exp(- x(k) ./ b) ./ gamma (a) ./ (b .^ a);
    else
        pdf(k) = (x(k) .^ (a(k) - 1)) ...
            .* exp(- x(k) ./ b(k)) ./ gamma (a(k)) ./ (b(k) .^ a(k));
    end
end

k = find ((x > 0) & (a > 1) & (b > 0));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = exp (- a .* log (b) + (a-1) .* log (x(k)) ...
            - x(k) ./ b - gammaln (a));
    else
        pdf(k) = exp (- a(k) .* log (b(k)) + (a(k)-1) .* log (x(k)) ...
            - x(k) ./ b(k) - gammaln (a(k)));
    end
end
end

function cdf = gamcdf (x, a, b)
if (nargin ~= 3)
    error ('gamcdf: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gamcdf: x, a and b must be of common size or scalars');
    end
end

sz = size (x);
cdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    cdf (k) = NaN;
end

k = find ((x > 0) & (a > 0) & (b > 0));
if (any (k))
    if (isscalar (a) && isscalar(b))
        cdf (k) = gammainc (x(k) ./ b, a);
    else
        cdf (k) = gammainc (x(k) ./ b(k), a(k));
    end
end
end
