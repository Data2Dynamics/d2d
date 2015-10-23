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
% takeY: if data or condition is used for integration of profile (default is true)
%
% gammas is vector of correction strengths for computation of integration
% steps for the vector of states ix (default is 0.1 for all ix)
%
% stepsize: Integration steps in time, e.g. 1, 0.5 ... (default is tFine stepsize)
%
% onlyProfile: Don't integrate prediction/validation confidence line over
% time, just compute profiles for specified t (default is false)
%
% tEnd: optional stopping time for the confidence band integration
% (default is tLim(end))
%
% rel_increase: increase of xTrial in computation of full profile for t's
% (default is 0.15)
%
% ed_steps: take integration steps or go parallel to orginial state solution (default is true)
%
% fineInt: do precise integration with a correction cycle at every
% integration step (default is false)
%
% backward: perform integration backwards in time (default is false)
%
% xstd: standard deviation of additional data point. If profile on data
% and ar.config.fiterrors=1, the fitted errors of the respective data is
% used (default is 0.1)
%
% doPPL: Alternative approach to integrate prediction confidence
% bands (default is false)
%
% trust_radius: If parameter change in integration step > boundaries,
% either whole step is downscaled or respective dimensions are downscaled
% (default is all downscaled)
%
% n_start: total steps in computation of full profile per side (default is 100)
%
% whichT: see t
%
% dir: 1 for only upper CI profile, -1 for lower, 0 for both (default is 0)
%
% constr_method: true means, constrained optimization is used for
% re-optimization of integration steps, false for consecutive parameter and
% data point optimization (more robust, default)
%
% Helge Hass, 2014 (helge.hass@fdm.uni-freiburg.de)

function doPPL(m, c, ix, t, takeY, gammas, stepsize, onlyProfile, tEnd, rel_increase, ed_steps, fineInt, backward, xstd, doPPL, trust_radius, n_start, whichT, dir, constr_method) % model, condition, states of interest, 

global ar
ar.config.SimuPPL=1;
ar.ppl.qLog10=0;

if(~exist('n_start','var') || isempty(n_start))
    n_start = 100;
end

ar.ppl.n=n_start;
n = n_start;

if(~exist('doPPL','var') || isempty(doPPL))
    doPPL = false;
end

if(~exist('onlyProfile','var') || isempty(onlyProfile))
    onlyProfile = false;
end

if(~exist('fineInt','var') || isempty(fineInt))
    fineInt = true;
end

if(~exist('backward','var') || isempty(backward))
    backward = false;
end

if(~exist('rel_increase','var') || isempty(rel_increase))
    rel_increase = 0.15;
end

if(~exist('trust_radius','var') || isempty(trust_radius))
    trust_radius = true;
end

ar.ppl.rel_increase=rel_increase;
if(~exist('xstd','var') || isempty(xstd))
    xstd = 0.1;
end

if(~exist('ed_steps','var') || isempty(ed_steps))
    ed_steps = true;
end

if(~exist('whichT','var') || isempty(whichT))
    whichT = 1;
end

if(~exist('tEnd','var') || isempty(tEnd))
    if(takeY)
       tEnd = ar.model(m).data(c).tFine(end);
    else
       tEnd = ar.model(m).condition(c).tFine(end); 
    end
end

if(~exist('dir','var') || isempty(dir))
    dir=0;
end

if(~exist('constr_method','var') || isempty(constr_method))
    constr_method = false;
end
ar.ppl.constr_method=constr_method;
if(~exist('takeY','var') || isempty(takeY))
    takeY = true;
end

if(~exist('stepsize','var') || isempty(stepsize))
    if(takeY)
        stepsize = 1/(size(ar.model(m).data(c).tFine,1)/(ar.model(m).data(c).tLim(end)-ar.model(m).data(c).tLim(1)));
    else
        stepsize = 1/(size(ar.model(m).condition(c).tFine,1)/(ar.model(m).condition(c).tFine(end)-ar.model(m).condition(c).tstart));
    end
end
nsteps = abs(floor((tEnd-t(whichT)) / stepsize));
ar.ppl.nsteps=nsteps;

if(~isempty(gammas) && length(gammas) ~= length(ix))
    gammas = repmat(gammas,1,length(ix)/length(gammas));
end

if(~exist('gammas','var') || isempty(gammas))
    gammas = ones(size(t))*1./stepsize;
end

% optimizer settings (set only once)
if(ar.config.fiterrors==1)
    ar.config.fiterrors=0;
    ar.ppl.fittederrors=1;
    ar.qFit(strncmp(ar.pLabel,'sd_',3))=2;
end
if(ar.config.useSensis)
    ar.config.optim.Jacobian = 'on';
else
    ar.config.optim.Jacobian = 'off';
end
ar.config.optim.Algorithm = 'trust-region-reflective';
%ar.config.optim.MaxIter = 200;

ar.ppl.dchi2 = chi2inv(1-ar.ppl.alpha_level, ar.ppl.ndof);
ar.ppl.dchi2;

arChi2();

ar.ppl.chi2_95 = nansum(ar.res.^2)+ar.ppl.dchi2 + 0.5;

% chi2 threshold
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
    ar.model(m).condition(c).xFineLB = nan(length(ar.model(m).condition(c).tFine),size(ar.model(m).condition(c).xExpSimu,2));
    ar.model(m).condition(c).xFineUB = nan(length(ar.model(m).condition(c).tFine),size(ar.model(m).condition(c).xExpSimu,2));

    ar.model(m).condition(c).ppl.t = nan(nsteps+1,size(ar.model(m).condition(c).xExpSimu,2));
    ar.model(m).condition(c).ppl.tstart = nan(length(t), size(ar.model(m).condition(c).xExpSimu,2));
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
        ar.model(m).data(c).yFineLB = nan(length(ar.model(m).data(c).tFine),size(ar.model(m).data(c).yExpSimu,2));
        ar.model(m).data(c).yFineUB = nan(length(ar.model(m).data(c).tFine),size(ar.model(m).data(c).yExpSimu,2));

        
        ar.model(m).data(c).ppl.t = nan(nsteps+1,size(ar.model(m).data(c).yExpSimu,2));
        ar.model(m).data(c).ppl.tstart = nan(length(t), size(ar.model(m).data(c).yExpSimu,2));
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

pReset = ar.p;
chi2start = nansum(ar.res.^2) + 0 ;

tic;

for jx = 1:length(ix)
    if(~takeY)
        %ar.ppl.qLog10=1;
        %ar.model(m).condition(c).ppl.x_orig(:,ix(jx))=log10(ar.model(m).condition(c).ppl.x_orig(:,ix(jx)));
    end
    %get starting value for jx
    if((~takeY && ((doPPL && isnan(ar.model(m).condition(c).ppl.kind_high(whichT,ix(jx)))) ...
            || (~doPPL && isnan(ar.model(m).condition(c).ppl.kind_high_vpl(whichT,ix(jx)))))) ...
            || (takeY && doPPL && isnan(ar.model(m).data(c).ppl.kind_high(whichT,ix(jx)))) ...
            || (takeY && ~doPPL && (isnan(ar.model(m).data(c).ppl.kind_high_vpl(whichT,ix(jx))))))
            
        if(ar.config.fiterrors~= -1 && takeY && ~isnan(ar.model(m).data(c).ystdFineSimu(1,ix(jx))))
            xstd = ar.model(m).data(c).ystdFineSimu(1,ix(jx));
           
        end
        xstart_ppl(m, c, ix(jx), t, doPPL, xstd, pReset, chi2start, whichT, takeY, true, 1, [], onlyProfile);
        if(onlyProfile)
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
    else
        if(~takeY)
            ar.model(m).condition(c).ppl.t(1,ix(jx)) = t(whichT);
        else
            ar.model(m).data(c).ppl.t(1,ix(jx)) = t(whichT);
        end
    end
    %loop through time  
    if((takeY && ((doPPL && (ar.model(m).data(c).ppl.lb_fit(whichT,ix(jx))==-Inf || ar.model(m).data(c).ppl.ub_fit(whichT,ix(jx))==Inf)) || ...
            (~doPPL && (ar.model(m).data(c).ppl.lb_fit_vpl(whichT,ix(jx))==-Inf || ar.model(m).data(c).ppl.ub_fit_vpl(whichT,ix(jx))==Inf))))...
            || (~takeY && ((doPPL && (ar.model(m).condition(c).ppl.lb_fit(whichT,ix(jx))==-Inf || ar.model(m).condition(c).ppl.ub_fit(whichT,ix(jx))==Inf)) ...
            || (~doPPL && (ar.model(m).condition(c).ppl.lb_fit_vpl(whichT,ix(jx))==-Inf || ar.model(m).condition(c).ppl.ub_fit_vpl(whichT,ix(jx))==Inf)))))
        fprintf('Starting points for Profile at t=%d not defined, check PPL computation!', t(whichT));
        if(takeY)
            arLink(true,0.,true,ix(jx), c, m,NaN);
        end
        ar.p = pReset;
        break;
    else
        if(takeY)        
            if(doPPL && ~isnan(ar.model(m).data(c).ppl.kind_high(whichT, ix(jx))))
                ar.model(m).data(c).ppl.x_high(1, ix(jx))=ar.model(m).data(c).ppl.ub_fit(whichT,ix(jx));
                ar.model(m).data(c).ppl.ps_high(1, ix(jx),:)=squeeze(ar.model(m).data(c).ppl.ps(whichT,ix(jx),ar.model(m).data(c).ppl.kind_high(whichT, ix(jx)),:));
            elseif(~doPPL  && ~isnan(ar.model(m).data(c).ppl.kind_high_vpl(whichT, ix(jx))))
                ar.model(m).data(c).ppl.x_high_vpl(1, ix(jx))=ar.model(m).data(c).ppl.ub_fit_vpl(whichT,ix(jx));
                ar.model(m).data(c).ppl.ps_high(1, ix(jx),:)=squeeze(ar.model(m).data(c).ppl.ps(whichT,ix(jx),ar.model(m).data(c).ppl.kind_high_vpl(whichT, ix(jx)),:));
            end
            
            if(doPPL && ~isnan(ar.model(m).data(c).ppl.kind_low(whichT, ix(jx))))
                ar.model(m).data(c).ppl.x_low(1, ix(jx))=ar.model(m).data(c).ppl.lb_fit(whichT, ix(jx));
                ar.model(m).data(c).ppl.ps_low(1, ix(jx),:)=squeeze(ar.model(m).data(c).ppl.ps(whichT, ix(jx),ar.model(m).data(c).ppl.kind_low(whichT, ix(jx)),:));
            elseif(~doPPL  && ~isnan(ar.model(m).data(c).ppl.kind_low_vpl(whichT, ix(jx))))
                ar.model(m).data(c).ppl.x_low_vpl(1, ix(jx))=ar.model(m).data(c).ppl.lb_fit_vpl(whichT, ix(jx));
                ar.model(m).data(c).ppl.ps_low(1, ix(jx),:)=squeeze(ar.model(m).data(c).ppl.ps(whichT, ix(jx),ar.model(m).data(c).ppl.kind_low_vpl(whichT, ix(jx)),:));
            end
            
            if(dir==1 && (doPPL && ~isnan(ar.model(m).data(c).ppl.kind_high(whichT, ix(jx)))))
                ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_high(1, ix(jx)), ar.model(m).data(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            elseif(dir==1 && ~doPPL && ~isnan(ar.model(m).data(c).ppl.kind_high_vpl(whichT, ix(jx))))
                ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_high_vpl(1, ix(jx)), ar.model(m).data(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            elseif(dir==-1 && (doPPL && ~isnan(ar.model(m).data(c).ppl.kind_low(whichT, ix(jx)))))
                ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_low(1, ix(jx)), ar.model(m).data(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            elseif(dir==-1 && (~doPPL && ~isnan(ar.model(m).data(c).ppl.kind_low_vpl(whichT, ix(jx)))))
                ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_low_vpl(1, ix(jx)), ar.model(m).data(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            else
                if(doPPL && ~isnan(ar.model(m).data(c).ppl.kind_high(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_high(1, ix(jx)), ar.model(m).data(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                elseif(~doPPL && ~isnan(ar.model(m).data(c).ppl.kind_high_vpl(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_high_vpl(1, ix(jx)), ar.model(m).data(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                end
                if(doPPL && ~isnan(ar.model(m).data(c).ppl.kind_low(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_low(1, ix(jx)), ar.model(m).data(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                elseif(~doPPL && ~isnan(ar.model(m).data(c).ppl.kind_low_vpl(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).data(c).ppl.x_low_vpl(1, ix(jx)), ar.model(m).data(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                end
            end
        else
            if(doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_high(whichT, ix(jx))))
                ar.model(m).condition(c).ppl.x_high(1, ix(jx))=ar.model(m).condition(c).ppl.ub_fit(whichT,ix(jx));
                ar.model(m).condition(c).ppl.ps_high(1, ix(jx),:)=squeeze(ar.model(m).condition(c).ppl.ps(whichT,ix(jx),ar.model(m).condition(c).ppl.kind_high(whichT, ix(jx)),:));
            elseif(~doPPL  && ~isnan(ar.model(m).condition(c).ppl.kind_high_vpl(whichT, ix(jx))))
                ar.model(m).condition(c).ppl.x_high_vpl(1, ix(jx))=ar.model(m).condition(c).ppl.ub_fit_vpl(whichT,ix(jx));
                ar.model(m).condition(c).ppl.ps_high(1, ix(jx),:)=squeeze(ar.model(m).condition(c).ppl.ps(whichT,ix(jx),ar.model(m).condition(c).ppl.kind_high_vpl(whichT, ix(jx)),:));
            end
            
            if(doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_low(whichT, ix(jx))))
                ar.model(m).condition(c).ppl.x_low(1, ix(jx))=ar.model(m).condition(c).ppl.lb_fit(whichT, ix(jx));
                ar.model(m).condition(c).ppl.ps_low(1, ix(jx),:)=squeeze(ar.model(m).condition(c).ppl.ps(whichT, ix(jx),ar.model(m).condition(c).ppl.kind_low(whichT, ix(jx)),:));
            elseif(~doPPL  && ~isnan(ar.model(m).condition(c).ppl.kind_low_vpl(whichT, ix(jx))))
                ar.model(m).condition(c).ppl.x_low_vpl(1, ix(jx))=ar.model(m).condition(c).ppl.lb_fit_vpl(whichT, ix(jx));
                ar.model(m).condition(c).ppl.ps_low(1, ix(jx),:)=squeeze(ar.model(m).condition(c).ppl.ps(whichT, ix(jx),ar.model(m).condition(c).ppl.kind_low_vpl(whichT, ix(jx)),:));
            end
            
            if(dir==1 && (doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_high(whichT, ix(jx)))))
                ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_high(1, ix(jx)), ar.model(m).condition(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            elseif(dir==1 && ~doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_high_vpl(whichT, ix(jx))))
                ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_high_vpl(1, ix(jx)), ar.model(m).condition(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            elseif(dir==-1 && (doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_low(whichT, ix(jx)))))
                ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_low(1, ix(jx)), ar.model(m).condition(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            elseif(dir==-1 && (~doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_low_vpl(whichT, ix(jx)))))
                ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_low_vpl(1, ix(jx)), ar.model(m).condition(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
            else
                if(doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_high(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_high(1, ix(jx)), ar.model(m).condition(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                elseif(~doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_high_vpl(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_high_vpl(1, ix(jx)), ar.model(m).condition(c).ppl.ps_high(1, ix(jx),:), t(whichT), doPPL, takeY, 1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                end
                if(doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_low(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_low(1, ix(jx)), ar.model(m).condition(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                elseif(~doPPL && ~isnan(ar.model(m).condition(c).ppl.kind_low_vpl(whichT, ix(jx))))
                    ppl_calc(m, c, ix(jx), ar.model(m).condition(c).ppl.x_low_vpl(1, ix(jx)), ar.model(m).condition(c).ppl.ps_low(1, ix(jx),:), t(whichT), doPPL, takeY, -1, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt);
                end
            end
        end
    end
    arWaitbar(-1);
    toc;
    if(takeY)
        arLink(true,0.,true,ix(jx), c, m,NaN);
    end
    ar.p = pReset;
        
end
    if(ar.config.fiterrors==0 && ar.ppl.fittederrors==1)
        ar.config.fiterrors=1;
        ar.qFit(strncmp(ar.pLabel,'sd_',3))=1;
    end
    ar.config.SimuPPL=0;
    arChi2(); 
    if(~takeY)
        ar.model(m).qPlotXs(ar.model(m).condition(c).dLink(1))=1;        
    else
        ar.model(m).qPlotYs(c) = 1;        
    end
    %ploterror_tmp = ar.config.ploterrors;
    %ar.config.ploterrors = -1;
    %arPlot2();
    %ar.config.ploterrors = ploterror_tmp;
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

arWaitbar(0);
ntot = ar.ppl.n;
tcount = 1;
for ts = 1:length(t)
    t_tmp = t(ts);
    fprintf('Caliculating PPL for t=%d and state x=%i \n',t(ts),jx);
    %if(save)
        if(takeY)
            if(ts == whichT && save)
                ar.model(m).data(c).ppl.t(1,jx) = t(ts);
            end
            if(save)
                ar.model(m).data(c).ppl.tstart(ts, jx) = t(ts);
            end
            [~,it_first] = min(abs(ar.model(m).data(c).tExp-t_tmp));            
            arLink(true, t_tmp, true, jx, c, m, ar.model(m).data(c).ppl.x_orig(it_first,jx), xstd);
            arChi2(true);
            [~,it] = min(abs(ar.model(m).data(c).tExp-t_tmp));
            if(length(find(ar.model(m).data(c).tExp==t_tmp))>1)
                it = it+1;               
            end
            xSim = ar.model(m).data(c).yExpSimu(it,jx);
            arLink(true, t_tmp, true, jx, c, m, xSim, xstd);
        else
            if(ts == whichT && save)
                ar.model(m).condition(c).ppl.t(1,jx) = t(ts);
            end
            if(save)
                ar.model(m).condition(c).ppl.tstart(ts, jx) = t(ts);
            end
            arLink(true, t_tmp);
            arChi2(true);
            [~,it] = min(abs(ar.model(m).condition(c).tExp-t_tmp));
            xSim = ar.model(m).condition(c).xExpSimu(it,jx);
            
        end
        if(ar.ppl.qLog10)
            xSim = log10(xSim);
        end
    % go up
    if(save || dir==1)
        npre = 0;
        [xtrial_up, xfit_up, ppl_up, vpl_up, ps_up, tcount] = ppl(m, c, jx, it, t_tmp, doPPL, xSim, ...
            chi2start, 1, ar.ppl.qLog10, ar.ppl.n, xstd, npre, ntot, tcount, takeY);
    
        % go down        
        ar.p = pReset;  
        arChi2(true);
        
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
            chi2start, -1, ar.ppl.qLog10, ar.ppl.n, xstd, npre, ntot, tcount, takeY);

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
    arChi2(true);
    xtrial_tmp = [fliplr(xtrial_down) xSim xtrial_up];
    xfit_tmp = [fliplr(xfit_down) xSim xfit_up];
    ppl_tmp = [fliplr(ppl_down) chi2start ppl_up];
    vpl_tmp = [fliplr(vpl_down) chi2start vpl_up];
    ps_tmp = [flipud(ps_down); pReset; ps_up];
    
    q_chi2good = ppl_tmp <= chi2start+ar.ppl.dchi2;
    q_nonnan = ~isnan(ppl_tmp); 
    
    vpl_chi2good = vpl_tmp <= chi2start+ar.ppl.dchi2;
    vpl_nonnan = ~isnan(vpl_tmp); 

    % calculate CI point-wise fit
    lb_tmp = min(xfit_tmp(q_chi2good));
    ub_tmp = max(xfit_tmp(q_chi2good));
    
    lb_tmp_vpl = min(xtrial_tmp(vpl_chi2good));
    ub_tmp_vpl = max(xtrial_tmp(vpl_chi2good));
    if(doPPL)
        if(lb_tmp==min(xfit_tmp(q_nonnan)))
            lb_tmp = -Inf;                
            fprintf('No -95 PPL for t=%d\n',t(ts));
        else
            kind_low_tmp = find(xfit_tmp==lb_tmp);
            if(length(kind_low_tmp)>1)
                kind_low_tmp = kind_low_tmp(1);
            end
            lb_tmp = interp1(ppl_tmp([kind_low_tmp kind_low_tmp-1]), ...
            xfit_tmp([kind_low_tmp kind_low_tmp-1]), chi2start+ar.ppl.dchi2);
        end
        if(ub_tmp==max(xfit_tmp(q_nonnan)))
            ub_tmp = Inf;
            fprintf('No +95 PPL for t=%d\n',t(ts));
        else
            kind_high_tmp = find(xfit_tmp==ub_tmp);   
            if(length(kind_high_tmp)>1)
                kind_high_tmp = kind_high_tmp(end);
            end
            ub_tmp = interp1(ppl_tmp([kind_high_tmp kind_high_tmp+1]), ...
            xfit_tmp([kind_high_tmp kind_high_tmp+1]), chi2start+ar.ppl.dchi2);
        end        
    else
        if(lb_tmp_vpl==min(xtrial_tmp(vpl_nonnan)))
            lb_tmp_vpl = -Inf;                
            fprintf('No -95 VPL for t=%d\n',t(ts));
            if(dir==-1 || dir==0)
%                 return;
            end
        else
            kind_low_tmp_vpl = find(xtrial_tmp==lb_tmp_vpl);
            if(length(kind_low_tmp_vpl)>1)
                kind_low_tmp_vpl = kind_low_tmp_vpl(1);
            end
            lb_tmp_vpl = interp1(vpl_tmp([kind_low_tmp_vpl kind_low_tmp_vpl-1]), ...
            xtrial_tmp([kind_low_tmp_vpl kind_low_tmp_vpl-1]), chi2start+ar.ppl.dchi2);
        end
        if(ub_tmp_vpl==max(xtrial_tmp(vpl_nonnan)))
            ub_tmp_vpl = Inf;
            fprintf('No +95 VPL for t=%d\n',t(ts));
            if(dir==1 || dir==0)
%                 return;
            end
        else
            kind_high_tmp_vpl = find(xtrial_tmp==ub_tmp_vpl); 
            if(length(kind_high_tmp_vpl)>1)
                kind_high_tmp_vpl = kind_high_tmp_vpl(end);
            end
            ub_tmp_vpl = interp1(vpl_tmp([kind_high_tmp_vpl kind_high_tmp_vpl+1]), ...
            xtrial_tmp([kind_high_tmp_vpl kind_high_tmp_vpl+1]), chi2start+ar.ppl.dchi2);
        end
    end
    if ~onlyProfile
        if(dir==1 && doPPL)        
            xFit = ub_tmp;
            ps = ps_tmp(kind_high_tmp,:);
        elseif(dir==-1 && doPPL)        
            xFit = lb_tmp;
            ps = ps_tmp(kind_low_tmp,:);        
        elseif(dir == 1 && ~doPPL)
            xFit = ub_tmp_vpl;
            ps = ps_tmp(kind_high_tmp_vpl,:);    
        else
            xFit = lb_tmp_vpl;
            ps = ps_tmp(kind_low_tmp_vpl,:);               
        end
    end
    if(takeY && save)   
        ar.model(m).data(c).ppl.xtrial(ts, jx,:) = xtrial_tmp;
        ar.model(m).data(c).ppl.xfit(ts, jx,:) = xfit_tmp;
        ar.model(m).data(c).ppl.ppl(ts, jx,:) = ppl_tmp;
        ar.model(m).data(c).ppl.vpl(ts, jx,:) = vpl_tmp;
        ar.model(m).data(c).ppl.ps(ts, jx,:,:) = ps_tmp;
        if ~onlyProfile
            if(doPPL)
                ar.model(m).data(c).ppl.lb_fit(ts, jx) = lb_tmp;
                ar.model(m).data(c).ppl.ub_fit(ts, jx) = ub_tmp;
                ar.model(m).data(c).ppl.kind_high(ts, jx) = kind_high_tmp;
                ar.model(m).data(c).ppl.kind_low(ts, jx) = kind_low_tmp;
            else
                ar.model(m).data(c).ppl.lb_fit_vpl(ts, jx) = lb_tmp_vpl;
                ar.model(m).data(c).ppl.ub_fit_vpl(ts, jx) = ub_tmp_vpl;
                ar.model(m).data(c).ppl.kind_high_vpl(ts, jx) = kind_high_tmp_vpl;
                ar.model(m).data(c).ppl.kind_low_vpl(ts, jx) = kind_low_tmp_vpl;
            end
        end
    elseif(~takeY && save)
        ar.model(m).condition(c).ppl.xtrial(ts, jx,:) = xtrial_tmp;
        ar.model(m).condition(c).ppl.xfit(ts, jx,:) = xfit_tmp;
        ar.model(m).condition(c).ppl.ppl(ts, jx,:) = ppl_tmp;
        ar.model(m).condition(c).ppl.vpl(ts, jx,:) = vpl_tmp;
        ar.model(m).condition(c).ppl.ps(ts, jx,:,:) = ps_tmp;
        if ~onlyProfile
            if(doPPL)            
                ar.model(m).condition(c).ppl.lb_fit(ts, jx) = lb_tmp;
                ar.model(m).condition(c).ppl.ub_fit(ts, jx) = ub_tmp;
                ar.model(m).condition(c).ppl.kind_high(ts, jx) = kind_high_tmp;
                ar.model(m).condition(c).ppl.kind_low(ts, jx) = kind_low_tmp;
            else
                ar.model(m).condition(c).ppl.lb_fit_vpl(ts, jx) = lb_tmp_vpl;
                ar.model(m).condition(c).ppl.ub_fit_vpl(ts, jx) = ub_tmp_vpl;
                ar.model(m).condition(c).ppl.kind_high_vpl(ts, jx) = kind_high_tmp_vpl;
                ar.model(m).condition(c).ppl.kind_low_vpl(ts, jx) = kind_low_tmp_vpl;  
            end
        end
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

dx = sqrt(ar.ppl.dchi2*ar.ppl.rel_increase) * xstd;

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
    if(takeY)
        arLink(true,t,true,ix, c,m,xExp,xstd);
    else
        arLink(true,t);
    end
    xtrial(j) = xExp;    
    arChi2(ar.config.useSensis, ar.p(ar.qFit==1))
    if(takeY)
        xSim = ar.model(m).data(c).yExpSimu(it,ix); 
    else
        xSim = ar.model(m).condition(c).xExpSimu(it,ix); 
    end
    xfit(j) = xSim;
    if(takeY)
        ppl(j) = nansum(ar.res.^2) - (xtrial(j)-xfit(j)).^2/xstd.^2;
        vpl(j) = nansum(ar.res.^2);
        
    else
        ppl(j) = nansum(ar.res.^2);
        vpl(j) = nansum(ar.res.^2) + (xtrial(j)-xfit(j)).^2/xstd.^2;
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
        if(takeY)
            arLink(true,t,true,ix, c,m,xExp,xstd);
        else
            arLink(true,t);
        end
        arChi2(ar.config.useSensis, pTrial)
        
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

function ppl_calc(m, c, jx, xFit, p, t, doPPL, takeY, dir, stepsize, xstd, ed_steps, constr_method, pReset, chi2start, trust_radius, backward, fineInt)
    global ar
    ar.ppl.xFit_tmp = xFit;
    chi2 = chi2start + ar.ppl.dchi2;
    xSim2 = NaN;
    xSim = NaN;
    xSim3 = NaN;
    npre=0;
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
        
        if(takeY)
            [~,it_orig] = min(abs(ar.model(m).data(c).tFine-t_tmp-t_dir*stepsize));
            it_getfRHS = it_orig;
            x1_orig=ar.model(m).data(c).ppl.x_orig(it_orig,jx);

            if(ar.model(m).data(c).tFine(it_orig)-t_tmp-t_dir*stepsize > 0)
                x2_orig=x1_orig;   
                if(it_orig>1);
                    x1_orig=ar.model(m).data(c).ppl.x_orig(it_orig-1,jx);
                end
            else
                it_orig=it_orig+1;
                if(size(ar.model(m).data(c).ppl.x_orig,1)>=it_orig);
                    x2_orig=ar.model(m).data(c).ppl.x_orig(it_orig,jx);
                else
                    x2_orig = x1_orig;
                    it_orig=it_orig-1;
                end
            end
            if(it_orig>1)
                x_orig = x1_orig + (x2_orig-x1_orig)/(ar.model(m).data(c).tFine(it_orig)-ar.model(m).data(c).tFine(it_orig-1))*(t_tmp+t_dir*stepsize-ar.model(m).data(c).tFine(it_orig-1));
            else
               x_orig = x1_orig; 
            end
        else
            [~,it_orig] = min(abs(ar.model(m).condition(c).tFine-t_tmp-t_dir*stepsize));
            it_getfRHS = it_orig;
            x1_orig=ar.model(m).condition(c).ppl.x_orig(it_orig,jx);

            if(ar.model(m).condition(c).tFine(it_orig)-t_tmp-t_dir*stepsize > 0)
                x2_orig=x1_orig;    
                if(it_orig>1);
                    x1_orig=ar.model(m).condition(c).ppl.x_orig(it_orig-1,jx);
                end
            else
                              
                if(size(ar.model(m).condition(c).ppl.x_orig,1)>it_orig)
                    it_orig=it_orig+1;
                    x2_orig=ar.model(m).condition(c).ppl.x_orig(it_orig,jx);
                else                    
                    x2_orig = x1_orig;
                end
            end
            if(it_orig>1)
                x_orig = x1_orig + (x2_orig-x1_orig)/(ar.model(m).condition(c).tFine(it_orig)-ar.model(m).condition(c).tFine(it_orig-1))*(t_tmp+t_dir*stepsize-ar.model(m).condition(c).tFine(it_orig-1));
            else
               x_orig = x1_orig; 
            end
        end
        if(toc>tcount)        
            if(dir==1)
                arWaitbar((jt+npre), nsteps, sprintf('PPL (up) for %s at t=%g %i/%i', xLabel, t_tmp, jt, nsteps));
            else
                arWaitbar((jt+npre), nsteps, sprintf('PPL (down) for %s at t=%g %i/%i', xLabel, t_tmp, jt, nsteps));
            end
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
                if(trust_radius)
                    region_fac_p = min(abs((ar.p(ar.qFit==1) - ar.ub(ar.qFit==1))./ (2*dps'*stepsize)));
                else
                    region_fac_p = abs( (ar.p(ar.qFit==1) - ar.ub(ar.qFit==1)) ./ (2*dps'*stepsize));
                    region_fac_p(region_fac_p>1) = 1;
                end
                    still_bad = abs((ar.p(ar.qFit==1) - ar.lb(ar.qFit==1)) ./ (2*dps'.*region_fac_p*stepsize)) < region_fac_p;
               if(sum( still_bad > 0))
                    region_low = abs( (ar.p(ar.qFit==1) - ar.lb(ar.qFit==1)) ./ (2*dps'*stepsize));
                    if(trust_radius)
                        region_fac_p = min(region_low);
                    else
                        region_fac_p(still_bad) = region_low(still_bad);
                    end
               end
            end
            
            if(max(region_fac_p) < 0.1)
                region_fac_p = zeros(1,length(dps));
            end     
            ar.p(ar.qFit==1)=ar.p(ar.qFit==1) + dps'*stepsize.*region_fac_p;
            ar.ppl.xFit_tmp = ar.ppl.xFit_tmp + dx(1)*stepsize*region_fac;
            ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
            ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);
%             end
            t_tmp = t_tmp + t_dir*stepsize;  
            if(doPPL)
                [chi2, xSim] = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, doPPL, t_dir*stepsize, ar.ppl.xFit_tmp, xstd);
                ar.ppl.xFit_tmp = xSim;
            end
            curve_xSim = 0;
            if(doPPL && jt>2 && it_orig<length(ar.model(m).condition(c).tFine))   
                %curve_xSim = (ar.model(m).condition(c).ppl.x_orig(it_orig+1,jx) - 2*ar.model(m).condition(c).ppl.x_orig(it_orig,jx) ...
                %            + ar.model(m).condition(c).ppl.x_orig(it_orig-1,jx)) / (ar.model(m).condition(c).tFine(it_orig) - ar.model(m).condition(c).tFine(it_orig-1))^2;
                               
            end
            last_corr=0;
            if(takeY && jt>10)
                if(dir==1)
                    last_corr = sum(ar.model(m).data(c).ppl.corr_high(jt-10:jt,jx));
                else
                    last_corr = sum(ar.model(m).data(c).ppl.corr_low(jt-10:jt,jx));
                end
            elseif(jt>10)
                if(dir==1)
                    last_corr = sum(ar.model(m).condition(c).ppl.corr_high(jt-10:jt,jx));
                else
                    last_corr = sum(ar.model(m).condition(c).ppl.corr_low(jt-10:jt,jx));
                end
            end
            if(((fineInt || (max(abs(dps))<1.e-4 || ~isempty(find(~region_fac_p,1))) ...
                    || (ar.ppl.xFit_tmp<x_orig && dir==1) || (ar.ppl.xFit_tmp>x_orig && dir==-1)) ...
                    && last_corr==0) || mod(jt/(floor(nsteps/10)),1)==0)%  || max(abs(dps))==0 || (doPPL && ((curve_xSim>0 && dir ==1) || (curve_xSim<0 && dir==-1))))
                %chi2 = PPL_doSim_calc(chi2, x_orig, m, c, jt, jx, xstd, t_tmp, qLog10, dir, takeY, doPPL);
                %[chi2, xSim] = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, doPPL, stepsize, ar.ppl.xFit_tmp, xstd);
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
%                     if(dir==1)
%                         ar.ppl.xFit_tmp =  ar.ppl.xFit_tmp + 3*xstd;
%                     else
%                         ar.ppl.xFit_tmp =  ar.ppl.xFit_tmp - 3*xstd;
%                     end
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
                arLink(true,0.,true,jx, c, m,NaN);
                [ar.ppl.xFit_tmp, ar.p] = xstart_ppl(m, c, jx, t_tmp, doPPL, xstd, pReset, chi2start, 10, takeY, false, dir, ar.ppl.xFit_tmp);                
                chi2 = PPL_chi2(t_tmp,false, m, c, jx, takeY, qLog10, doPPL, stepsize, ar.ppl.xFit_tmp, xstd);               
                if((chi2 - ar.ppl.chi2_95 + 0.5) > 0.2 || (chi2 - ar.ppl.chi2_95 + 0.5) < -0.2 )
                    fprintf('Check the Profile at t=%d for inconsistencies, diff in chi2 is %d \n',t_tmp,abs(chi2 - ar.ppl.chi2_95 + 0.5));
                end
            end
        end
        
        
%         %calculate gamma change
%         if(jt>2 && ~doPPL)
%             if(dir==1)
%                 if(takeY)                             
%                    %chi2_change = chi2 - ar.ppl.chi2_tmp-(grad*squeeze(ar.model(m).data(c).ppl.ps_high(jt,jx,ar.qFit==1)));
%                    ps_prod = (ar.p-squeeze(ar.model(m).data(c).ppl.ps_high(jt,jx,:))')/norm(ar.p-squeeze(ar.model(m).data(c).ppl.ps_high(jt,jx,:))') * ...
%                        (squeeze(ar.model(m).data(c).ppl.ps_high(jt,jx,:))-squeeze(ar.model(m).data(c).ppl.ps_high(jt-1,jx,:)))/norm(squeeze(ar.model(m).data(c).ppl.ps_high(jt,jx,:))-squeeze(ar.model(m).data(c).ppl.ps_high(jt-1,jx,:)));
%                 else
%                    %chi2_change = chi2 - ar.ppl.chi2_tmp-(grad*squeeze(ar.model(m).condition(c).ppl.ps_high(jt,jx,ar.qFit==1)));
%                    ps_prod = (ar.p-squeeze(ar.model(m).condition(c).ppl.ps_high(jt,jx,:))')/norm(ar.p-squeeze(ar.model(m).condition(c).ppl.ps_high(jt,jx,:))') * ...
%                        (squeeze(ar.model(m).condition(c).ppl.ps_high(jt,jx,:))-squeeze(ar.model(m).condition(c).ppl.ps_high(jt-1,jx,:)))/norm(squeeze(ar.model(m).condition(c).ppl.ps_high(jt,jx,:))-squeeze(ar.model(m).condition(c).ppl.ps_high(jt-1,jx,:)));
%                 end
%             else
%                 if(takeY)                             
%                    %chi2_change = chi2 - ar.ppl.chi2_tmp-(grad*squeeze(ar.model(m).data(c).ppl.ps_low(jt,jx,ar.qFit==1)));
%                    ps_prod = (ar.p-squeeze(ar.model(m).data(c).ppl.ps_low(jt,jx,:))')/norm(ar.p-squeeze(ar.model(m).data(c).ppl.ps_low(jt,jx,:))') * ...
%                        (squeeze(ar.model(m).data(c).ppl.ps_low(jt,jx,:))-squeeze(ar.model(m).data(c).ppl.ps_low(jt-1,jx,:)))/norm(squeeze(ar.model(m).data(c).ppl.ps_low(jt,jx,:))-squeeze(ar.model(m).data(c).ppl.ps_low(jt-1,jx,:)));
%                 else
%                    %chi2_change = chi2 - ar.ppl.chi2_tmp-(grad*squeeze(ar.model(m).condition(c).ppl.ps_low(jt,jx,ar.qFit==1)));
%                    ps_prod = (ar.p-squeeze(ar.model(m).condition(c).ppl.ps_low(jt,jx,:))')/norm(ar.p-squeeze(ar.model(m).condition(c).ppl.ps_low(jt,jx,:))') * ...
%                        (squeeze(ar.model(m).condition(c).ppl.ps_low(jt,jx,:))-squeeze(ar.model(m).condition(c).ppl.ps_low(jt-1,jx,:)))/norm(squeeze(ar.model(m).condition(c).ppl.ps_low(jt,jx,:))-squeeze(ar.model(m).condition(c).ppl.ps_low(jt-1,jx,:)));
%                 end
%             end
%             if((ps_prod>-0.9 && ps_prod<-0.1) || (ps_prod>0.1 && ps_prod<0.9))
%                  gamma_tmp = 1/2*gamma_tmp;                 
%             end
%         end
        
        ar.ppl.chi2_tmp = chi2;
        if(dir==1)
            if(takeY)
                ar.model(m).data(c).ppl.corr_high(jt+t_dir*1, jx)=corr_tmp;
                if(doPPL)
                    ar.model(m).data(c).ppl.x_high(jt+t_dir*1, jx)=xSim;
                    ar.model(m).data(c).ppl.ppl_high(jt+t_dir*1, jx)=chi2;
                else
                    ar.model(m).data(c).ppl.x_high_vpl(jt+t_dir*1, jx)=ar.ppl.xFit_tmp;
                    ar.model(m).data(c).ppl.vpl_high(jt+t_dir*1, jx)=chi2;
                end
                ar.model(m).data(c).ppl.t(jt+t_dir*1,jx)=t_tmp;
                ar.model(m).data(c).ppl.ps_high(jt+t_dir*1, jx,:)=ar.p;
            else
                ar.model(m).condition(c).ppl.corr_high(jt+t_dir*1, jx)=corr_tmp;
                if(doPPL)
                    ar.model(m).condition(c).ppl.x_high(jt+t_dir*1, jx)=xSim;
                    ar.model(m).condition(c).ppl.ppl_high(jt+t_dir*1, jx)=chi2;
                else
                    ar.model(m).condition(c).ppl.x_high_vpl(jt+t_dir*1, jx)=ar.ppl.xFit_tmp;
                    ar.model(m).condition(c).ppl.vpl_high(jt+t_dir*1, jx)=chi2;
                end
                ar.model(m).condition(c).ppl.t(jt+t_dir*1,jx)=t_tmp;
                ar.model(m).condition(c).ppl.ps_high(jt+t_dir*1, jx,:)=ar.p;
            end
        elseif(dir==-1)
            
            if(takeY)
                ar.model(m).data(c).ppl.corr_low(jt+t_dir*1, jx)=corr_tmp;
                if(doPPL)
                    ar.model(m).data(c).ppl.x_low(jt+t_dir*1, jx)=xSim;
                    ar.model(m).data(c).ppl.ppl_low(jt+t_dir*1, jx)=chi2;
                else
                    ar.model(m).data(c).ppl.x_low_vpl(jt+t_dir*1, jx)=ar.ppl.xFit_tmp;
                    ar.model(m).data(c).ppl.vpl_low(jt+t_dir*1, jx)=chi2;
                end
                ar.model(m).data(c).ppl.t(jt+t_dir*1,jx)=t_tmp;
                ar.model(m).data(c).ppl.ps_low(jt+t_dir*1, jx,:)=ar.p;
            else
                ar.model(m).condition(c).ppl.corr_low(jt+t_dir*1, jx)=corr_tmp;
                if(doPPL)
                    ar.model(m).condition(c).ppl.x_low(jt+t_dir*1, jx)=xSim;
                    ar.model(m).condition(c).ppl.ppl_low(jt+t_dir*1, jx)=chi2;
                else
                    ar.model(m).condition(c).ppl.x_low_vpl(jt+t_dir*1, jx)=ar.ppl.xFit_tmp;
                    ar.model(m).condition(c).ppl.vpl_low(jt+t_dir*1, jx)=chi2;
                end
                ar.model(m).condition(c).ppl.t(jt+t_dir*1,jx)=t_tmp;
                ar.model(m).condition(c).ppl.ps_low(jt+t_dir*1, jx,:)=ar.p;
            end
            
        end       
        
%         if(t_dir==-1 && jt>2 && ((doPPL && ((dir==1 && ((ar.model(m).data(c).ppl.x_high(jt, jx)-ar.model(m).data(c).ppl.x_high(jt-1, jx))/(ar.model(m).data(c).ppl.x_high(jt-1, jx)-ar.model(m).data(c).ppl.x_high(jt-2, jx))>1 && ...
%                 (ar.model(m).data(c).ppl.x_high(jt, jx)-ar.model(m).data(c).ppl.x_high(jt-1, jx))/(ar.model(m).data(c).ppl.x_high(jt-1, jx)-ar.model(m).data(c).ppl.x_high(jt-2, jx))<-0.5)) ...
%                 || (dir==-1 && ((ar.model(m).data(c).ppl.x_low(jt, jx)-ar.model(m).data(c).ppl.x_low(jt-1, jx))/(ar.model(m).data(c).ppl.x_low(jt-1, jx)-ar.model(m).data(c).ppl.x_low(jt-3, jx))>1 && ...
%                 (ar.model(m).data(c).ppl.x_low(jt, jx)-ar.model(m).data(c).ppl.x_low(jt-1, jx))/(ar.model(m).data(c).ppl.x_low(jt-1, jx)-ar.model(m).data(c).ppl.x_low(jt-3, jx))<-0.5 ))))...
%                 || (~doPPL && ((dir==1 && ((ar.model(m).data(c).ppl.x_high_vpl(jt, jx)-ar.model(m).data(c).ppl.x_high_vpl(jt-1, jx))/(ar.model(m).data(c).ppl.x_high_vpl(jt-1, jx)-ar.model(m).data(c).ppl.x_high_vpl(jt-2, jx))>1 && ...
%                 (ar.model(m).data(c).ppl.x_high_vpl(jt, jx)-ar.model(m).data(c).ppl.x_high_vpl(jt-1, jx))/(ar.model(m).data(c).ppl.x_high_vpl(jt-1, jx)-ar.model(m).data(c).ppl.x_high_vpl(jt-2, jx))<-0.5)) ...
%                 || (dir==-1 && ((ar.model(m).data(c).ppl.x_low_vpl(jt, jx)-ar.model(m).data(c).ppl.x_low_vpl(jt-1, jx))/(ar.model(m).data(c).ppl.x_low_vpl(jt-1, jx)-ar.model(m).data(c).ppl.x_low_vpl(jt-2, jx))>1 && ...
%                 (ar.model(m).data(c).ppl.x_low_vpl(jt, jx)-ar.model(m).data(c).ppl.x_low_vpl(jt-1, jx))/(ar.model(m).data(c).ppl.x_low_vpl(jt-1, jx)-ar.model(m).data(c).ppl.x_low_vpl(jt-2, jx))<-0.5))))))
%            t_dir = 1;
%            jt = ar.ppl.bkp_jt;
%            t_tmp = t+jt*stepsize;
%            fprintf('stopping backward integration, onwards from t=%f \n',t+jt*stepsize);
%            continue
%         end
%         
        if(t_dir==-1 && jt>1)
            jt=jt-2;
        elseif(t_dir==-1 && jt<=1)
            fprintf('Backward integration stopped because lower time bound hit, proceeding with normal integration \n');
            break
%             t_dir = 1;  
%             jt = ar.ppl.bkp_jt;
%             t_tmp = t+jt*stepsize;
        end
%         %backward integration if step in integration band
%         if(t_dir==1 && jt>2 && ((doPPL && ((dir==1 && ((xSim-ar.model(m).data(c).ppl.x_high(jt, jx))/(ar.model(m).data(c).ppl.x_high(jt, jx)-ar.model(m).data(c).ppl.x_high(jt-1, jx))>3 || ...
%                 (xSim-ar.model(m).data(c).ppl.x_high(jt, jx))/(ar.model(m).data(c).ppl.x_high(jt, jx)-ar.model(m).data(c).ppl.x_high(jt-1, jx))<-0.5)) ...
%                 || (dir==-1 && ((xSim-ar.model(m).data(c).ppl.x_low(jt, jx))/(ar.model(m).data(c).ppl.x_low(jt, jx)-ar.model(m).data(c).ppl.x_low(jt-1, jx))>3 || ...
%                 (xSim-ar.model(m).data(c).ppl.x_low(jt, jx))/(ar.model(m).data(c).ppl.x_low(jt, jx)-ar.model(m).data(c).ppl.x_low(jt-1, jx))<-0.5 ))))...
%                 || (~doPPL && ((dir==1 && ((ar.ppl.xFit_tmp-ar.model(m).data(c).ppl.x_high_vpl(jt, jx))/(ar.model(m).data(c).ppl.x_high_vpl(jt, jx)-ar.model(m).data(c).ppl.x_high_vpl(jt-1, jx))>3 || ...
%                 (ar.ppl.xFit_tmp-ar.model(m).data(c).ppl.x_high_vpl(jt, jx))/(ar.model(m).data(c).ppl.x_high_vpl(jt, jx)-ar.model(m).data(c).ppl.x_high_vpl(jt-1, jx))<-0.5)) ...
%                 || (dir==-1 && ((ar.ppl.xFit_tmp-ar.model(m).data(c).ppl.x_low_vpl(jt, jx))/(ar.model(m).data(c).ppl.x_low_vpl(jt, jx)-ar.model(m).data(c).ppl.x_low_vpl(jt-1, jx))>3 || ...
%                 (ar.ppl.xFit_tmp-ar.model(m).data(c).ppl.x_low_vpl(jt, jx))/(ar.model(m).data(c).ppl.x_low_vpl(jt, jx)-ar.model(m).data(c).ppl.x_low_vpl(jt-1, jx))<-0.5))))))
%            fprintf('discovered jump at t=%f , going backwards \n',t_tmp);
%             t_dir = -1;
%            ar.ppl.bkp_jt = jt;
%         end
       
    end
    %write LB/UB in ar struct
    if(takeY)
        if(doPPL)
            if(dir==1)
                ar.model(m).data(c).yFineUB(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_high(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_high(:,jx)),jx),...
                    ar.model(m).data(c).ppl.x_high(~isnan(ar.model(m).data(c).ppl.x_high(:,jx)),jx),...
                    ar.model(m).data(c).tFine(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_high(:,jx)),jx))),...
                    'pchip',NaN);
            elseif(dir==-1)
                ar.model(m).data(c).yFineLB(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_low(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_low(:,jx)),jx),...
                    ar.model(m).data(c).ppl.x_low(~isnan(ar.model(m).data(c).ppl.x_low(:,jx)),jx),...
                    ar.model(m).data(c).tFine(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_low(:,jx)),jx))),...
                    'pchip',NaN);
            end
        else
            if(dir==1)
                ar.model(m).data(c).yFineUB(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_high_vpl(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_high_vpl(:,jx)),jx),...
                    ar.model(m).data(c).ppl.x_high_vpl(~isnan(ar.model(m).data(c).ppl.x_high_vpl(:,jx)),jx),...
                    ar.model(m).data(c).tFine(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_high_vpl(:,jx)),jx))),...
                    'pchip',NaN);
            elseif(dir==-1)
                ar.model(m).data(c).yFineLB(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_low_vpl(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_low_vpl(:,jx)),jx),...
                    ar.model(m).data(c).ppl.x_low_vpl(~isnan(ar.model(m).data(c).ppl.x_low_vpl(:,jx)),jx),...
                    ar.model(m).data(c).tFine(ar.model(m).data(c).tFine<=max(ar.model(m).data(c).ppl.t(~isnan(ar.model(m).data(c).ppl.x_low_vpl(:,jx)),jx))),...
                    'pchip',NaN);
            end
        end                                
    else
        if(doPPL)
            if(dir==1)
                ar.model(m).condition(c).xFineUB(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_high(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_high(:,jx)),jx),...
                    ar.model(m).condition(c).ppl.x_high(~isnan(ar.model(m).condition(c).ppl.x_high(:,jx)),jx),...
                    ar.model(m).condition(c).tFine(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_high(:,jx)),jx))),...
                    'pchip',NaN);
            elseif(dir==-1)
                ar.model(m).condition(c).xFineLB(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_low(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_low(:,jx)),jx),...
                    ar.model(m).condition(c).ppl.x_low(~isnan(ar.model(m).condition(c).ppl.x_low(:,jx)),jx),...
                    ar.model(m).condition(c).tFine(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_low(:,jx)),jx))),...
                    'pchip',NaN);
            end
            else
            if(dir==1)
                ar.model(m).condition(c).xFineUB(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_high_vpl(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_high_vpl(:,jx)),jx),...
                    ar.model(m).condition(c).ppl.x_high_vpl(~isnan(ar.model(m).condition(c).ppl.x_high_vpl(:,jx)),jx),...
                    ar.model(m).condition(c).tFine(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_high_vpl(:,jx)),jx))),...
                    'pchip',NaN);
            elseif(dir==-1)
                ar.model(m).condition(c).xFineLB(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_low_vpl(:,jx)),jx)),jx) = ...
                    interp1(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_low_vpl(:,jx)),jx),...
                    ar.model(m).condition(c).ppl.x_low_vpl(~isnan(ar.model(m).condition(c).ppl.x_low_vpl(:,jx)),jx),...
                    ar.model(m).condition(c).tFine(ar.model(m).condition(c).tFine<=max(ar.model(m).condition(c).ppl.t(~isnan(ar.model(m).condition(c).ppl.x_low_vpl(:,jx)),jx))),...
                    'pchip',NaN);
            end
        end 
    end
    
    
    function [xFit_par] = getxFit(xFit)
        xFit_par = xFit;
        if(takeY)
            x1=ar.model(m).data(c).ppl.x_orig(it_orig,jx);            
            xFit_par = xFit + (x1-xFit)/(ar.model(m).data(c).tFine(it_orig)-t_tmp)*t_dir*stepsize;                
        else
            x1=ar.model(m).condition(c).ppl.x_orig(it_orig,jx);            
            xFit_par = xFit + (x1-xFit)/(ar.model(m).condition(c).tFine(it_orig)-t_tmp)*t_dir*stepsize;
        end
        if(isnan(xFit_par))
           ar.p=pReset;
           if(takeY)
                arLink(true,ar.model(m).data(c).tExp(1),true,jx, c, m,NaN);
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
