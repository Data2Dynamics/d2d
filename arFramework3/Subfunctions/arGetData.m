function [t, y, ystd, tExp, yExp, yExpStd, lb, ub, ...
    yExpHl, dydt, y_ssa, y_ssa_lb, y_ssa_ub, qFit, t_ppl, y_ppl_ub, y_ppl_lb] = arGetData(jm, jd, jtype)

global ar

% SSA
if(jtype==1)
    if(isfield(ar.model(jm).data(jd),'yFineSSA'))
        y_ssa = ar.model(jm).data(jd).yFineSSA;
    else
        y_ssa = nan;
    end
    if(isfield(ar.model(jm).data(jd),'yFineSSA_lb'))
        y_ssa_lb = ar.model(jm).data(jd).yFineSSA_lb;
    else
        y_ssa_lb = nan;
    end
    if(isfield(ar.model(jm).data(jd),'yFineSSA_ub'))
        y_ssa_ub = ar.model(jm).data(jd).yFineSSA_ub;
    else
        y_ssa_ub = nan;
    end
    
elseif(jtype==2)
    % TODO implement x and z SSA
    y_ssa = nan;
    y_ssa_lb = nan;
    y_ssa_ub = nan;
    
elseif(jtype==3)
    y_ssa = nan;
    y_ssa_lb = nan;
    y_ssa_ub = nan;
    
end

if(isfield(ar.model(jm), 'data'))
    jc = ar.model(jm).data(jd).cLink;
else
    jc = 1;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end
t_ppl = [];
y_ppl_ub = [];
y_ppl_lb = [];

% trajectories and error bands
if(jtype==1 && isfield(ar.model(jm), 'data') && isfield(ar.model(jm).data(jd),'tFine'))
    t = ar.model(jm).data(jd).tFine;
    y = ar.model(jm).data(jd).yFineSimu;
    ystd = ar.model(jm).data(jd).ystdFineSimu;
    %Get data points of model profile likelihood
     if(isfield(ar.model(jm).data(jd),'ppl') && ~isempty(ar.model(jm).data(jd).ppl))
        t_ppl = ar.model(jm).data(jd).ppl.tstart;
        if(nansum(nansum(~isnan(ar.model(jm).data(jd).ppl.ub_fit))>0))
            y_ppl_ub = ar.model(jm).data(jd).ppl.ub_fit;
            y_ppl_lb = ar.model(jm).data(jd).ppl.lb_fit;
        else
            y_ppl_ub = ar.model(jm).data(jd).ppl.ub_fit_vpl;
            y_ppl_lb = ar.model(jm).data(jd).ppl.lb_fit_vpl;
        end
    end
elseif(jtype==2)
    t = ar.model(jm).condition(jc).tFine;
    y = [ar.model(jm).condition(jc).uFineSimu ar.model(jm).condition(jc).xFineSimu ...
        ar.model(jm).condition(jc).zFineSimu];
    ystd = [];
    %Get data points of model profile likelihood
    if(isfield(ar.model(jm).condition(jc),'ppl') && ~isempty(ar.model(jm).condition(jc).ppl))
        trows = size(ar.model(jm).condition(jc).ppl.tstart,1);
        t_ppl = [nan(trows,length(ar.model(1).qPlotU)) ar.model(jm).condition(jc).ppl.tstart nan(trows,length(ar.model(1).qPlotZ))];
        if(nansum(nansum(~isnan(ar.model(jm).condition(jc).ppl.ub_fit)))>0)
            y_ppl_ub = [nan(trows,length(ar.model(1).qPlotU)) ar.model(jm).condition(jc).ppl.ub_fit nan(trows,length(ar.model(1).qPlotZ))];
            y_ppl_lb = [nan(trows,length(ar.model(1).qPlotU)) ar.model(jm).condition(jc).ppl.lb_fit nan(trows,length(ar.model(1).qPlotZ))];
        else
            y_ppl_ub = [nan(trows,length(ar.model(1).qPlotU)) ar.model(jm).condition(jc).ppl.ub_fit_vpl nan(trows,length(ar.model(1).qPlotZ))];
            y_ppl_lb = [nan(trows,length(ar.model(1).qPlotU)) ar.model(jm).condition(jc).ppl.lb_fit_vpl nan(trows,length(ar.model(1).qPlotZ))];
        end
    end
elseif(jtype==3)
    t = ar.model(jm).condition(jc).tFine;
    y = ar.model(jm).condition(jc).vFineSimu;
    ystd = [];
    
else
    t = nan;
    y = nan;
    ystd = nan;
end

% data
if(jtype==1 && isfield(ar.model(jm), 'data') && isfield(ar.model(jm).data(jd), 'yExp') && ...
        ~isempty(ar.model(jm).data(jd).yExp))
    tExp = ar.model(jm).data(jd).tExp;
    yExp = ar.model(jm).data(jd).yExp;
    if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==-1) || ...
                        (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==-1) )
        yExpStd = ar.model(jm).data(jd).yExpStd;
    else
        if(isfield(ar.model(jm).data(jd),'ystdExpSimu'))
            yExpStd = ar.model(jm).data(jd).ystdExpSimu;
        else
            yExpStd = nan;
        end
    end
    if(isfield(ar.model(jm).data(jd),'highlight'))
        hl = ar.model(jm).data(jd).highlight;
    else
        hl = zeros(size(yExp));
    end
    yExpHl = yExp;
    yExpHl(hl==0) = NaN;
    qFit = ~isfield(ar.model(jm).data(jd),'qFit') | ...
        ar.model(jm).data(jd).qFit;
else
    tExp = [];
    yExp = [];
    yExpStd = [];
    yExpHl = [];
    qFit = [];
end

% confidence bands
lb = [];
ub = [];
if(jtype==1 && isfield(ar.model(jm), 'data') && isfield(ar.model(jm).data(jd),'tFine'))
    if(isfield(ar.model(jm).data(jd), 'yFineLB'))
        lb = ar.model(jm).data(jd).yFineLB;
        ub = ar.model(jm).data(jd).yFineUB;
    end
    
elseif(jtype==2)
    if(isfield(ar.model(jm).condition(jc), 'xFineLB'))
        lb = [ar.model(jm).condition(jc).xFineLB];
        ub = [ar.model(jm).condition(jc).xFineUB];
    end
    if(isfield(ar.model(jm).condition(jc), 'uFineLB'))
        lb = [ar.model(jm).condition(jc).uFineLB lb];
        ub = [ar.model(jm).condition(jc).uFineUB ub];
    end
    if(isfield(ar.model(jm).condition(jc), 'zFineLB'))
        lb = [lb ar.model(jm).condition(jc).zFineLB];
        ub = [ub ar.model(jm).condition(jc).zFineUB];
    end
elseif(jtype==3)
    if(isfield(ar.model(jm).condition(jc), 'vFineLB'))
        lb = ar.model(jm).condition(jc).vFineLB;
        ub = ar.model(jm).condition(jc).vFineUB;
    end
    
end

% steady state
if(jtype==2)
    dydt = ar.model(jm).condition(jc).dxdt;
    dydt = dydt(ar.model(jm).condition(jc).qSteadyState==1);
    if(~isempty(dydt))
        dydt = [nan(size(ar.model(jm).u)) dydt nan(size(ar.model(jm).z))];
    end
else
    dydt = [];
end
    
