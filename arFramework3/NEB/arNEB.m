function [spring_out, q_pathfound] = arNEB(i_pstart, i_pend, springs, qplot)
% performs NEB method fitting between 'i_pstart' and 'i_pend' (fit indicies)
% for a vector of spring constants 'springs' until a connecting path is found.
% Plots if flag 'qplot' is active.

global ar

% initialization
arInitNEBPath
ar.merger.neb.pstart = ar.mergerlhs.ps(i_pstart,:);
ar.merger.neb.pend = ar.mergerlhs.ps(i_pend,:);

spring_out = nan;

% caluclate inital path
r =  ar.merger.neb.pend - ar.merger.neb.pstart;
r_norm =  norm(r);
ar.merger.neb.ps_init = [];
ar.merger.neb.ps_init(1,:) = ar.merger.neb.pstart;
for i = 1:(ar.merger.neb.steps+1)
    ar.merger.neb.ps_init(i+1,:) = r_norm / (ar.merger.neb.steps+1) * i * r./r_norm + ar.merger.neb.pstart;
end


ar.merger.neb.spr = [];
ar.merger.neb.spring_out = [];

q_intipathfound = arNEBCheckInitPath;

%check again with very strict integration tolerances
if q_intipathfound == 0
    ar.bkptols.atol = ar.config.atol;
    ar.bkptols.rtol = ar.config.rtol;
    ar.bkptols.msteps = ar.config.maxsteps;
    
    ar.config.atol = 1e-12;
    ar.config.rtol = 1e-12;
    ar.config.maxsteps = 1e6;

    q_intipathfound = arNEBCheckInitPath;

    ar.config.atol = ar.bkptols.atol;
    ar.config.rtol = ar.bkptols.rtol;
    ar.config.maxsteps = ar.bkptols.msteps;
end


% NEB
for ispr =1:length(springs)
    
    % break if initial path is already suffcient
    if ispr == 1 && q_intipathfound == 1
        %for function output
        spring_out = 0;
        q_pathfound = q_intipathfound;
        break
    end
    
    ar.merger.neb.spr(ispr).springconst = springs(ispr);
    
    fprintf(['fitting NEB with ' num2str(ar.merger.neb.steps)...
        ' steps and spring constant ' ...
        num2str(ar.merger.neb.spr(ispr).springconst) '\n'] )
    
    tic 
    
    arNEBFit(springs(ispr)) % do the actual fitting  / relaxing fo the band
    fprintf('\n')
    
    ar.merger.neb.spr(ispr).time = toc;
    
    % save & store current spring results
    ar.merger.neb.spr(ispr).ps_result = ar.merger.neb.ps_result;
    ar.merger.neb.spr(ispr).chi2s = ar.merger.neb.chi2_step;
    
    ar.merger.neb.spr(ispr).res = ...
        [ar.merger.neb.res_start, ...
        ar.merger.neb.res_spring(3:end), ...
        ar.merger.neb.res_end];
    
    ar.merger.neb.spr(ispr).sres = ...
        [ar.merger.neb.sres_start,
        ar.merger.neb.sres_spring(3:end,:),
        ar.merger.neb.sres_end];
    
    q_pathfound = arNEBCheckPath;
    ar.merger.neb.spr(ispr).q_pathfound = q_pathfound;
    
    if q_pathfound
        spring_out = springs(ispr); 
        ar.merger.neb.spring_out = spring_out;
        break          %%%%%%%%%%%%%%% comment here for k profiles
    end
    
end


if qplot
    
    arNEBPlotpair
    drawnow
    if q_pathfound
        print(gcf,sprintf([ar.merger.neb.savefolder '/NEB_LHS%i_to%i_PARS.png'],i_pstart,i_pend), '-dpng');
        saveas(gcf,sprintf([ar.merger.neb.savefolder '/fig/NEB_LHS%i_to%i_PARS'],i_pstart,i_pend));
    else
        print(gcf,sprintf([ar.merger.neb.savefolder '/NEB_LHS%i_to%i_FAILED_PARS.png'],i_pstart,i_pend), '-dpng');
        saveas(gcf,sprintf([ar.merger.neb.savefolder '/fig/NEB_LHS%i_to%i_FAILED_PARS'],i_pstart,i_pend));
    end
    tmpgcf = gcf;
    close(tmpgcf.Number)
    
    
    arNEBAnalyPlot;
    drawnow
    if q_pathfound
        print(gcf,sprintf([ar.merger.neb.savefolder '/NEB_LHS%i_to%i_ANALY.png'],i_pstart,i_pend), '-dpng');
        saveas(gcf,sprintf([ar.merger.neb.savefolder '/fig/NEB_LHS%i_to%i_ANALY'],i_pstart,i_pend));
    else
        print(gcf,sprintf([ar.merger.neb.savefolder '/NEB_LHS%i_to%i_FAILED_ANALY.png'],i_pstart,i_pend), '-dpng');
        saveas(gcf,sprintf([ar.merger.neb.savefolder '/fig/NEB_LHS%i_to%i_FAILED_ANALY'],i_pstart,i_pend));
    end
    tmpgcf = gcf;
    close(tmpgcf.Number)

    
end


end



