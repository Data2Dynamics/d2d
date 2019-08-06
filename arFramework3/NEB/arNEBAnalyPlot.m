function q_pathfound = arNEBAnalyPlot
% Analysis plot and quality control polts for the current NEB for all
% spring constants


global ar

figure

nrow = 5;
ncol = 4;

legcell = {};

legcell{1} =  ['direct path ' ];

subplot(nrow,ncol,1)

plot( 1, 0, '-ok' )
hold on

tmpchi = ar.merger.neb.ini.chi2_step;

subplot(nrow, ncol, [2 3])
plot(tmpchi, '-ok')
hold on
plot([1 length(tmpchi)],[tmpchi(end) tmpchi(end) ] , '.k')
plot([1 length(tmpchi)],[tmpchi(1) tmpchi(1) ] , '.k')
ylabel( 'chi 2' )
xlabel( 'steps' )
title('paths')
 

tmpresult =  ar.merger.neb.ps_init;
for i = 1:size(tmpresult,1)-1
    tmpstepsize(i) = norm(tmpresult(i,:) - tmpresult(i+1,:));
end

subplot(nrow, ncol, 5)
plot(tmpstepsize, '-xk')
hold on
ylabel('stepsize')
xlabel('step')

subplot(nrow, ncol, 13)
plot(tmpstepsize-mean(tmpstepsize), '-xk')
hold on
ylabel('stepsize')
xlabel('step')


% below start and end?   
init_q_below_chi2 = sum(tmpchi > ...
    max(tmpchi(1),tmpchi(end))) < 1;

% tolerance for monotocitiy test
init_q_monotonic = check_monotonicity_path(tmpchi, ...
    max( abs(tmpchi(1)-tmpchi(end))/50, 1e-4));


subplot(nrow, ncol, 17)
plot(-10, init_q_below_chi2, '-xk')
hold on

   subplot(nrow, ncol, 18)
plot(-10, init_q_monotonic, '-xk')
hold on



for i_spr = 1:length(ar.merger.neb.spr)
    
    legcell{i_spr+1} =  ['k = ' num2str(ar.merger.neb.spr(i_spr).springconst)];
    subplot(nrow,ncol,1)
    plot( 1, ar.merger.neb.spr(i_spr).springconst,'-o' )
    hold on
    legend(legcell,'Location','southoutside')  
    
    tmpchi = ar.merger.neb.spr(i_spr).chi2s;
    
    subplot(nrow, ncol, [2 3])
    plot(tmpchi, '-*')
    hold on
    plot([1 length(tmpchi)],[tmpchi(end) tmpchi(end) ] , '--k')
    plot([1 length(tmpchi)],[tmpchi(1) tmpchi(1) ] , '--k')
    ylabel( 'chi 2' )
    xlabel( 'steps' )
    title('paths')
 
    chi2_sum(i_spr) = sum(tmpchi);
    chi2_max(i_spr) = max(tmpchi);
    chi2_sum_norm(i_spr) = sum(tmpchi)/length(tmpchi);
    chi2_toolarge(i_spr) = sum(tmpchi > max(tmpchi(1),tmpchi(end)));
    
    
    subplot(nrow, ncol, 4)
    plot(ar.merger.neb.spr(i_spr).res./ar.merger.neb.spr(i_spr).springconst,'-o')
    hold on
    title('residuals / springsize')
    
    subplot(nrow, ncol, 12)
    plot(abs(ar.merger.neb.spr(i_spr).res)./ar.merger.neb.spr(i_spr).springconst,'-o')
    hold on
    title('residuals / springsize')
    
    tmpresult =  ar.merger.neb.spr(i_spr).ps_result;
    for i = 1:size(tmpresult,1)-1
        tmpstepsize(i) = norm(tmpresult(i,:) - tmpresult(i+1,:));
    end
    subplot(nrow, ncol, 5)
    plot(tmpstepsize, '-x')
    hold on
    ylabel('stepsize')
    xlabel('step')
    
    subplot(nrow, ncol, 13)
    plot(tmpstepsize-mean(tmpstepsize), '-x')
    hold on
    ylabel('stepsize')
    xlabel('step')
    
    subplot(nrow, ncol, 6)
    hist(tmpstepsize)
    hold on
    xlabel('stepsize')

    springconst(i_spr) = ar.merger.neb.spr(i_spr).springconst;
    step_mean(i_spr) = mean(tmpstepsize);
    step_std(i_spr) = std(tmpstepsize);
    step_cov(i_spr) = std(tmpstepsize)/mean(tmpstepsize);
    
    
    % below start and end?   
    q_below_chi2(i_spr) = sum(tmpchi > ...
        max(tmpchi(1),tmpchi(end))) < 1;
    
    % tolerance for monoticitiy test
    q_monotonic(i_spr) = check_monotonicity_path(tmpchi, ...
        max( abs(tmpchi(1)-tmpchi(end))/50, 1e-4));
    
    % stepsize not too extreme
    q_stepextreme(i_spr) = ...
    sum( tmpstepsize > 0.3 * ( norm(tmpresult(1,:) - tmpresult(end,:))) )  < 1;
    
    
end


if exist('springconst')

    subplot(nrow, ncol, 7)
    plot(springconst, step_mean, 'o-k')
    hold on
    plot(springconst, step_std, 'x:k')
    plot(springconst, step_cov, '*-k')
    
    xlabel('step')
    title('step sizes : o mean, x, std, * c.o.v')
    
    subplot(nrow, ncol, 8)
    plot(springconst, chi2_sum, 'o-k')
    hold on
    title('chi2 summe')
    
    subplot(nrow, ncol, 9)
    plot(springconst, chi2_max, 'x:k')
        title('chi2 max')
    
    subplot(nrow, ncol, 10)
    plot(springconst, chi2_sum_norm, '*-k')
        title('chi2 normierte summe')
        
    subplot(nrow, ncol, 11)
    plot(springconst, chi2_toolarge, '*-k')
        title('chi2 groesser endpunkte')

        
        
    subplot(nrow, ncol, 17)
    plot(-10, init_q_below_chi2, '-xk')
    hold on
    plot(springconst, q_below_chi2, '-o')
    hold on 
    title('below chi2')
    
    subplot(nrow, ncol, 18)
    plot(-10, init_q_monotonic, '-xk')
    hold on
    plot(springconst, q_monotonic, '-o')
    hold on     
    title('monotonic')
    
    subplot(nrow, ncol, 19)   
    plot(springconst, q_stepextreme, '-o')
    hold on     
    title('stepsize not too extreme')
    
   	subplot(nrow, ncol, 20)   
    plot(springconst, q_stepextreme .* q_monotonic .* q_below_chi2, '-o')
    hold on     
    title('total good path')
 
    
    q_pathfound = sum(q_stepextreme .* q_monotonic .* q_below_chi2) > 0;
end


end

