function Plot2DPL(plotmode)
% Plot2DPL(plotmode)
%
% Plots the 2D-Profile in 2D as a scatter plot (1) or a contour plot (2) 
% depending on flag plotmode. 

global ar

if ~exist('plotmode','var') 
    fprintf('\n ERROR Plot2DPL: Please specify argument plotmode. \n')
    return
end

if plotmode == 1 %Scatter Plot
    %Unprocessed data points from ar.ple2d:
    x = ar.ple2d.raw.plpar;
    y = ar.ple2d.raw.predsteps;
    z = ar.ple2d.raw.chi2;
    
    y = repmat(y,size(x,1),1);
    y = y(:);
    %Values in vector y now exactly fit to the corresponding values in x
    x = x(:);
    z = z(:);
    
    scatter3(x(:),y(:),z(:),ar.ple2d.config.plot.thicknessscatter,z(:),'Marker','.')
    set(gca,'GridAlpha',0.3)
    

else %Contour Plot
    %Get smooth 2d-profile:
    if ~isfield(ar.ple2d,'smooth') 
    	disp('ERROR Plot2DPL: Run smooth2d to generate the smooth 2D-landscape')
        return
    end
    xq = ar.ple2d.smooth.xq;
    yq = ar.ple2d.smooth.yq;
    zq = ar.ple2d.smooth.zq;

    %Trajectory of the validation profile:     
    [x_vpl,y_vpl,z_vpl] = vplfrom2d;

    %Minimum-index:
    [~,ind_minval] = min(z_vpl);
    
    % Plot everything:
    hold on
    contour(xq,yq,zq,ar.ple2d.config.plot.ncontourlines) % Basic contour plot
    line([min(min(xq)),max(max(xq))],[y_vpl(ind_minval),...
        y_vpl(ind_minval)],'Color','black','LineStyle','--','LineWidth',1.5)
    % Indicates optimal prediction
    plot(x_vpl,y_vpl,'Color','black','LineStyle',':','LineWidth',1.5);
    % Validation profile
    scatter([min(ar.ple2d.raw.plpar),max(ar.ple2d.raw.plpar)],...
        [ar.ple2d.raw.predsteps,ar.ple2d.raw.predsteps],'black','x')
    % Hard boundaries of the computed parameter profiles
    
    % Trajectory of boundary values:
    if ~isfield(ar.ple2d,'bounds') || ~isfield(ar.ple2d.bounds,'bounds')
        fprintf(['WARNING Plot2DPL: Run bounds2d to generate',...
            ' \n confidence intervals with specified confidence levels. \n']);
    else
        try
            plot(ar.ple2d.bounds.bounds,ar.ple2d.bounds.pred_points,'ro')
        catch exception
            fprintf(['ERROR Plot2DPL: Try calculating the bounds by using bounds2d ',...
                'again with ar.ple2d.config.bounds.save_mode = 2.',... 
                '\n Error message: \n %s \n Line: %s \n'] ,...
                exception.message, sprintf('%i, ',exception.stack.line));
        end
        %Add boundaries of original profile:
        try
            opt_pred = find(ar.ple2d.raw.levels == 0);
            if ar.ple2d.bounds.alpha_bounds == 1-ar.ple.alpha_level 
                %Check whether 1D and 2D confidence levels were chosen consistently
                scatter([ar.ple.conf_lb_point(ar.ple2d.general.idpar),...
                    ar.ple.conf_ub_point(ar.ple2d.general.idpar)],...
                    [ar.ple2d.raw.predsteps(opt_pred),ar.ple2d.raw.predsteps(opt_pred)],...
                    80,'black','o','filled')
            else
                fprintf(['WARNING Plot2DPL: Original boundaries are not consistent with ',...
                    'the chosen 2D-boundaries, \n thus they are not plotted. \n'])
            end
        catch exception
            fprintf(['WARNING Plot2DPL: Original boundaries could not be plotted.\n',... 
                '\n Error message: \n %s \n Line: %s \n'] ,...
                exception.message, sprintf('%i, ',exception.stack.line));
        end
    end
    
    hold off
    
end

set(get(gca,'XAxis'),'LineWidth',0.9);
set(get(gca,'YAxis'),'LineWidth',0.9);
set(get(gca,'ZAxis'),'LineWidth',0.9);
xlabel('par','FontSize',12);
ylabel('pred','FontSize',12);
zlabel('\chi^2','FontSize',12);
end
