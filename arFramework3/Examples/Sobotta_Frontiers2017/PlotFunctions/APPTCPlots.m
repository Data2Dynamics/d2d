function APPTCPlots(png, printFigs)
    
    close all;
    PlotSettings;
    global ar;

    mkdir figures;

    if ~exist('png', 'var')
        png = 0;
    end
    if ~exist('printFigs', 'var')
        printFigs = 0;
    end

	 % Anonymous functions which control export
    if ( printFigs )
        switch( png )
            case 0
                directory = 'figures/Calibration';
                printfig = @(filename)print('-depsc2', '-cmyk', '-loose', '-r300', '-painters', sprintf('%s/%s.eps', directory, filename) );
            case 1
                directory = 'figures/Calibration_png';
                printfig = @(filename)print('-dpng',   '-cmyk', '-loose', '-r300', '-painters', sprintf('%s/%s.png', directory, filename) );
            case 2
                directory = 'figures/Calibration';
                printfig = @(filename)printTikz(sprintf('%s/%s.pdf', directory, filename));
        end
        mkdir( directory );
    else
        printfig = @(filename)disp('No file export!');
    end
    
    % APP Time Course
    clear a;
    for replicate = 1 : 2
        figure('units','normalized','outerposition',[0 0 1 1]);
        xt = [0, 480, 960, 1440];
        
        axes; %#ok
        a(1) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 0), [0 0 0]/255, 'CXCL10_mrna' );
        a(2) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 40), [0 127 255]/255, 'CXCL10_mrna' );
        a(3) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 100), [1, 0.5, 0.0], 'CXCL10_mrna' );
        qLabels( 'CXCL10\_mrna', 'Time (hours)', 'Concentration (log_{10} a.u.)', xt, [-1.5 -1 -0.5 0 .5], 1, 2 );
        xTickHours(xt);
        
        axes; %#ok
        a(1) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 0), [0 0 0]/255, 'FGG_mrna' );
        a(2) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 40), [0 127 255]/255, 'FGG_mrna' );
        a(3) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 100), [1, 0.5, 0.0], 'FGG_mrna' );
        qLabels( 'FGG\_mrna', 'Time (hours)', 'Concentration (log_{10} a.u.)', xt, [-1.5 -1 -0.5 0 .5, 1.0], 1, 1 );
        xTickHours(xt);
        
        axes; %#ok
        a(1) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 0), [0 0 0]/255, 'HAMP_mrna' );
        a(2) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 40), [0 127 255]/255, 'HAMP_mrna' );
        a(3) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 100), [1, 0.5, 0.0], 'HAMP_mrna' );
        qLabels( 'HAMP\_mrna', 'Time (hours)', 'Concentration (log_{10} a.u.)', xt, [-1.5 -1 -0.5 0 .5, 1.0, 1.5], 2, 1 );
        ylim([-1.5, 1.5]);
        xTickHours(xt);
        
        axes; %#ok
        a(1) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 0), [0 0 0]/255, 'IL33_mrna' );
        a(2) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 40), [0 127 255]/255, 'IL33_mrna' );
        a(3) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 100), [1, 0.5, 0.0], 'IL33_mrna' );
        qLabels( 'IL33\_mrna', 'Time (hours)', 'Concentration (log_{10} a.u.)', xt, [-1.5 -1 -0.5 0 .5], 3, 1 );
        xTickHours(xt);
        
        axes; %#ok
        a(1) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 0), [0 0 0]/255, 'APCS_mrna' );
        a(2) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 40), [0 127 255]/255, 'APCS_mrna' );
        a(3) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 100), [1, 0.5, 0.0], 'APCS_mrna' );
        qLabels( 'APCS\_mrna', 'Time (hours)', 'Concentration (log_{10} a.u.)', xt, [-1.2, -1.0 -.8, -.6, -.4, -.2, 0, .2, .4, .6], 1, 0 );
        ylim([-1, .4]);
        xTickHours(xt);
        
        axes; %#ok
        a(1) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 0), [0 0 0]/255, 'HP_mrna' );
        a(2) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 40), [0 127 255]/255, 'HP_mrna' );
        a(3) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 100), [1, 0.5, 0.0], 'HP_mrna' );
        qLabels( 'HP\_mrna', 'Time (hours)', 'Concentration (log_{10} a.u.)', xt, [-1.5 -1 -0.5 0 .25 .5], 2, 0 );
        xTickHours(xt);
        
        axes; %#ok
        a(1) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 0), [0 0 0]/255, 'HPX_mrna' );
        a(2) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 40), [0 127 255]/255, 'HPX_mrna' );
        a(3) = qPlotTC( ar, arFindData( sprintf('xiaoyun_calibration_APP_tc_replicate%d', replicate ), 'input_il6', 100), [1, 0.5, 0.0], 'HPX_mrna' );
        qLabels( 'HPX\_mrna', 'Time (hours)', 'Concentration (log_{10} a.u.)', xt, [-1.5 -1 -0.5 0 .25 .5], 3, 0 );
        xTickHours(xt);
        
        qLegend( a, {'Untreated', 'IL6 (40 ng/mL)', 'IL6 (100 ng/mL)'} );

        printfig( sprintf( 'TC Replicate %d', replicate ) );
    end
end

function xTickHours(xt)
    set(gca, 'XTickLabels', xt/60);
end

% Move legend outside of the figure
function h = qLegend( varargin )
    h = legend( varargin{:}, 'Box', 'off', 'Location', 'Best' );
    k = get( h, 'position' );
    set( h, 'position', [k(1) .2*k(2) k(3:4)] );
end
