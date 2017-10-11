function ValidationPlot(png, printFigs)

    PlotSettings;
    close all;
    global ar;

    rux_color       = [255 0 255]./255;
    sta_color       = [175 36 126]./255;
    control_color   = [0 0 0]./255;

    if ~exist('png', 'var')
        png = 0;
    end
    if ~exist('printFigs', 'var')
        printFigs = 0;
    end
    
    DMSO    = control_color;
    rux     = rux_color;
    sta     = sta_color;

    if ( printFigs )
        switch( png )
            case 1
                dir = 'figures/Validation_png';
                printfig = @(filename)print('-dpng',   '-loose', '-cmyk', '-r300', '-painters', sprintf('%s/%s.png', dir, filename) );
            case 0
                dir = 'figures/Validation';
                printfig = @(filename)print('-depsc2', '-loose', '-cmyk', '-r300', '-painters', sprintf('%s/%s.eps', dir, filename) );
            case 2
                dir = 'Validation';
                printfig = @(filename)printTikz(sprintf('figures/%s/%s.pdf', dir, filename));
        end
        mkdir( dir );
    else
        printfig = @(filename)disp( 'No file export!' );
    end
    
    %% Validation early
    %% Early response Replicate 1
    socs3_rux_wt        = arFindData('braun_validation_Ruxolitinib_1h_hep_2014_04_28_IL6DR_Rux_1h_CXCL10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_rux           = arFindData('braun_validation_Ruxolitinib_1h_hep_2014_04_28_IL6DR_Rux_1h_CXCL10', 'input_ruxolitinib', 500, 'input_stattic', 0);
    socs3_stattic_wt    = arFindData('braun_validation_Stattic_1h_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Cxcl10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_stattic       = arFindData('braun_validation_Stattic_1h_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Cxcl10', 'input_ruxolitinib', 0, 'input_stattic', 60);

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    plh(1) = qPlotDR( ar, socs3_stattic_wt, 'input_il6', 60, DMSO );
    plh(2) = qPlotDR( ar, socs3_stattic, 'input_il6', 60, sta );
    labels( '{\it Cxcl10} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0 0.5], 1 );
    ylim([-1.0 0.5]);

    h(2) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, socs3_rux_wt, 'input_il6', 60, DMSO )
    plh(3) = qPlotDR( ar, socs3_rux, 'input_il6', 60, rux );
    labels( '{\it Cxcl10} mRNA', [-1, 0, 1, 2, 3], [-2.0 -1.5, -1, -0.5, 0 0.5], 2 );
    setLegend(plh, {'DMSO', 'Stattic (60 \muM)', 'Ruxolitinib (500 nM)'} );
    ylim([-1.8 0.5]);
    xlim([-2 3]);
    fixLabels;
    printfig( 'Validation early Replicate 1' );

    %% Early response Replicate 2
    socs3_rux_wt        = arFindData('braun_validation_Ruxolitinib_1h_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Cxcl10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_rux           = arFindData('braun_validation_Ruxolitinib_1h_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Cxcl10', 'input_ruxolitinib', 500, 'input_stattic', 0);
    socs3_stattic_wt    = arFindData('braun_validation_Stattic_1h_hep_2012_04_10_qPCR_140224_IL6DR_Inh_1h_Cxcl10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_stattic       = arFindData('braun_validation_Stattic_1h_hep_2012_04_10_qPCR_140224_IL6DR_Inh_1h_Cxcl10', 'input_ruxolitinib', 0, 'input_stattic', 60);

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    plh(1) = qPlotDR( ar, socs3_stattic_wt, 'input_il6', 60, DMSO );
    plh(2) = qPlotDR( ar, socs3_stattic, 'input_il6', 60, sta );
    labels( '{\it Cxcl10} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0 0.5], 1 );
    ylim([-1.0 0.5]);

    h(2) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, socs3_rux_wt, 'input_il6', 60, DMSO )
    plh(3) = qPlotDR( ar, socs3_rux, 'input_il6', 60, rux );
    labels( '{\it Cxcl10} mRNA', [-1, 0, 1, 2, 3], [-2.0 -1.5, -1, -0.5, 0 0.5], 2 );
    setLegend(plh, {'DMSO', 'Stattic (60 \muM)', 'Ruxolitinib (500 nM)'} );
    ylim([-1.8 0.5]);
    xlim([-2 3]);
    fixLabels;

    printfig( 'Validation early Replicate 2' );

    %% Validation intermediate
    app_rux_wt      = arFindData('braun_validation_Ruxolitinib_6h_hep_2014_04_22_qPCR_140526_IL6DR_Rux_6h_APP', 'input_ruxolitinib', '0');
    app_rux         = arFindData('braun_validation_Ruxolitinib_6h_hep_2014_04_22_qPCR_140526_IL6DR_Rux_6h_APP', 'input_ruxolitinib', '500');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    plh(1) = qPlotDR( ar, app_rux_wt, 'input_il6', 360, DMSO, 'FGG_qpcr' );
    plh(2) = qPlotDR( ar, app_rux, 'input_il6', 360, rux, 'FGG_qpcr' );
    ylim([-1.25 0.25])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 1 );

    h(2) = axes;
    plh(1) = qPlotDR( ar, app_rux_wt, 'input_il6', 360, DMSO, 'HAMP_qpcr' );
    plh(2) = qPlotDR( ar, app_rux, 'input_il6', 360, rux, 'HAMP_qpcr' );
    ylim([-1.5 0.25])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    plh(1) = qPlotDR( ar, app_rux_wt, 'input_il6', 360, DMSO, 'IL33_qpcr' );
    plh(2) = qPlotDR( ar, app_rux, 'input_il6', 360, rux, 'IL33_qpcr' );
    ylim([-1 0.25])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 3 );

    setLegend(plh, {'DMSO', 'Ruxolitinib (500 nM)'} );

    fixLabels;
    printfig( 'Validation intermediate Replicate 1' );

    %% Replicate 2
    app_rux_wt      = arFindData('braun_validation_Ruxolitinib_6h_hep_2014_05_19_qPCR_140526_IL6DR_Rux_6h_APP', 'input_ruxolitinib', '0');
    app_rux         = arFindData('braun_validation_Ruxolitinib_6h_hep_2014_05_19_qPCR_140526_IL6DR_Rux_6h_APP', 'input_ruxolitinib', '500');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    plh(1) = qPlotDR( ar, app_rux_wt, 'input_il6', 360, DMSO, 'FGG_qpcr' );
    plh(2) = qPlotDR( ar, app_rux, 'input_il6', 360, rux, 'FGG_qpcr' );
    ylim([-1.25 0.25])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 1 );

    h(2) = axes;
    plh(1) = qPlotDR( ar, app_rux_wt, 'input_il6', 360, DMSO, 'HAMP_qpcr' );
    plh(2) = qPlotDR( ar, app_rux, 'input_il6', 360, rux, 'HAMP_qpcr' );
    ylim([-1.5 0.25])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    plh(1) = qPlotDR( ar, app_rux_wt, 'input_il6', 360, DMSO, 'IL33_qpcr' );
    plh(2) = qPlotDR( ar, app_rux, 'input_il6', 360, rux, 'IL33_qpcr' );
    ylim([-1 0.25])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 3 );

    setLegend(plh, {'DMSO', 'Ruxolitinib (500 nM)'} );
    fixLabels;
    printfig( 'Validation intermediate Replicate 2' );

    %% Validation late
    %% Replicate 1
    app_rux_wt      = arFindData('braun_validation_Ruxolitinib_24h_hep_2014_04_22_qPCR_140604_IL6DR_Rux_24h_APP', 'input_ruxolitinib', '0');
    app_rux         = arFindData('braun_validation_Ruxolitinib_24h_hep_2014_04_22_qPCR_140604_IL6DR_Rux_24h_APP', 'input_ruxolitinib', '500');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'FGG_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'FGG_qpcr' );
    ylim([-1.4 0.5])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 1 );

    h(2) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'HAMP_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'HAMP_qpcr' );
    ylim([-1.5 0.35])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'IL33_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'IL33_qpcr' );
    ylim([-.7 0.1])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 3 );

    h(1) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'APCS_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'APCS_qpcr' );
    ylim([-1.2 0.1])
    labels( '{\it apcs} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1, 2 );

    h(2) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'HP_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'HP_qpcr' );
    ylim([-.6 0.1])
    labels( '{\it hp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 2, 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'HPX_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'HPX_qpcr' );
    ylim([-0.5 0.1])
    labels( '{\it hpx} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 3, 2 );    

    fixLabels;
    printfig( 'Validation late Replicate 1' );

    %% Replicate 2
    app_rux_wt      = arFindData('braun_validation_Ruxolitinib_24h_hep_2014_05_19_qPCR_140604_IL6DR_Rux_24h_APP', 'input_ruxolitinib', '0');
    app_rux         = arFindData('braun_validation_Ruxolitinib_24h_hep_2014_05_19_qPCR_140604_IL6DR_Rux_24h_APP', 'input_ruxolitinib', '500');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'FGG_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'FGG_qpcr' );
    ylim([-1.4 0.5])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 1 );

    h(2) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'HAMP_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'HAMP_qpcr' );
    ylim([-1.5 0.35])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0, 0.5], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'IL33_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'IL33_qpcr' );
    ylim([-.5 0.1])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 3 );

    h(1) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'APCS_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'APCS_qpcr' );
    ylim([-1.2 0.1])
    labels( '{\it apcs} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1, 2 );

    h(2) = axes;
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'HP_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'HP_qpcr' );
    ylim([-.6 0.1])
    labels( '{\it hp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 2, 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, app_rux_wt, 'input_il6', 1440, DMSO, 'HPX_qpcr' );
    qPlotDR( ar, app_rux, 'input_il6', 1440, rux, 'HPX_qpcr' );
    ylim([-0.5 0.1])
    labels( '{\it hpx} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 3, 2 );    

    fixLabels;
    printfig( 'Validation late Replicate 2' );
end

% Export setup, 20x20 cm, painters/vector 300 dpi RGB
function setLegend( handles, legendItems )
    legHandle           = legend(handles, legendItems);
    legPos              = get( legHandle, 'position' );
    set( legHandle, 'position', [0.1 0.1 legPos(3:4)] );
end

function fixLabels( )
    global allTicks;
    allAx = findall(gcf,'type','axes');
    
    for a = 1 : length( allAx )
        allAx(a)
        axes(allAx(a)); %#ok
        xt = get( gca, 'XTick' );

        XL = get( gca, 'xlim' );
        xt(xt<XL(1)) = [];
        xt(xt>XL(end)) = [];
        
        l = {};
        if ( allTicks )
            for b = 1 : length(xt)
                if ( xt(b) == floor(xt(b)) )
                    l{b} = sprintf('10^{%d}',xt(b) ); %#ok
                else
                    l{b} = sprintf('0',xt(b) ); %#ok
                end
            end
            try
                xticklabel_rotate(xt, 0, l );
            catch
            end
        end
    end
end

function labels( lab, xt, yt, N, N2, xM )

    if (nargin < 5)
        N2 = 1;
    end
    
    xL = xt;
    if ( nargin > 5 )
        xL( xt == xM ) = -inf;
    end

    ylabel('log_{10}(conc) [au]');
    xlabel('IL-6 [ng/mL]');
        
    box off;
    set( gca, 'XTickMode', 'manual' );
    set( gca, 'YTickMode', 'manual' );    
    set( gca, 'XTick', xt );
    set( gca, 'XTickLabel', 10.^xL );    
    set( gca, 'YTick', yt );
    set( gca, 'YTickLabel', yt );    
    set( gca, 'TickDir', 'out' );
    set( gca, 'TickLength', [0.025, 0.025] );
    set( gca, 'Position', [0.1 + (N-1) * 0.3, 0.4 + (N2-1)*0.25, 0.2, 0.15] );
    title( lab );
    
    set(gcf, 'Clipping', 'off' );
    set(gcf, 'Renderer', 'painters' );
    set(gcf, 'RendererMode', 'manual' );
    set(gcf, 'PaperType', '<custom>' );
    set(gcf, 'PaperUnits', 'centimeters' );
    set(gcf, 'PaperSize', [20, 20] )
    set(gcf, 'Position', [0 0 20 20] );
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 20 20] );
    set(gcf, 'Resize', 'On' );
end
