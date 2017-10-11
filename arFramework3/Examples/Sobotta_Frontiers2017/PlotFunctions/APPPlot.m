function APPPlot(png, printfigs)
    
    close all;
    PlotSettings;
    global ar;

    if ~exist('png', 'var')
        png = 0;
    end
    if ~exist('printfigs', 'var')
        printfigs = 0;
    end
    
    % Set colors
    WT_IL6      = [0 127 255]./255;
    rux         = [255 0 255]./255;
    sta         = [175 36 126]./255;
    DMSO        = [0 0 0]./255;

    if ( printfigs )
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
        printfig = @(filename)disp('No file export!'); %#ok
    end
     
    
    %% Early response
    socs3_rux_wt        = arFindData('braun_app_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_rux           = arFindData('braun_app_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 500, 'input_stattic', 0);
    socs3_stattic_wt    = arFindData('braun_app_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_stattic       = arFindData('braun_app_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 60);

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    plh(1) = qPlotDR( ar, socs3_stattic_wt, 'input_il6', 60, DMSO );
    plh(2) = qPlotDR( ar, socs3_stattic, 'input_il6', 60, sta );
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, socs3_rux_wt, 'input_il6', 60, DMSO )
    plh(3) = qPlotDR( ar, socs3_rux, 'input_il6', 60, rux );
    ylim([-1.8, 0.25]);
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 2 );
    setLegend(plh, {'DMSO', 'Stattic (60 \muM)', 'Ruxolitinib (500 nM)'} );

    fixLabels;
    printfig( 'Early response [old]' );

    %% Early response Replicate 1
    socs3_rux_wt        = arFindData('braun_app_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_rux           = arFindData('braun_app_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 500, 'input_stattic', 0);
    socs3_stattic_wt    = arFindData('braun_app_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_stattic       = arFindData('braun_app_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 60);

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    plh(1) = qPlotDR( ar, socs3_stattic_wt, 'input_il6', 60, DMSO, [], -1 );
    plh(2) = qPlotDR( ar, socs3_stattic, 'input_il6', 60, sta, [], -1 );
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1, 1, -1 );
    ylim( [-1.65, .25] );

    h(2) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, socs3_rux_wt, 'input_il6', 60, DMSO )
    plh(3) = qPlotDR( ar, socs3_rux, 'input_il6', 60, rux );
    ylim([-1.8, 0.25]);
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 2 );
    setLegend(plh, {'DMSO', 'Stattic (60 \muM)', 'Ruxolitinib (500 nM)'} );

    cxcl10_wt        = arFindData('braun_calibration_hep_2013_10_14_qPCR_140109_IL6DR_1h_CXCL10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, cxcl10_wt, 'input_il6', 60, WT_IL6, [] );
    if ( ~allTicks )
        labels( '{\it CXCL10} mRNA', [-3, -1,  1, 3], [-1.5, -1, -0.5, 0], 3, 1 );
    else
        labels( '{\it CXCL10} mRNA', [-3, -2, -1,  0, 1, 2, 3], [-1.5, -1, -0.5, 0], 3, 1 );
    end
    ylim([-1.3, 0.25]);
    setLegend(plh, {'DMSO', 'Stattic', 'Ruxolitinib', 'WT'} );

    text(h(1),0,-10,'Early response 10_14 02_14 and 10_14');

    fixLabels;
    printfig( 'Early response 10_14 02_14 and 10_14');

    %% Early response Replicate 2
    socs3_rux_wt        = arFindData('braun_app_hep_2014_04_28_IL6DR_Rux_1h_SOCS3', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_rux           = arFindData('braun_app_hep_2014_04_28_IL6DR_Rux_1h_SOCS3', 'input_ruxolitinib', 500, 'input_stattic', 0);
    socs3_stattic_wt    = arFindData('braun_app_hep_2012_04_10_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 0);
    socs3_stattic       = arFindData('braun_app_hep_2012_04_10_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 60);

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    plh(1) = qPlotDR( ar, socs3_stattic_wt, 'input_il6', 60, DMSO );
    plh(2) = qPlotDR( ar, socs3_stattic, 'input_il6', 60, sta );
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, socs3_rux_wt, 'input_il6', 60, DMSO )
    plh(3) = qPlotDR( ar, socs3_rux, 'input_il6', 60, rux );
    ylim([-2.0, 0.25]);
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-2.0, -1.5, -1, -0.5, 0], 2 );
    setLegend(plh, {'DMSO', 'Stattic (60 \muM)', 'Ruxolitinib (500 nM)'} );

    cxcl10_wt        = arFindData('braun_calibration_hep_2011_06_06_qPCR_140109_IL6DR_1h_CXCL10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, cxcl10_wt, 'input_il6', 60, WT_IL6, [], -4.5 );
    if ( ~allTicks )
        labels( '{\it CXCL10} mRNA', [-4, -2,  1, 2], [-1.5, -1, -0.5, 0], 3, 1, -4.5 );
    else
        labels( '{\it CXCL10} mRNA', [-4, -2, 0, 2], [-1.5, -1, -0.5, 0], 3, 1, -4.5 );
    end
    ylim([-1, 0.25]);
    setLegend(plh, {'DMSO', 'Stattic', 'Ruxolitinib', 'WT'} );

    fixLabels;
    printfig( 'Early response 04_28 04_10 and 06_06');      

    %% Early response (2)
    %% Replicate 1
    figure;
    wt = arFindData('braun_calibration_hep_2011_06_06_qPCR_140109_IL6DR_1h_CXCL10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, wt, 'input_il6', 60, WT_IL6, [], -4.5 );
    labels( '{\it CXCL10} mRNA', [-4 -2, 0, 2], [-1.5, -1, -0.5, 0], 3, 1, -4.5 );
    wt = arFindData('braun_app_hep_2011_06_06_qPCR_140109_IL6DR_1h', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, wt, 'input_il6', 60, WT_IL6, [], -4.5 );
    labels( '{\it Socs3} mRNA', [-4 -2, 0, 2], [-2.0, -1.5, -1, -0.5, 0], 2 );
    printfig( 'Early response replicate 1');  
    
    figure;
    wt = arFindData('braun_calibration_hep_2013_10_14_qPCR_140109_IL6DR_1h_CXCL10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, wt, 'input_il6', 60, WT_IL6, [], -4.5 );
    labels( '{\it CXCL10} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 3, 1, -4.5 );
    wt = arFindData('braun_app_hep_2013_10_14_qPCR_140109_IL6DR_1h', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, wt, 'input_il6', 60, WT_IL6, [], -4.5 );
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-2.0, -1.5, -1, -0.5, 0], 2 );
    printfig( 'Early response replicate 2');  
    
    figure;
    wt = arFindData('braun_calibration_hep_2013_10_21_qPCR_140109_IL6DR_1h_CXCL10', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, wt, 'input_il6', 60, WT_IL6, [], -4.5 );
    labels( '{\it CXCL10} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 3, 1, -4.5 );
    wt = arFindData('braun_app_hep_2013_10_21_qPCR_140109_IL6DR_1h', 'input_ruxolitinib', 0, 'input_stattic', 0);
    h(3) = axes;
    plh(4) = qPlotDR( ar, wt, 'input_il6', 60, WT_IL6, [], -4.5 );
    labels( '{\it Socs3} mRNA', [-1, 0, 1, 2, 3], [-2.0, -1.5, -1, -0.5, 0], 2 );
    printfig( 'Early response replicate 3');      
    
    %% Intermediate response
    %% Replicate 1
    APP = arFindData('braun_calibration_hep_2014_01_13_IL6DR_6h_Replicate1');
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'FGG_qpcr' );
    ylim([-1 0.1])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'HAMP_qpcr' );
    ylim([-.9 0.1])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.66, -0.33, 0], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'IL33_qpcr' );
    ylim([-.53 0.05])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25 0], 3 );

    fixLabels;
    printfig( 'Intermediate response replicate 1');      

    %% Replicate 2
    APP = arFindData('braun_calibration_hep_2014_01_13_IL6DR_6h_Replicate2');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'FGG_qpcr' );
    ylim([-1 0.1])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'HAMP_qpcr' );
    ylim([-.8 0.1])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.66, -0.33, 0], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'IL33_qpcr' );
    ylim([-.65 0.05])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25 0], 3 );

    fixLabels;
    printfig( 'Intermediate response replicate 2');  

    %% Replicate 3
    APP = arFindData('braun_calibration_hep_2014_01_13_IL6DR_6h_Replicate3');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'FGG_qpcr' );
    ylim([-1.2 0.1])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'HAMP_qpcr' );
    ylim([-.8 0.1])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.66, -0.33, 0], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 360, WT_IL6, 'IL33_qpcr' );
    ylim([-.65 0.05])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25 0], 3 );

    fixLabels;
    printfig( 'Intermediate response replicate 3');    

    %% Late response
    %% Replicate 1
    APP = arFindData('braun_calibration_hep_2014_01_27_IL6DR_24h_Replicate1');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'FGG_qpcr' );
    ylim([-1.4 0.25])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HAMP_qpcr' );
    ylim([-1.2 0.1])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'IL33_qpcr' );
    ylim([-0.5 0.1])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 3 );

    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'APCS_qpcr' );
    ylim([-1.2 0.1])
    labels( '{\it apcs} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1, 2 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HP_qpcr' );
    ylim([-.6 0.1])
    labels( '{\it hp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 2, 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HPX_qpcr' );
    ylim([-0.5 0.1])
    labels( '{\it hpx} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 3, 2 );    

    fixLabels;
    printfig( 'Late response Replicate 1');      

    %% Replicate 2
    APP = arFindData('braun_calibration_hep_2014_01_27_IL6DR_24h_Replicate2');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'FGG_qpcr' );
    ylim([-1.6 0.25])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HAMP_qpcr' );
    ylim([-1.01 0.1])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'IL33_qpcr' );
    ylim([-0.5 0.1])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 3 );

    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'APCS_qpcr' );
    ylim([-1 0.1])
    labels( '{\it apcs} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1, 2 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HP_qpcr' );
    ylim([-.6 0.1])
    labels( '{\it hp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 2, 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HPX_qpcr' );
    ylim([-0.4 0.1])
    labels( '{\it hpx} mRNA', [-1, 0, 1, 2, 3], [-.8, -.6, -0.4, -0.2, 0], 3, 2 );    

    fixLabels;
    printfig( 'Late response Replicate 2');
    
    %% Replicate 3
    APP = arFindData('braun_calibration_hep_2014_01_27_IL6DR_24h_Replicate3');

    figure('units','normalized','outerposition',[0 0 1 1]);
    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'FGG_qpcr' );
    ylim([-1.4 0.25])
    labels( '{\it fgg} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HAMP_qpcr' );
    ylim([-1.2 0.1])
    labels( '{\it hamp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'IL33_qpcr' );
    ylim([-0.5 0.1])
    labels( '{\it il33} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 3 );

    h(1) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'APCS_qpcr' );
    ylim([-1.2 0.1])
    labels( '{\it apcs} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, 0], 1, 2 );

    h(2) = axes;
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HP_qpcr' );
    ylim([-.6 0.1])
    labels( '{\it hp} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 2, 2 );

    h(3) = axes;
    linkaxes(h, 'x');
    qPlotDR( ar, APP, 'input_il6', 1440, WT_IL6, 'HPX_qpcr' );
    ylim([-0.5 0.1])
    labels( '{\it hpx} mRNA', [-1, 0, 1, 2, 3], [-1.5, -1, -0.5, -0.25, 0], 3, 2 );    

    fixLabels;
    printfig( 'Late response Replicate 3');  
   
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
    
    ylabel('log_{10}(conc) [au]');
    xlabel('IL-6 [ng/mL]');
end
