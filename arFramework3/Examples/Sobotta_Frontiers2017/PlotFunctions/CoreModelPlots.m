function CoreModelPlots( png, printFigs )
   
    PlotSettings;
    close all;
    global ar;
    
    % Plot colors
    WT_color        = [0 127 255]./255;
    rux_color       = [255 0 255]./255;
    sta_color       = [175 36 126]./255;
    control_color   = [0 0 0]./255;

    % Output directory
    directory = 'Core';
    mkdir figures;

    % Output format
    if ~exist( 'png', 'var' )
        png = 0;
    end
    
    % Export figures?
    if ~exist( 'printFigs', 'var' )
        printFigs = 0;
    end    

	 % Set up export function
    if ( printFigs )
        switch( png )
            case 0
                % eps
                mkdir( sprintf( 'figures/%s', directory ) );
                printfig = @(filename)print('-depsc2', '-loose', '-cmyk', '-r300', '-painters', sprintf('figures/%s/%s.eps', directory, filename) );
            case 1
                %png
                directory = [directory '_png'];
                mkdir( sprintf( 'figures/%s', directory ) );
                printfig = @(filename)print('-dpng', '-loose', '-cmyk', '-r300', '-painters', sprintf('figures/%s/%s.png', directory, filename) );
            case 2
                mkdir( sprintf( 'figures/%s', directory ) );
                printfig = @(filename)printTikz(sprintf('figures/%s/%s.pdf', directory, filename));
        end
    else
        printfig = @(filename)disp('No file export!'); %#ok
    end

    % Core model figures pt. 1
    clear a;
    figure('units','normalized','outerposition',[0 0 1 1]);
    a = axes; qPlotTC( ar, arFindData( '2005_02_03'), WT_color, 'pJAK1_wb' );
    qLabels( 'phospho-JAK1', '', '', [0 20 40 60 80 100], [-1 -0.5 0], 1, 2 );
    qLegend( a, {'Untreated', 'IL6 (40 ng/mL)'} );
    axes; qPlotTC( ar, arFindData( '2005_02_03'), WT_color, 'pgp130_wb' );
    qLabels( 'phospho-gp130', '', '', [0 20 40 60 80 100], [-1 -0.5 0], 2, 2 );
    axes; qPlotTC( ar, arFindData( '2005_02_03'), WT_color, 'pSTAT3_wb' );
    qLabels( 'phospho-STAT3', '', '', [0 20 40 60 80 100], [-0.8, -0.6, -0.4, -0.2, 0 0.2], 1, 1 );
    axes; a(1) = qPlotTC( ar, arFindData( '2013_07_08', 'input_il6', 0), control_color, 'SOCS3_qpcr' );
    a(2) = qPlotTC( ar, arFindData( '2013_07_08', 'input_il6', 40), WT_color, 'SOCS3_qpcr' );
    qLabels( 'Socs3 mRNA', '', '', [0 20 40 60 80 100 120], [-2 -1.5 -1 -0.5 0], 2, 1 );
    pyData = arFindData( '2011_04_04');
    axes; b(1)=qPlotTC( ar, pyData(1), control_color, 'pSTAT3_ms' );
    hold on; b(2)=qPlotTC( ar, pyData(2), WT_color, 'pSTAT3_ms' );
    qLabels( 'phospho-STAT3 (Tyr-705)', 'Time (min)', 'Degree of phosphorylation (%)', [0 10 20], [0, 20, 40, 60], 3, 2 );
    axes; qPlotTC( ar, arFindData( 'xiaoyun_nucSTAT3_Ratio' ), WT_color, 'ntSTAT3_ratio' );
    qLabels( 'mKate2-STAT3', 'Time (min)', 'Nuc/cyt ratio', [-10 0 10 20 30], [0.3 1 3.2 10], 3, 1, [0,1] );
    printfig( 'Core Model A' );

    % Core model figures pt 2
    clear a;
    figure('units','normalized','outerposition',[0 0 1 1]);
    axes; qPlotDR( ar, arFindData( '2011_12_12'), 'input_stattic', 20, WT_color, 'pSTAT3_wb', -1 ); hold on;
    qPlotDR( ar, arFindData( '2011_12_12'), 'input_stattic', 0, control_color, 'pSTAT3_wb', -1 ); ylim([-2.25 .75]);
    qLabels( 'phospho-STAT3', 'Stattic (\mu M)', 'Concentration (log_{10} a.u.)', [0, 1, 5, 10, 20, 50], [-2. -1. 0], 1, 2, [1, 0] );
    axes; a(2) = qPlotDR( ar, arFindData( '2013_06_03', 'input_il6', 40), 'input_ruxolitinib', [], WT_color, 'pSTAT3_wb', 0 ); hold on;
    a(1) = qPlotDR( ar, arFindData( '2013_06_03', 'input_il6', 0), 'input_ruxolitinib', [], control_color, 'pSTAT3_wb', 0 ); ylim([-4 1]);
    qLabels( 'phospho-STAT3', 'Ruxolitinib (\mu M)', 'Concentration (log_{10} a.u.)', [0, 10, 100, 1000, 10000], [-4, -2, 0], 1, 1, [1, 0] );
    xo = get(gca,'XLim'); xlim([xo(1), 4.1]);
    qLegend( a, {'Untreated', 'IL6 (40 ng/mL), 20 min'} );
    printfig( 'Core Model B' );
    
    % Core model figures pt 3
    clear a;
    figure('units','normalized','outerposition',[0 0 1 1]);
    axes; a(1) = qPlotTC( ar, arFindData( 'braun_hep_2012_04_10_Stattic_Stat3_Inhibitor_TC', 'input_stattic', 0, 'input_il6', 50), control_color, 'pSTAT3_wb' );
    a(2) = qPlotTC( ar, arFindData( 'braun_hep_2012_04_10_Stattic_Stat3_Inhibitor_TC', 'input_stattic', 60), sta_color, 'pSTAT3_wb' );
    qLabels( 'phospho-STAT3', 'Time (min)', 'Concentration (log_{10} a.u.)', [-60 -40 -20 0 20 40], [-1.5, -1, -0.5, 0], 1, 2 ); ylim( [-1.7, 0.3] );

    axes; a(1) = qPlotTC( ar, arFindData( '2013_09_30', 'input_ruxolitinib', 0), control_color, 'SOCS3_qpcr' );
    a(2) = qPlotTC( ar, arFindData( '2013_09_30', 'input_ruxolitinib', 500), rux_color, 'SOCS3_qpcr' );
    qLabels( 'Socs3 mRNA', 'Time (min)', 'Concentration (log_{10} a.u.)', [-60 -30 0 30 60 90 120], [-2 -1 0], 2, 2 ); ylim( [-2.5, 0.5] );

    axes; a(1) = qPlotDR( ar, arFindData( 'braun_app_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 60), 'input_il6', [], sta_color, 'SOCS3_qpcr', -1 ); hold on;
    a(2) = qPlotDR( ar, arFindData( 'braun_app_hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 0), 'input_il6', [], control_color, 'SOCS3_qpcr', -1 );
    qLabels( 'Socs3 mRNA (1 h)', 'IL6 (ng/mL)', 'Concentration (log_{10} a.u.)', [0, 1, 10, 100, 1000], [-1.5, -1, -0.5, 0], 1, 1, [1, 0] );
    qLegend( a, {'DMSO', 'Stattic (60 \mu M)'} );
    xlim([-1.3, 3.1]);

    axes; a(1) = qPlotDR( ar, arFindData( 'braun_app_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 500, 'input_stattic', 0), 'input_il6', [], rux_color, 'SOCS3_qpcr', -1 ); hold on;
    a(2) = qPlotDR( ar, arFindData( 'braun_app_hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 'input_ruxolitinib', 0, 'input_stattic', 0), 'input_il6', [], control_color, 'SOCS3_qpcr', -1 );
    qLabels( 'Socs3 mRNA (1 h)', 'IL6 (ng/mL)', 'Concentration (log_{10} a.u.)', [0.1, 1, 10, 100, 1000], [-1.5, -1, -0.5, 0], 2, 1, [1, 0] ); ylim([-1.8, 0.25]);
    qLegend( a, {'DMSO', 'Ruxolitinib (500 nM)'} );
    xlim([-1.3, 3.1]);

    axes; a(1) = qPlotDR( ar, arFindData( 'app_hep_2014_05_19_qPCR_140526', 'input_ruxolitinib', 500, 'input_stattic', 0), 'input_il6', 360, rux_color, 'SOCS3_qpcr', -1 ); hold on;
    a(2) = qPlotDR( ar, arFindData( 'app_hep_2014_05_19_qPCR_140526', 'input_ruxolitinib', 0, 'input_stattic', 0), 'input_il6', 360, control_color, 'SOCS3_qpcr', -1 );
    qLabels( 'Socs3 mRNA (6 h)', 'IL6 (ng/mL)', 'Concentration (log_{10} a.u.)', [0.1, 1, 10, 100, 1000], [-1.5, -1, -0.5, 0], 3, 2, [1, 0] ); ylim([-1.45, 0.25]);
    xlim([-1.3, 3.1]);     

    axes; a(1) = qPlotDR( ar, arFindData( 'app_hep_2014_04_22_qPCR_140604', 'input_ruxolitinib', 500, 'input_stattic', 0), 'input_il6', 1440, rux_color, 'SOCS3_qpcr', -1 ); hold on;
    a(2) = qPlotDR( ar, arFindData( 'app_hep_2014_04_22_qPCR_140604', 'input_ruxolitinib', 0, 'input_stattic', 0), 'input_il6', 1440, control_color, 'SOCS3_qpcr', -1 );
    qLabels( 'Socs3 mRNA (24 h)', 'IL6 (ng/mL)', 'Concentration (log_{10} a.u.)', [0.1, 1, 10, 100, 1000], [-1.5, -1, -0.5, 0], 3, 1, [1, 0] ); ylim([-1.51, 0.1]);
    xlim([-1.3, 3.1]);  
    printfig( 'Core Model C' );

    % Core model figures pt 4
    clear a;
    figure('units','normalized','outerposition',[0 0 1 1]);
    axes; a(1) = qPlotDR( ar, arFindData( '2013_11_04', 'input_stattic', 0, 'input_ruxolitinib', 0), 'input_il6', [], control_color, 'pSTAT3_lumi', -1, 1 ); hold on;
    a(2) = qPlotDR( ar, arFindData( '2013_11_04', 'input_stattic', 60, 'input_ruxolitinib', 0), 'input_il6', [], sta_color, 'pSTAT3_lumi', -1 );
    a(3) = qPlotDR( ar, arFindData( '2013_11_04', 'input_stattic', 0, 'input_ruxolitinib', 500), 'input_il6', [], rux_color, 'pSTAT3_lumi', 0 );
    qLabels( 'phospho-STAT3 (lumi)', 'Inhibitor (\mu M)', 'Concentration (log_{10} a.u.)', [0.1, 1, 10, 100, 1000], [-4, -2, -1.5, -1, -0.5, 0], 2, 1, [1, 0] ); xlim([-1.3, 3]); ylim([-1.5 .4]);
    qLegend( a, {'DMSO', 'Stattic (60 \mu M)', 'Ruxolitinib (500 nM)'} );
    printfig( 'Core Model D' );

    % Core model figures pt 5
    clear a;
    figure('units','normalized','outerposition',[0 0 1 1]);
    axes; a(1) = qPlotTC( ar, arFindData( '2011_09_06', 'input_il6', 0.5), [0 0 0]/255, 'pSTAT3_wb' ); hold on;
    a(2) = qPlotTC( ar, arFindData( '2011_09_06', 'input_il6', 5), [0 144 54]/255, 'pSTAT3_wb' );
    a(3) = qPlotTC( ar, arFindData( '2011_09_06', 'input_il6', 50), [0 127 255]/255, 'pSTAT3_wb' );
    a(4) = qPlotTC( ar, arFindData( '2011_09_06', 'input_il6', 500), [226 0 26]/255, 'pSTAT3_wb' );
    qLabels( 'phospho-STAT3 (wb)', 'Time (min)', 'Concentration (log_{10} a.u.)', [0 10 20 30 40 50 60 70], [-1.5 -1.0, -0.5, 0], 1, 1, [0, 0] );
    qLegend( a, {'IL6 (0.5 ng/mL)', 'IL6 (5 ng/mL)', 'IL6 (50 ng/mL)', 'IL6 (500 ng/mL)'} ); ylim([-1.9, 0.5]); xlim([-4 74]);
    printfig( 'Core Model E' );

end
