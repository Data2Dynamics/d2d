function xiaoPlot(png, printFigs)

    PlotSettings;
    global ar;

    if ~exist('png', 'var')
        png = 0;
    end
    if ~exist('printFigs', 'var')
        printFigs = 0;
    end    
    
    WT_color        = [0 127 255]./255;
    gray_color      = [127 127 127]/255;
    rux_color       = [255 0 255]./255;

    data{1} = 'xiaoyun_validation_nucSTAT3_validation_20150428_nExpID1';
    data{2} = 'xiaoyun_validation_nucSTAT3_validation_20150428_nExpID2';
    data{3} = 'xiaoyun_validation_nucSTAT3_validation_20150508_nExpID1';
    data{4} = 'xiaoyun_validation_nucSTAT3_validation_20150508_nExpID2';
    data{5} = 'xiaoyun_validation_nucSTAT3_validation_20150528_nExpID1';
    data{6} = 'xiaoyun_validation_nucSTAT3_validation_20150528_nExpID2';

    name{1} = '20150428 rep 1';
    name{2} = '20150428 rep 2';
    name{3} = '20150508 rep 1';
    name{4} = '20150508 rep 2';
    name{5} = '20150528 rep 1';
    name{6} = '20150528 rep 2';

    tp = [1200, 1215, 1230, 1245, 1260, 1275, 1290, 1305, 1320, 1335, 1350, 1395, 1410, 1425, 1440];

    if ( printFigs )
        switch( png )
            case 1
                directory = 'Validation_png';
                printFigure = @(filename)print('-dpng', '-cmyk', '-loose', '-r300', '-painters', sprintf('figures/%s/%s.png', directory, filename) );
            case 0
                directory = 'Validation';
                printFigure = @(filename)print('-depsc2', '-cmyk', '-loose', '-r300', '-painters', sprintf('figures/%s/%s.eps', directory, filename) );
            case 2
                directory = 'Validation';
                printFigure = @(filename)printTikz(sprintf('figures/%s/%s.pdf', directory, filename));
        end
        mkdir( sprintf( 'figures/%s', directory ) );
    else
        printFigure = @(filename)disp('No file export!');
    end

    obsName         =   'nucSTAT3';
    meanSimu        =   @(id)(ar.model.data(id).yExpSimu(ismember(ar.model.data(id).tExp, tp), ismember(ar.model.data(id).y, obsName) ) ).';
    meanExp         =   @(id)(ar.model.data(id).yExp(ismember(ar.model.data(id).tExp, tp), ismember(ar.model.data(id).y, obsName) ) ).';
    meanStd         =   @(id)(ar.model.data(id).ystdExpSimu(ismember(ar.model.data(id).tExp, tp), ismember(ar.model.data(id).y, obsName) ) ).';
    meanMean        =   @(fun, ids)mean( cell2mat( arrayfun(fun, ids, 'UniformOutput', false) ) );

    figure('units','normalized','outerposition',[0 0 1 1]);

    for replicate = 1 : 6
        dose{1}{1} = arFindData(ar, data(replicate), 'input_il6', 0,    'input_rux1', 0,    'isTriple', 0 );
        dose{1}{2} = arFindData(ar, data(replicate), 'input_il6', 0,    'input_rux1', 500,  'isTriple', 0 );
        dose{1}{3} = arFindData(ar, data(replicate), 'input_il6', 0,    'input_rux1', 500,  'isTriple', 1 );

        dose{2}{1} = arFindData(ar, data(replicate), 'input_il6', 7.5,  'input_rux1', 0,    'isTriple', 0 );
        dose{2}{2} = arFindData(ar, data(replicate), 'input_il6', 7.5,  'input_rux1', 500,  'isTriple', 0 );
        dose{2}{3} = arFindData(ar, data(replicate), 'input_il6', 7.5,  'input_rux1', 500,  'isTriple', 1 );

        dose{3}{1} = arFindData(ar, data(replicate), 'input_il6', 100,  'input_rux1', 0,    'isTriple', 0 );
        dose{3}{2} = arFindData(ar, data(replicate), 'input_il6', 100,  'input_rux1', 500,  'isTriple', 0 );
        dose{3}{3} = arFindData(ar, data(replicate), 'input_il6', 100,  'input_rux1', 500,  'isTriple', 1 );

        for d = 1 : length( dose )
            q = dose{d};
            y(d,:)              = [ meanMean( meanSimu, q{1} ), meanMean( meanSimu, q{2} ), meanMean( meanSimu, q{3} ) ];   %#ok
            yExp(d,:)           = [ meanMean( meanExp, q{1} ), meanMean( meanExp, q{2} ), meanMean( meanExp, q{3} ) ];      %#ok
            ystdExpSimu(d,:)    = [ meanMean( meanStd, q{1} ), meanMean( meanStd, q{2} ), meanMean( meanStd, q{3} ) ];      %#ok
        end       

        qBar( {'0', '7.5', '100'}, y, yExp, ystdExpSimu, {gray_color, rux_color, WT_color}, 1.7 );
        grid off;
        ylim([1.7, 2.3]);
        xlabel('IL6 dose (ng/mL)');
        qLabels( 'ntSTAT3', 'IL6 dose (ng/mL)', 'log_{10}( mean nuclear STAT3 ) [a.u.]', [], [], replicate-2*floor(replicate/3), floor(replicate/3),  [0,0] )
        title(name{replicate});
    end

	 printFigure( 'nuclearSTAT3' );
end
