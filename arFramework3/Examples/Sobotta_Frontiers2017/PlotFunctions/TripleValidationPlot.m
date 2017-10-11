function TripleValidationPlot( png, printFigs )

	close all;
    PlotSettings;
	global ar;

	WT_color        = [0 127 255]./255;
	gray_color      = [127 127 127]./255;
	rux_color       = [255 0 255]./255;

	if ~exist( 'png', 'var' )
		 png = 0;
    end
    if ~exist( 'printFigs', 'var' )
        printFigs = 0;
    end
    
	if ( printFigs )
        switch( png )
             case 0
    		     directory = 'Validation';
    		     printfig = @(filename)print('-depsc2', '-loose', '-cmyk', '-r300', '-painters', sprintf('figures/%s/%s.eps', directory, filename) );
             case 1
                directory = 'Validation_png';
                printfig = @(filename)print('-dpng',   '-loose', '-cmyk', '-r300', '-painters', sprintf('figures/%s/%s.png', directory, filename) );
             case 2
                directory = 'Validation';
                printfig = @(filename)printTikz(sprintf('figures/%s/%s.pdf', directory, filename));
        end
        mkdir figures;      
        dir = sprintf( 'figures/%s', directory );
        mkdir( dir );
	else
		 printfig = @(filename)disp('No file export!');
	end



	% Triple validation experiment
     for replicate = 1 : 3
         close all;
         figure('units','normalized','outerposition',[0 0 1 1]);

         dose{1}(1) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 0,    'input_rux1', 0,    'input_rux2', 0);
         dose{1}(2) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 0,    'input_rux1', 500,  'input_rux2', 0);
         dose{1}(3) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 0,    'input_rux1', 500,  'input_rux2', 190.9888);

         dose{2}(1) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 7.5,  'input_rux1', 0,    'input_rux2', 0);
         dose{2}(2) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 7.5,  'input_rux1', 500,  'input_rux2', 0);
         dose{2}(3) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 7.5,  'input_rux1', 500,  'input_rux2', 190.9888);

         dose{3}(1) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 100,  'input_rux1', 0,    'input_rux2', 0);
         dose{3}(2) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 100,  'input_rux1', 500,  'input_rux2', 0);
         dose{3}(3) = arFindData(ar, sprintf('2015_04_13_triple_treatment_replicate%d', replicate), 'input_il6', 100,  'input_rux1', 500,  'input_rux2', 190.9888);

         APPnames = { '{\itSocs3} mRNA', '{\itFgg} mRNA', '{\itHamp} mRNA', '{\itIl33} mRNA', '{\itApcs} mRNA', '{\itHp} mRNA', '{\itHpx} mRNA' };
         APPs = { 'SOCS3_qpcr', 'FGG_qpcr', 'HAMP_qpcr', 'IL33_qpcr', 'APCS_qpcr', 'HP_qpcr', 'HPX_qpcr' };

         for a = 1 : length( APPs )
            for d = 1 : length( dose )
                q = dose{d};
                y(d,:) = [ ar.model.data(q(1)).yExpSimu( strcmp( ar.model.data(q(1)).y, APPs{a} ) ), ...
                           ar.model.data(q(2)).yExpSimu( strcmp( ar.model.data(q(2)).y, APPs{a} ) ), ...
                           ar.model.data(q(3)).yExpSimu( strcmp( ar.model.data(q(3)).y, APPs{a} ) ) ]; %#ok
                yExp(d, :) = [ ...   
                                ar.model.data(q(1)).yExp( strcmp( ar.model.data(q(1)).y, APPs{a} ) ), ...
                                ar.model.data(q(2)).yExp( strcmp( ar.model.data(q(2)).y, APPs{a} ) ), ...
                                ar.model.data(q(3)).yExp( strcmp( ar.model.data(q(3)).y, APPs{a} ) ) ]; %#ok
                ystdExpSimu(d, :) = [ ...
                                ar.model.data(q(1)).ystdExpSimu( strcmp( ar.model.data(q(1)).y, APPs{a} ) ), ...
                                ar.model.data(q(1)).ystdExpSimu( strcmp( ar.model.data(q(2)).y, APPs{a} ) ), ...
                                ar.model.data(q(1)).ystdExpSimu( strcmp( ar.model.data(q(3)).y, APPs{a} ) ) ]; %#ok
            end

            qBar( {'0', '7.5', '100'}, y, yExp, ystdExpSimu, {gray_color, rux_color, WT_color}, [] );
            title(APPnames{a});
            xlabel('IL6 dose (ng/mL)');
            qLabels( APPnames{a}, 'IL6 dose (ng/mL)', 'log_{10}( Concentration ) [a.u.]', [], [], a-2*floor(a/3), floor(a/3), [0,0] );
            grid off;
         end

         printfig( sprintf( 'Triple treatment replicate %d',  replicate ) );
     end
end

