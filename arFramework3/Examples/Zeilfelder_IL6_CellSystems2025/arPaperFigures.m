% This function reproduces the figures from the paper.
%
%  arPaperFigures(CORE, SMAD, PHH, HAMP, (fn), (annotateString))
%
% Load the respective models and plot with that flag active.
%
% Example:
%   SetupControlIL6;
%   loadParameters;
%   arPaperFigures(1, 0, 0, 0);
% 

function arPaperFigures(CORE, SMAD, PHH, HAMP, fn, annotateString)
    global ar;
        
    if ( nargin < 1 )
        disp('function arPaperFigures(CORE, SMAD, PHH, HAMP, fn, annotateString)');
        return;
    end
    
    if ~exist('fn', 'var')
        fn = [];
    end
    
    % Any annotations requested?
    if ~exist('annotateString', 'var')
        annotate = @()1;
    else
        annotate = @()annotation('textbox', [0, 0, 1, 1], 'string', annotateString);
    end
    
    % Make sure simulation is up to date
    arCalcMerit(false, ar.p(ar.qFit==1), true);
    arSimu(false, true, true);
    
    % How to display model errors
    errmode = 1;
    
    % Font size
    fs = 12;
    
    %% HAMP
    if ( HAMP )
             
        condis    = {   {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '500',  'input_il6',  '0',   'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '500',  'input_il6',  '10',  'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '500',  'input_il6',  '10',  'input_apap', '0',   'input_noggin', 500, 'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '10',  'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '10',  'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 0  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '10',  'input_noggin', 500, 'input_bmp2', 0, 'input_phh', 0  }, ...
                       ...
                        {'input_dcf', '0',    'input_il6',  '0',  'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 25, 'input_phh', 0  }, ...
                        ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '500',  'input_il6',  '0',   'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '500',  'input_il6',  '10',  'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '500',  'input_il6',  '10',  'input_apap', '0',   'input_noggin', 500, 'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '10',  'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',   'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '10',  'input_noggin', 0,   'input_bmp2', 0, 'input_phh', 1  }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '10',  'input_noggin', 500, 'input_bmp2', 0, 'input_phh', 1  }, ...
                        };

        colors    = { [.5, .5, .5], [1, .5, .5], [0, 0, 0], [1, 0, 0], [.5, 0, .6], [.5, .5, .5], [.3, .5, .3], [0, 0, 0], [0, .5, 0], [0, .2, .6], [1, .5, 1] };
       
        hasNoggin = 1;
        
        % STATs
            figure('Position', [0,0,1400,1000])
            subplot(3,2,1);
            h = plotTogether( 1, 'tpSTAT3_wb_scaled', condis(1:4+hasNoggin), colors(1:4+hasNoggin), errmode );
            ylabel( 'total pSTAT3 [au]' );
            subplot(3,2,2);
            h2 = plotTogether( 1, 'tpSTAT3_wb_scaled', condis(6:9+hasNoggin), colors(6:9+hasNoggin), errmode );
            ylabel( 'cytoplasmic pSTAT3 [au]' );
            
            subplot(3,2,3);
            plotTogether( 1, 'cpSTAT3_wb_scaled', condis(1:4+hasNoggin), colors(1:4+hasNoggin), errmode );
            ylabel( 'cytoplasmic pSTAT3 [au]' );
            subplot(3,2,4);
            plotTogether( 1, 'cpSTAT3_wb_scaled', condis(6:9+hasNoggin), colors(6:9+hasNoggin), errmode );
            ylabel( 'cytoplasmic pSTAT3 [au]' );

            subplot(3,2,5);
            plotTogether( 1, 'npSTAT3_wb_scaled', condis(1:4+hasNoggin), colors(1:4+hasNoggin), errmode );
            ylabel( 'nuclear pSTAT3 [au]' );
            subplot(3,2,6);
            plotTogether( 1, 'npSTAT3_wb_scaled', condis(6:9+hasNoggin), colors(6:9+hasNoggin), errmode );
            ylabel( 'nuclear pSTAT3 [au]' );
        
            subplot(3,2,1);
            str = { 'Control', 'DCF', 'IL6', 'DCF + IL6', 'DCF + IL6 + Noggin', 'BMP2' };
            legend( [h{:}], str{1:5+hasNoggin} );
            
            subplot(3,2,2);
            str = { 'Control', 'APAP', 'IL6', 'APAP + IL6', 'APAP + IL6 + Noggin', 'BMP2' };
            legend( [h2{:}], str{1:5+hasNoggin} );        
            drawnow;
            
        % SMAD, SOCS3 and HAMP
            figure('Position', [0,0,1400,1000]); drawnow;
            subplot(3,3,1);
            plotTogether( 1, 'pSMAD_wb_scaled', condis([1:4+hasNoggin, 11]), colors([1:4+hasNoggin, 11]), errmode );
            ylabel( 'pSMAD [au]' );
            
            subplot(3,3,2);
            plotTogether( 1, 'pSMAD_wb_scaled', condis([6:9+hasNoggin, 11]), colors([6:9+hasNoggin,11]), errmode );
            ylabel( 'pSMAD [au]' );
            
            condisA = [ condis(1:4+hasNoggin), condis(11) ];
            colorsA = [ colors(1:4+hasNoggin), colors(11) ];
            
            condisB = [ condis(6:9+hasNoggin), condis(11) ];
            colorsB = [ colors(6:9+hasNoggin), colors(11) ];
            
            subplot(3,3,4);
            plotTogether( 1, 'SOCS3_qpcr_scaled', condisA, colorsA, errmode );
            ylabel( 'log_{10} SOCS3 mRNA [au]' );
            xlabel( 'Time [min]' );
            ylim([-2,2]);

            subplot(3,3,5);
            plotTogether( 1, 'SOCS3_qpcr_scaled', condisB, colorsB, errmode );
            ylabel( 'log_{10} SOCS3 mRNA [au]' );
            xlabel( 'Time [min]' );  
            ylim([-2,2]);
        
            subplot(3,3,7);
            h = plotTogether( 1, 'HAMP_qpcr_scaled', condisA, colorsA, errmode );
            ylabel( 'log_{10} HAMP mRNA [au]' );
            xlabel( 'Time [min]' );
            ylim([-1, 2]);

            subplot(3,3,8);
            h2 = plotTogether( 1, 'HAMP_qpcr_scaled', condisB, colorsB, errmode );
            ylabel( 'log_{10} HAMP mRNA [au]' );
            xlabel( 'Time [min]' );    
            ylim([-1, 2]);
        
            subplot(3,3,7);
            str = { 'Control', 'DCF', 'IL6', 'DCF + IL6', 'DCF + IL6 + Noggin', 'BMP2' };
            legend( [h{:}], str );
            
            subplot(3,3,8);
            str = { 'Control', 'APAP', 'IL6', 'APAP + IL6', 'APAP + IL6 + Noggin', 'BMP2' };
            legend( [h2{:}], str );
            drawnow;
    end
            
    %% IL-6 module plots - only control, no APAP/DCF
    if ( CORE )
        condis    = {   {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '500',  'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '500',  'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '10', 'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '10', 'input_noggin', 0, 'input_phh', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0, 'input_phh', 1 }, ...
                        {'input_dcf', '500',  'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0, 'input_phh', 1 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0, 'input_phh', 1 }, ...
                        {'input_dcf', '500',  'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0, 'input_phh', 1 }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0, 'input_phh', 1 }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '10', 'input_noggin', 0, 'input_phh', 1 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0, 'input_phh', 1 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '10', 'input_noggin', 0, 'input_phh', 1 }, ...
                        };

        colors    = { [0, 0, 0], [1, .5, .5], [0, 0, 0], [1, 0, 0], [0, 0, 0], [1, .5, .5], [0,0,0], [0, .5, 0] };

        vc = str2num(ar.model.fp{strcmp(ar.model.p, 'vol_cyt')}); %#OK str2double doesn't handle brackets
        vn = str2num(ar.model.fp{strcmp(ar.model.p, 'vol_nuc')}); %#OK
        mc = arGetPars('mcyto', 0);
        mn = arGetPars('mnuc', 0) ;
            
        figure('Position', [0, 0, 800, 1100]);

        subplot(4,2,1);
        plotTogether( 1, 'cpSTAT3_wb_scaled', condis(1:4), colors(1:4), errmode, @(x)1e-3*x/(vc)  );
        ylabel( 'cytoplasmic pSTAT3 [\muM]', 'FontSize', fs );
        ylim([0, .8]);
        subplot(4,2,2);
        plotTogether( 1, 'cpSTAT3_wb_scaled', condis(5:8), colors(5:8), errmode, @(x)1e-3*x/(vc) );
        ylabel( 'cytoplasmic pSTAT3 [\muM]', 'FontSize', fs );
        ylim([0, .8]);

        subplot(4,2,3);
        plotTogether( 1, 'npSTAT3_wb_scaled', condis(1:4), colors(1:4), errmode, @(x)1e-3*x/(vn*(mc/mn)) );
        ylabel( 'nuclear pSTAT3 [\muM]', 'FontSize', fs );
        ylim([0,1])
        subplot(4,2,4);
        plotTogether( 1, 'npSTAT3_wb_scaled', condis(5:8), colors(5:8), errmode, @(x)1e-3*x/(vn*(mc/mn)) );
        ylabel( 'nuclear pSTAT3 [\muM]', 'FontSize', fs );
        ylim([0,1])

        subplot(4,2,5);
        plotTogether( 1, 'SOCS3_qpcr_scaled', condis(1:4), colors(1:4), errmode, @(x)x );
        ylabel( 'log_{10} SOCS3 mRNA [au]', 'FontSize', fs );
        xlabel( 'Time [min]', 'FontSize', fs );
        ylim([-2,2.5]);

        subplot(4,2,6);
        plotTogether( 1, 'SOCS3_qpcr_scaled', condis(5:8), colors(5:8), errmode, @(x)x );
        ylabel( 'log_{10} SOCS3 mRNA [au]', 'FontSize', fs );
        xlabel( 'Time [min]', 'FontSize', fs );
        ylim([-2,2.5]);
        
        subplot(4,2,7);
        plotTogether( 1, 'cSOCS3_wb', condis(1:4), colors(1:4), errmode );
        ylabel( 'SOCS3 [nM]', 'FontSize', fs );
        ylim([0, 30]);
        subplot(4,2,8);
        plotTogether( 1, 'cSOCS3_wb', condis(5:8), colors(5:8), errmode );
        ylabel( 'SOCS3 [nM]', 'FontSize', fs );
        ylim([0, 30]);

        annotate();
        if ~isempty(fn)
            saveas(gcf, sprintf('Figures/%s_HepG2_SOCS3RNA_STAT3.fig', fn) );
        end
    end
    
    if ( SMAD )
        condis    = {   {'input_dcf', '0',    'input_apap', '0',   'input_bmp2', 0,     'input_noggin', 0  }, ...
                        {'input_dcf', '500',  'input_apap', '0',   'input_bmp2', 0,     'input_noggin', 0  }, ...
                        {'input_dcf', '0',    'input_apap', '0',   'input_bmp2', 25,    'input_noggin', 0 }, ...
                        {'input_dcf', '0',    'input_apap', '10',  'input_bmp2', 0,     'input_noggin', 0  }, ...
                        {'input_dcf', '0',    'input_apap', '0',   'input_bmp2', 25,    'input_noggin', 500 }, ...
                        };
                    
        colors    = { [0, 0, 0], [1, .5, .5], [1, .5, 1], [.5, .5, 1], [1, .3, 0] };

        figure('Position', [0,0,1400,1100]);
        subplot(2,3,1);
        h = plotTogether( 1, 'pSMAD_wb_scaled', condis(1:5), colors(1:5), errmode );
        ylabel( 'pSMAD [au]' );
        
        subplot(2,3,2);
        plotTogether( 1, 'SMAD6_qpcr_scaled', condis(1:5), colors(1:5), errmode );
        ylabel( 'log_{10} SMAD6 qPCR [au]' );
        drawLegend( h, {'Control', 'DCF', 'BMP2', 'APAP', 'BMP2 + Noggin' } );      
        
        subplot(2,3,3);
        offs = max(10.^ar.p(arFindPar('offset_ID3')))/mean(10.^ar.p(arFindPar('scale_ID3')));
        plotTogether( 1, 'ID3_qpcr_scaled', condis(1:5), colors(1:5), errmode, @(x)x, [], [], offs );
        ylabel( 'log_{10} ID3 qPCR [au]' );
        drawLegend( h, {'Control', 'DCF', 'BMP2', 'APAP', 'BMP2 + Noggin' } );
        %ylim([-2.5,1.5])

        subplot(2,3,4);
        plotTogether( 1, 'pSMAD_DR_wb_scaled', {{'input_bmp2', 0, 'input_apap', 0, 'input_noggin', 0}},  colors([5]), errmode, [], 60, 'input_dcf' );
        xlabel( 'DCF' );
        ylabel( 'pSMAD 1/5/8' );
        ax1 = gca;
        ax1.XColor = colors{5};
        ax1_pos = ax1.Position;
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right', ...
            'XColor', colors{3}, ...
            'Color','none');
        
        hold on; plotTogether( 1, 'pSMAD_DR_wb_scaled', {{'input_dcf', 0, 'input_noggin', 0, 'input_apap', 0}},  colors([3]), errmode, [], 60, 'input_bmp2' );
        xlabel( 'BMP2' );
        ylabel( 'pSMAD 1/5/8 [au]' );
        title( 't = 60');
        linkaxes([ax1, ax2], 'y')
        
        subplot(2,3,5);
        plotTogether( 1, 'pSMAD_DR_wb_scaled', {{'input_bmp2', 0, 'input_dcf', 0, 'input_noggin', 0}},  colors([4]), errmode, [], 60, 'input_apap' );
        xlabel( 'APAP' );
        ylabel( 'pSMAD 1/5/8 [au]' );
        ax1 = gca;
        ax1.XColor = colors{4};
        ax1_pos = ax1.Position;
        ax2 = axes('Position',ax1_pos,...
            'XAxisLocation','top',...
            'YAxisLocation','right', ...
            'XColor', colors{3}, ...
            'Color','none');
        
        hold on; plotTogether( 1, 'pSMAD_DR_wb_scaled', {{'input_dcf', 0, 'input_noggin', 0, 'input_apap', 0}},  colors([3]), errmode, [], 60, 'input_bmp2' );
        xlabel( 'BMP2' );
        ylabel( 'pSMAD 1/5/8' );
        title( 't = 60');
        linkaxes([ax1, ax2], 'y')

        annotate();
        if ~isempty(fn)
            saveas(gcf, sprintf('Figures/%s_figure.fig', fn) );
        end
    end
    

    if ( PHH )
        condis    = {   {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0 }, ...
                        {'input_dcf', '500',  'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0 }, ...
                        {'input_dcf', '500',  'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '0',  'input_noggin', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '0',   'input_apap', '10', 'input_noggin', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '0',  'input_noggin', 0 }, ...
                        {'input_dcf', '0',    'input_il6',  '10',  'input_apap', '10', 'input_noggin', 0 }, ...
                        };

        colors    = { [0, 0, 0], [1, .5, .5], [0, 0, 0], [1, 0, 0], [0, 0, 0], [.3, .5, .3], [0, 0, 0], [0, .5, 0] };
        
        vc = str2num(ar.model.fp{strcmp(ar.model.p, 'vol_cyt')});
        vn = str2num(ar.model.fp{strcmp(ar.model.p, 'vol_nuc')});
        mc = arGetPars('mcyto', 0);
        mn = arGetPars('mnuc', 0) ;

        figure('Position', [0, 0, 800, 1100]);
        subplot(3,2,1);
        plotTogether( 1, 'tpSTAT3_wb_scaled', condis([1:4]), colors(1:4), errmode, @(x)1e-3*x/((vc+vn)*(mc/(mn+mc))) );
        ylabel( 'total pSTAT3 [\muM]', 'FontSize', fs );
        
        subplot(3,2,2);
        plotTogether( 1, 'tpSTAT3_wb_scaled', condis([5:8]), colors(5:8), errmode, @(x)1e-3*x/((vc+vn)*(mc/(mn+mc))) );
        ylabel( 'total pSTAT3 [\muM]', 'FontSize', fs );

        subplot(3,2,3);
        plotTogether( 1, 'SOCS3_qpcr_scaled', condis([1:4]), colors(1:4), errmode, @(x)x );
        ylabel( 'log_{10} SOCS3 mRNA [au]', 'FontSize', fs );
        xlabel( 'Time [min]', 'FontSize', fs );
        ylim([-3, 0])

        subplot(3,2,4);
        plotTogether( 1, 'SOCS3_qpcr_scaled', condis([5:8]), colors(5:8), errmode, @(x)x );
        ylabel( 'log_{10} SOCS3 mRNA [au]', 'FontSize', fs );
        xlabel( 'Time [min]', 'FontSize', fs );
        ylim([-3, 0])
        
        subplot(3,2,5);
        % Offset to scale the other variables to
        refOffset = max(10.^ar.p(arFindPar('offset_hamp')));
        refOffset = refOffset / mean(10.^ar.p(arFindPar('scale_hamp')));
        
        plotTogether( 1, 'HAMP_qpcr_scaled', condis([1:4]), colors(1:4), errmode, @(x)x, [], [], refOffset );
        ylabel( 'log_{10} HAMP mRNA [au]', 'FontSize', fs );
        xlabel( 'Time [min]', 'FontSize', fs );
        ylim([-3, -1.8]);

        subplot(3,2,6);
        plotTogether( 1, 'HAMP_qpcr_scaled', condis([5:8]), colors(5:8), errmode, @(x)x, [], [], refOffset );
        ylabel( 'log_{10} HAMP mRNA [au]', 'FontSize', fs );
        xlabel( 'Time [min]', 'FontSize', fs );
        ylim([-3, -1.8]);

        figure;
        subplot(1,2,1);
        plotTogether( 1, 'stat3_abs', condis([1:4]), colors(1:4), errmode );
        ylabel( 'STAT3 [zmol/l]', 'FontSize', fs );
        ylim([0,1000]);

        subplot(1,2,2);
        plotTogether( 1, 'socs3_abs', condis([1:4]), colors(1:4), errmode );
        ylabel( 'SOCS3 [zmol/l]', 'FontSize', fs );     
        ylim([0, 3]);

        annotate();
        if ~isempty(fn)
            saveas(gcf, sprintf('Figures/%s_figure.fig', fn) );
        end
    end
end

function drawLegend( h, items )
    r = find( ~cellfun( @isempty, h ) );
    items = items(r);
    h = [h{r}];
    legend( h, items );
end

% Function to plot some data
% m             - Model
% observable    - Observable
% condi         - Cell array of conditions
% color         - Color [R,G,B]
% errmode       - Mean+SEM or replicates
% trafo         - Output transformation
% timePt        - Needed in case of dose responses (time point at which the dose response was taken)
% responseVar   - Dose response variable
% refOffset     - In data with offset, use this as target offset for
%                 plotting. Model errors are propagated using Gaussian
%                 Error Propagation.
function h = plotTogether( m, observable, condis, colors, errmode, trafo, timePt, responseVar, refOffset )
    
    % Storage to keep minimum and maximum value
    xmi = 1e9;
    xma = -1e9;
    ymi = 1e9;
    yma = -1e9;

    for j = 1 : numel( condis )
        if ( nargin > 6 ) && ~isempty( timePt ) && ~isempty( responseVar )
            [ ts, tfs, yExps, ySims, yExpStds, ySimStds, scales, offsets, scalesSim, offsetsSim, qLogExps, qLogSims ] = grabDR( m, observable, timePt, responseVar, condis{j} );
        else
            [ ts, tfs, yExps, ySims, yExpStds, ySimStds, scales, offsets, scalesSim, offsetsSim, qLogExps, qLogSims ] = grabData( m, observable, condis{j} );
        end
        
        % Non-logarithmic observables
        ySims(~qLogSims)        = ( ySims(~qLogSims) - offsetsSim(~qLogSims) ) ./ scalesSim(~qLogSims);
        ySimStds(~qLogSims)     = ySimStds(~qLogSims) ./ scalesSim(~qLogSims);
        yExps(~qLogExps)        = ( yExps(~qLogExps) - offsets(~qLogExps) ) ./ scales(~qLogExps);
        yExpStds(~qLogExps)     = yExpStds(~qLogExps) ./ scales(~qLogExps);
        
        % Logarithmic observables
        lID = find(qLogSims);
        lID2 = find(qLogExps);
        if ( nansum( offsetsSim ) == 0 )
            refOffsetSim = 0;
            refOffsetExp = 0;
            refOffset = 0;
        end
        
        % Gaussian error propagation for log transform
        x = ySims(qLogSims); gamma = refOffset; s = 10.^scalesSim(qLogSims); beta = offsetsSim(qLogSims);
        ySimStds(qLogSims)      = ySimStds(qLogSims) .* ( 10.^x ./ ( 10.^x + gamma.*s - beta ) );
        
        % Gaussian error propagation for log transform
        x = yExps(qLogExps); gamma = refOffset; s = 10.^scales(qLogExps); beta = offsets(qLogExps);
        yExpStds(qLogExps)      = yExpStds(qLogExps) .* ( 10.^x ./ ( 10.^x + gamma.*s - beta ) );
        
        % Undoing the observation function for logarithmically transformed
        % data. Note that this is only approximate if an offset is present.
        ySims(qLogSims)         = log10( ( 10.^ySims(qLogSims) - offsetsSim(qLogSims) ) ./ 10.^scalesSim(qLogSims) + refOffset );
        yExps(qLogExps)         = log10( ( 10.^yExps(qLogExps) - offsets(qLogExps) ) ./ 10.^scales(qLogExps) + refOffset );
        
        if ( errmode == 1 )
            [ ts, yExps, yExpStds ] = convertMeanStd(ts, yExps, yExpStds);
        end

        sU                      = ySims + ySimStds;
        sL                      = ySims - ySimStds;        
        eU                      = yExps + yExpStds;
        eL                      = yExps - yExpStds;
        
        [tfs, I]                = sort( tfs );
        tsim{j}                 = tfs;
        sim{j}                  = ySims(I);
        simL{j}                 = sL(I);
        simU{j}                 = sU(I);
        
        texp{j}                 = ts;
        exp{j}                  = yExps;
        expL{j}                 = eL;
        expU{j}                 = eU;
    end
    
    hold on;
        
    % Draw patches
    plotErrors = 0;
    if ( plotErrors )
        for j = 1 : numel( condis ) 
            nL = numel( simL{j} );
            if ( nL > 250 )
                idx = round([1:(nL-1)/250:nL]);               
                simL{j} = simL{j}(idx) + 0;
                simU{j} = simU{j}(idx) + 0;
                tsim{j} = tsim{j}(idx) + 0;
                sim{j} = sim{j}(idx) + 0;
            end
                
            time = [tsim{j}; flipud(tsim{j})];
            hc = patch( [tsim{j}; flipud(tsim{j})], [ simL{j}; flipud(simU{j}) ], zeros(size(time)), 'EdgeColor', lighten(colors{j}), 'FaceColor', lighten(colors{j}) );
        end
    end
    for j = 1 : numel( condis )
        mockpatch = patch( [0,0,0,0], [0,0,0,0], [0,0,0,0], 'EdgeColor', colors{j}, 'FaceColor', colors{j} );
        h{j} = mockpatch;
    end
    
    lw = 1.1;
    if ~exist( 'trafo', 'var') || isempty( trafo )
        trafo = @(x)x;
    end
    for j = 1 : numel( condis )
        plot( tsim{j}, trafo(sim{j}), '-', 'Color', colors{j}, 'LineWidth', lw )
        if ( errmode == 0 )
            if ~isreal(trafo(exp{j}))
                error('Fatal error: Result is not real');
            end
            plot( texp{j}, trafo(exp{j}), 'o', 'Color', colors{j}, 'MarkerSize', 2.5 );
        else
            errorbar( texp{j}, trafo(exp{j}), abs(trafo(expL{j})-trafo(exp{j})), abs(trafo(expU{j})-trafo(exp{j})), 'o', 'Color', colors{j}, 'MarkerSize', 2, 'LineWidth', lw );
        end
        xmi = min( [ xmi, min( texp{j} ) ] );
        xma = max( [ xma, max( texp{j} ) ] );
        ymi = min( [ ymi, min( trafo(expL{j}) ), min(trafo(simL{j})) ] );
        yma = max( [ yma, max( trafo(expU{j}) ), max(trafo(simU{j})) ] );
    end
    sz = xma - xmi;
    
    try
        xlim([xmi - 0.1*sz,xma + 0.1*sz]);
    catch
        warning( 'Invalid limit' );
    end
    
    sz = abs(yma - ymi);
    try
        ylim([ymi - 0.1*sz,yma + 0.1*sz]);
    catch
        warning( 'Invalid limit' );
    end
    
end

function c = lighten( c )
    c = c*.2 + .8;
    c(c>1.0) = 1.0;
end

function [ tList, yAvg, yStdAvg ] = convertMeanStd(t, y, yStd)
    tList   = t;
    yAvg    = [];
    yStdAvg = [];
    
    [tList, ~, I] = unique(t);
    for j = 1 : numel( tList )
        yAvg(j) = nanmean( y(I==j) ); %#ok<AGROW>
        stds = yStd(I==j);
        N = sum( ~isnan( stds ) );
        yStdAvg(j) = sqrt( nansum( stds.^2 ) / (N*N) ); %#ok<AGROW>
    end
end

function printCondition(cond)
    s = '';
    for i = 1 : numel( cond )
        s = sprintf( '%s %s=%s ', s, cond(i).parameter, cond(i).value );
    end
    fprintf( '%s\n', s );
end

% Grab dose responses
function [ ts, tfs, yExps, ySims, yExpStds, ySimStds, scales, offsets, scalesSim, offsetsSim, qLogExps, qLogSims ] = grabDR( m, observable, timePt, responseVar, condis )
    global ar;
    
    debugInfo = 0;
    ds = setdiff( setdiff( arFindData('_', condis{:} ), arFindData('decay') ), arFindData('steady' ) );
    
    ts          = [];
    tfs         = [];
    ySims       = [];
    yExps       = [];
    scales      = [];
    offsets     = [];
    scalesSim   = [];
    offsetsSim  = [];
    yExpStds    = [];
    ySimStds    = [];
    qLogExps    = [];
    qLogSims    = [];
    
    for jd = 1 : numel( ds )
        id = ds(jd);
        obsidx  = find( strcmp( ar.model(m).data(id).y, observable ) );
               
        if ~isempty( obsidx ) && ~isempty( ar.model(m).data(id).yExp )
            % Observable exists, check if our time point is present.
            selected_tExp = find( ar.model(m).data(id).tExp == timePt );
            selected_tSim = find( ar.model(m).data(id).tFine == timePt );
            
            if ~isempty( selected_tExp ) && ~isempty( selected_tSim )
                if ( debugInfo )
                    ar.model(m).data(id).name
                    printCondition(ar.model(m).data(id).condition);
                end
                
                % Find response variable
                for jc = 1:length(ar.model(m).data(id).condition)
                    if(strcmp(ar.model(m).data(id).condition(jc).parameter, responseVar))
                        jcondi = jc;
                        jval = str2num( ar.model(m).data(id).condition(jcondi).value );
                    end
                end
                
                ovars   = symvar(ar.model(m).data(id).fy{obsidx});

                qLog    = ar.model(m).data(id).logfitting(obsidx);
                scale   = contains(ovars, 'scale_');
                offset  = contains(ovars, 'offset_');

                scale   = arFindPar( ovars{scale}, 'exact' );
                offset  = arFindPar( ovars{offset}, 'exact' );

                if ( numel( scale ) > 1 )
                    error( 'Ambiguous scaling' );
                end
                if ( numel( offset ) > 1 )
                    error( 'Ambiguous offset' );
                end

                if ( numel(scale) == 1 )
                    scale = arGetPars( ar.pLabel{scale}, qLog );
                else
                    warning('NO SCALE FOUND');
                    scale = 1;
                end
                if ( numel(offset) == 1 )
                    offset = arGetPars( ar.pLabel{offset}, 0 );
                    if ( qLog > 0 && ( offset ~= 0 ) )
                        warning( 'This model has an offset in log-space, which means the error estimates will be biased' );
                    end
                else
                    offset = 0;
                end

                % Exp points
                yExp        = ar.model(m).data(id).yExp(selected_tExp, obsidx);
                t           = repmat( jval, size(yExp) );
                yExpStd     = ar.model(m).data(id).ystdExpSimu(selected_tExp, obsidx);
                scaleExp    = repmat( scale, size(yExp) );
                offsetExp   = repmat( offset, size( yExp ) );
                qLogExp     = repmat( qLog, size( yExp ) );

                ts          = [ts; t(:)];
                yExps       = [yExps; yExp(:)];
                yExpStds    = [yExpStds; yExpStd(:)];
                scales      = [scales; scaleExp(:)];
                offsets     = [offsets; offsetExp(:)];
                qLogExps    = [qLogExps; qLogExp(:)] == 1;

                % Simulations
                ySim        = ar.model.data(id).yFineSimu(selected_tSim, obsidx);
                ySimStd     = ar.model.data(id).ystdFineSimu(selected_tSim, obsidx);
                tFine       = repmat( jval, size(ySim) );
                scaleSim    = repmat( scale, size( ySim ) );
                offsetSim   = repmat( offset, size( ySim ) );
                qLogSim     = repmat( qLog, size( ySim ) );

                tfs         = [tfs; tFine(:)];
                ySims       = [ySims; ySim(:)];
                ySimStds    = [ySimStds; ySimStd(:)];
                scalesSim   = [scalesSim; scaleSim(:)];
                offsetsSim  = [offsetsSim; offsetSim(:)];
                qLogSims    = [qLogSims; qLogSim(:)] == 1;
            end
        end
    end
end

% Grab data corresponding to a specific observable in a specific condition
function [ ts, tfs, yExps, ySims, yExpStds, ySimStds, scales, offsets, scalesSim, offsetsSim, qLogExps, qLogSims ] = grabData( m, observable, condis )
    global ar;
    
    debugInfo = 0;
    ds = setdiff( setdiff( arFindData('_', condis{:} ), arFindData('decay') ), arFindData('steady' ) );
    
    ts          = [];
    tfs         = [];
    ySims       = [];
    yExps       = [];
    scales      = [];
    offsets     = [];
    scalesSim   = [];
    offsetsSim  = [];
    yExpStds    = [];
    ySimStds    = [];
    qLogExps    = [];
    qLogSims    = [];
    
    debugInfo = 0;
    
    for jd = 1 : numel( ds )
        id = ds(jd);
        obsidx  = find( strcmp( ar.model(m).data(id).y, observable ) );
        if ( ~ar.model(m).data(id).qFit(obsidx) )
            obsidx = [];
        end
        if ~isempty( obsidx ) && ~isempty( ar.model(m).data(id).yExp )
            
            if ( debugInfo )
                ar.model(m).data(id).name
                printCondition(ar.model(m).data(id).condition);
            end
            ovars   = symvar(ar.model(m).data(id).fy{obsidx});

            qLog    = ar.model(m).data(id).logfitting(obsidx);
            
            scale   = contains(ovars, 'scale_');
            offset  = contains(ovars, 'offset_');
            
            scale   = arFindPar( ovars{scale}, 'exact' );
            offset  = arFindPar( ovars{offset}, 'exact' );
            
            if ( numel( scale ) > 1 )
                error( 'Ambiguous scaling' );
            end
            if ( numel( offset ) > 1 )
                error( 'Ambiguous offset' );
            end
            
            if ( numel(scale) == 1 )
                scale = arGetPars( ar.pLabel{scale}, qLog );
            else
                scale = 1;
            end
            if ( numel(offset) == 1 )
                offsetidx = offset;
                offset = arGetPars( ar.pLabel{offset}, 0 );
            else
                offset = 0;
            end

            s = '';
            for jco = 1 : numel( ar.model(m).data(id).condition )
                fprintf( '%s %s=%s ', s, ar.model(m).data(id).condition(jco).parameter, ar.model(m).data(id).condition(jco).value );
            end
            fprintf( '\n' );
            
            % Exp points
            t           = ar.model(m).data(id).tExp;
            yExp        = ar.model(m).data(id).yExp(:, obsidx);
            yExpStd     = ar.model(m).data(id).yExpStd(:, obsidx);
            
            if isnan(nanmean(yExpStd))
                yExpStd     = ar.model(m).data(id).ystdExpSimu(:, obsidx);
            end
            
            scaleExp    = repmat( scale, size(yExp) );
            offsetExp   = repmat( offset, size( yExp ) );
            qLogExp     = repmat( qLog, size( yExp ) );
           
            ts          = [ts; t(:)];
            yExps       = [yExps; yExp(:)];
            yExpStds    = [yExpStds; yExpStd(:)];
            scales      = [scales; scaleExp(:)];
            offsets     = [offsets; offsetExp(:)];
            qLogExps    = [qLogExps; qLogExp(:)] == 1;
            
            % Simulations
            tFine       = ar.model(m).data(id).tFine;
            ySim        = ar.model.data(id).yFineSimu(:, obsidx);
            ySimStd     = ar.model.data(id).ystdFineSimu(:, obsidx);
            scaleSim    = repmat( scale, size( ySim ) );
            offsetSim   = repmat( offset, size( ySim ) );
            qLogSim     = repmat( qLog, size( ySim ) );
            
            tfs         = [tfs; tFine(:)];
            ySims       = [ySims; ySim(:)];
            ySimStds    = [ySimStds; ySimStd(:)];
            scalesSim   = [scalesSim; scaleSim(:)];
            offsetsSim  = [offsetsSim; offsetSim(:)];
            qLogSims    = [qLogSims; qLogSim(:)] == 1;
        end
    end
end

