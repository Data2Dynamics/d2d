function figExport( filename, S )

    if ~exist( 'filename', 'var' )
        filename = 'test';
    end
    
    if (~exist('S', 'var'))
        S = 20;
    end

    % Here we preserve the size of the image when we save it.
    AR = 1;
    set(gcf, 'PaperUnits', 'centimeters' );
    set(gcf, 'PaperSize', [S, AR*S] )
    set(gcf, 'Position', [0 0 S AR*S] );
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperType', '<custom>' );
    set(gcf, 'PaperPosition', [0 0 S AR*S] );
    
    print(filename,'-depsc2','-r300');
    
    if ispc % Use Windows ghostscript call
        system( sprintf( 'gswin64c -o -q -sDEVICE=png256 -dEPSCrop -r300 -o%s.png %s.eps', filename ) );
    else % Use Unix/OSX ghostscript call
        system( sprintf( 'gs -dSAFER -dBATCH -dNOPAUSE -dNOPROMPT -o -q -sDEVICE=pngalpha -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -dUseCIEColor -dEPSCrop -r300 -o%s.png %s.eps', filename, filename ) );
    end
