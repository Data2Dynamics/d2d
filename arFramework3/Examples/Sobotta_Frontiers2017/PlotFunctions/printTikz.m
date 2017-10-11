function printTikz( filename )

    delete('tikz.pdf');
    qLabels( 'Socs3 mRNA', '', '', [0 20 40 60 80 100 120], [-2 -1.5 -1 -0.5 0], 3, 3 );
    axoptions={'x tick label style={/pgf/number format/fixed}','y tick label style={/pgf/number format/fixed}'};
    matlab2tikz(['tikz/tikzdoc', '.tex'],'figurehandle',gcf,'showInfo', false, 'showWarnings',false, 'width','1.5\textwidth','automaticLabels',true,'extraAxisOptions',axoptions);
    eval(['!pdflatex -interaction=nonstopmode ' [arPathConvert('tikz/tikz') '.tex'] ' > log_pdflatex.txt']);
    copyfile( 'tikz.pdf', filename );
    
end