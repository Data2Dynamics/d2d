function savePath = arSaveFigure(h, name, subdir)

savePath = [arSave subdir];

if(~exist(savePath, 'dir'))
	mkdir(savePath)
end

if(length(name)>60)
    name = name([1:29 (end-29):end]);
end

savePath = arPathConvert([savePath '/' name]);

set(h,'Renderer','painters')

saveas(h, savePath, 'fig');

%% pdf generation
if(ispc)
    print('-dpdf', savePath);
    print('-dpng', savePath);
    print('-dmeta', savePath);
elseif(isunix)
    set(h,'Units','in')
    myaxes = findobj(h,'Type','axes');
    mypos = get(myaxes,'Position');
    if(iscell(mypos))
        mypos = cell2mat(mypos);
    end
    mycols = unique(round(mypos(:,1)*1e6)/1e6);
    myrows = unique(round(mypos(:,2)*1e6)/1e6);
    if(myrows(end) - max(mypos(:,3)) > 0)
        myrows = [0 myrows'];
    end
    nCols = length(mycols);
    nRows = length(myrows);
    
    myfigpos = get(h,'Position');
    set(h,'Position',[myfigpos(1:2) 4*nCols 4*nRows]);

    for ia = 1:length(myaxes)
        [tmp, indCol] = min(abs(mypos(ia,1) - mycols));
        [tmp, indRow] = min(abs(mypos(ia,2) - myrows));

        set(myaxes(ia),'Position',[(indCol-1)/nCols+(1/nCols)*0.15 (indRow-1)/nRows+(1/nRows)*0.1 (1/nCols)*0.7 (1/nRows)*0.8])
    end
    axoptions={'x tick label style={/pgf/number format/fixed}','y tick label style={/pgf/number format/fixed}'};    
    matlab2tikz([savePath '_Report.tex'],'figurehandle',h,'showInfo', false, 'showWarnings',false, 'width','0.9\textwidth','automaticLabels',true,'extraAxisOptions',axoptions)
    matlab2tikz([savePath '.tex'],'figurehandle',h,'showInfo', false, 'showWarnings',false, 'standalone', true,'automaticLabels',true,'extraAxisOptions',axoptions)
    
    % empty LD_LIBRARY_PATH (MATLAB-shipped libraries conflict with libs pdflatex is linked to)
    library_path = getenv('LD_LIBRARY_PATH');
    setenv('LD_LIBRARY_PATH', '');
    
    workingdir = cd;
    cd([arSave subdir]);
    
    if(ismac)
        eval(['!/usr/texbin/pdflatex ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
    else
        eval(['!pdflatex ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
    end
    
    cd(workingdir);
    setenv('LD_LIBRARY_PATH', library_path);
    
    for ia = 1:length(myaxes)
        set(myaxes(ia),'Position',mypos(ia,:))
    end
    set(h,'Position',myfigpos);
end

if(which('plot2svg'))
%     [filepath,filename] = strsplit(savePath,'/');
    plot2svg([savePath '.svg'],h)
end


