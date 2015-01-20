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
elseif(isunix)
    set(h,'Units','in')
    
    matlab2tikz([savePath '_Report.tex'],'figurehandle',h,'showInfo', false, 'width','0.9\textwidth')
    matlab2tikz([savePath '.tex'],'figurehandle',h,'showInfo', false, 'standalone', true)
    
    library_path = getenv('LD_LIBRARY_PATH');
    cd([arSave subdir]);
    
    if(ismac)
        eval(['!/usr/texbin/pdflatex ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
    else
        eval(['!pdflatex ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
    end
    
    cd('../../../..');
    setenv('LD_LIBRARY_PATH', library_path);
end