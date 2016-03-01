function [ savePath, nRows, nCols ] = arSaveFigure(h, name, subdir)

global ar

useNewExport = 0;
if ( isfield( ar, 'config' ) && isfield( ar.config, 'useNewExport' ) )
    useNewExport = ar.config.useNewExport;
end

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

%% tex generation
set(h,'Units','in')
if ( useNewExport )
    myaxes = unique( [ findobj(h, 'Type', 'Axes'); findobj(h,'Type','Legend') ] );
    
    % Remove axes that are not visible
    rmlist = [];
    for a = 1 : length( myaxes )
        if ( strcmp( get(myaxes(a), 'Visible'), 'off' ) )
            rmlist = [rmlist a]; %ok<AGROW>
        end
    end
    myaxes(rmlist) = [];
else
    myaxes = findobj(h, 'Type', 'axes');
end

axesList = 1 : length( myaxes );

% Manual legend positioning
if ( useNewExport )
    legendList = [];

    for a = 1 : length( myaxes )
        if( strcmp( get(myaxes(a), 'Tag'), 'legend' ) )
            legendList(end+1) = a; %#ok<AGROW>
        end
    end
    axesList(legendList) = [];
end

mypos = get(myaxes,'Position');
if(iscell(mypos))
    mypos = cell2mat(mypos);
end
mycols = unique(round(mypos(axesList,1)*1e6)/1e6);
myrows = unique(round(mypos(axesList,2)*1e6)/1e6);
if ( ~useNewExport )
    if(myrows(end) - max(mypos(:,3)) > 0)
        myrows = [0 myrows'];
    end
end
nCols = length(mycols);
nRows = length(myrows);

if ( useNewExport )
    % No room was allocated specifically for the legend
    if ( ( nRows * nCols ) == ( numel(myaxes) - length( legendList ) ) )
        if ( length(legendList) > 0 )
            nCols = nCols + 1;
        end
    end
end

myfigpos = get(h,'Position');
set(h,'Position',[myfigpos(1:2) 4*nCols 4*nRows]);

maxCol = zeros(nRows, 1);
myLegend = [];
for ia = 1:length(myaxes)
    [~, indCol] = min(abs(mypos(ia,1) - mycols));
    [~, indRow] = min(abs(mypos(ia,2) - myrows));
    
    % Only place the actual axes
    if ( ~useNewExport || ~strcmp( get(myaxes(ia), 'Tag'), 'legend' ) )
        maxCol(indRow) = max( [ maxCol(indRow), indCol ] );
        if ( ~useNewExport )
            set(myaxes(ia),'Position',[(indCol-1)/nCols+(1/nCols)*0.15 (indRow-1)/nRows+(1/nRows)*0.1 (1/nCols)*0.7 (1/nRows)*0.8])
        else
            set(myaxes(ia),'Position',[(indCol-1)/nCols+(1/nCols)*0.25 (indRow-1)/nRows+(1/nRows)*0.2 (1/nCols)*0.7 (1/nRows)*0.7])
        end
    else
        myLegend = myaxes(ia);
    end
end

if ( useNewExport )
    if ~isempty( myLegend )
        [indCol, indRow] = min(maxCol);
        indCol = indCol + 1;

        oldPos  = get(myLegend, 'Position');
        AR = oldPos(3) / oldPos(4);
        sc = .015*numel(get(myLegend, 'String'));
        
        set(myLegend,'Position',[(indCol-1)/nCols+(1/nCols)*0.25 (indRow-1)/nRows+(1/nRows)*0.85-sc sc*AR sc]);
    end
end

axoptions={'x tick label style={/pgf/number format/fixed}','y tick label style={/pgf/number format/fixed}'};
matlab2tikz([savePath '_Report.tex'],'figurehandle',h,'showInfo', false, 'showWarnings',false, 'width','0.9\textwidth','automaticLabels',true,'extraAxisOptions',axoptions)
matlab2tikz([savePath '.tex'],'figurehandle',h,'showInfo', false, 'showWarnings',false, 'standalone', true,'automaticLabels',true,'extraAxisOptions',axoptions)

%% pdf generation
if(ispc || ar.config.useNewPlots==0)
    print('-dpdf', savePath);
    print('-dpng', savePath);
    if ispc
        print('-dmeta', savePath);
    end
elseif(isunix)
  
    % empty LD_LIBRARY_PATH (MATLAB-shipped libraries conflict with libs pdflatex is linked to)
    library_path = getenv('LD_LIBRARY_PATH');
    setenv('LD_LIBRARY_PATH', '');
    
    workingdir = cd;
    cd([arSave subdir]);
    
    if(ismac)
        eval(['!/usr/texbin/pdflatex -interaction=nonstopmode ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
    else
        eval(['!pdflatex -interaction=nonstopmode ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
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
    try
        plot2svg([savePath '.svg'],h);
    end
end


