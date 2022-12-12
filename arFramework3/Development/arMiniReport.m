% Create customized report
%
% PLEASE DO NOT EDIT THIS FUNCTION YET. I am still actively working on it
% and I know it's still full of bugs (Joep)
%
% arMiniReport( reportName, option flags )
%
% Option flags are:
%   PlotAll           - Plot all Ys
%   PlotFitted        - Plot only those observables that are fitted.
%   NoExperiments     - Don't export any of the experiments.
%   OmitNonFitted     - Omit plotting the graph and displaying data for
%                       datasets that were not fitted (still shows condition
%                       data however. Useful for steady states)
%   OmitLikelihood    - Omit explicit values for -2 Log(L)
%   OmitNonPlotted    - Omit data/observables that were not plotted
%   KeepRandoms       - Currently, condition replacements based on RANDOMs
%                       are not explicitly mentioned in the tables. If you
%                       do want to see these, specify this keyword
%   KeepFilenames     - Do not shorten data parameters containing entire
%                       filenames (ones that use the _filename tag)
%   AlternateFont     - Use a different font for the report
%   TexPlots          - Use the original tex files for plotting rather than
%                       the pdfs
%   ExcludeDynPars    - Provide cell array of masks to exclude from
%                       parameter shortlist. Cell array must be provided
%                       as subsequent argument
%   TexPages          - Provide files with custom tex pages (provide cell
%                       array of files as next argument)
%   Bibs              - Provide files with additional bibliography items
%                       (provide cell array of files as next argument)
%   Figures           - Copy figure directories into tex output dir.
%                       These can be used in submitted tex files (provide 
%                       a directory as next argument)
%   OverridePlots     - Add custom experiment plots from own directory
%                       (provide a directory as next argument)

function arMiniReport(varargin)

global ar

warning( 'This report functionality is currently in beta status.' );

switches    = { 'PlotAll', 'PlotFitted', 'OmitNonFitted', 'OmitNonPlotted', 'OmitLikelihood', 'KeepRandoms', 'KeepFilenames', 'AlternateFont', 'TexPlots', 'ExcludeDynPars', 'TexPages', 'Figures', 'OverridePlots', 'Bibs', 'FluxNames', 'NoExperiments' };
extraArgs   = [         0,            0,               0,                0,                0,             0,               0,               0,         0,                 1,          1,         1,               1,      1,           0,               0 ];
descriptions = {    { 'Plotting all Ys', '' }, ...
                    { 'Plotting all fitted Ys', '' }, ...
                    { 'Omitting non fitted Ys from report', '' }, ...
                    { 'Omitting non plotted Ys from report', '' }, ...
                    { 'Omitting likelihood values from report', 'Including likelihood values in report' }, ...
                    { 'Displaying RANDOM transforms explicitly', 'Not displaying RANDOM transforms explicitly' }, ...
                    { 'Not shortening data parameters with filenames.', 'Shortening data parameters with filenames (ones that use the _filename tag)' }, ...
                    { 'Using alternate font', '' }, ...
                    { 'Using tex files for incorporating plots', 'Reverting to pre-generated PDFs for plots' }, ...
                    { 'Excluding specific parameters', '' }, ...
                    { 'Including tex pages', '' }, ...
                    { 'Including figure directories', '' }, ...
                    { 'Including custom plots', '' }, ...
                    { 'Including custom bibliography items', '' }, ...
                    { 'Including flux names instead of indices', '' }, ...
                    { 'Omitting experiments', '' } , ...
                    };

if( (nargin > 0) && max( strcmpi( varargin{1}, switches ) ) == 0 )
    project_name = varargin{1};
    varargin = varargin(2:end);
else
    project_name = 'Data 2 Dynamics Software -- Modeling Report';
end

opts = argSwitch( switches, extraArgs, descriptions, 1, varargin );
if(isempty(ar))
    error('please initialize by arInit')
end

opts.funcwrap = 100;

if ( opts.plotall && opts.plotfitted )
    error( 'Plotall and plotfitted are mutually exclusive options' );
end

if ( opts.plotall )
    disp( 'Plotting all Ys' );
    for jm = 1 : length( ar.model )
        ar.model(jm).qPlotYs = ar.model(jm).qPlotYs > -1;
        arPlotY(true, false, true, 'hidell', 'nameonly');
        close all;
    end
end

if ( opts.plotfitted )
    disp( 'Plotting fitted Ys' );
    for jm = 1 : length( ar.model )
        if ( isfield( ar.model(jm), 'plot' ) )
            for jp = 1 : length( ar.model(jm).plot )
                ar.model(jm).qPlotYs(jp) = isFitted(jm,jp);
            end
            arPlotY(true, false, true, 'hidell', 'nameonly');
            close all;
        end
    end
end

pti = 1 - opts.keepfilenames;

% Fetch MATLAB version
matVer = arVer;
ar.config.matlabVersion = str2double(matVer.Version);

savePath = [arSave '/Latex'];
if ~exist(savePath, 'dir')
    mkdir( savePath );
end

if opts.figures
    for a = 1 : length( opts.figures_args )
        copyfile( [opts.figures_args{a} '/*.pdf'], [savePath '/' opts.figures_args{a} '/'] );
    end
end

if(~exist([cd '/' savePath], 'dir'))
    mkdir([cd '/' savePath])
end

% empty LD_LIBRARY_PATH (MATLAB-shipped libraries conflict with libs pdflatex is linked to)
library_path = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', '');

% copy lib.bib
copyfile(which('lib.bib'), [savePath '/lib.bib']);

if (isfield(ar, 'additionalBib'))
    try
        txt = fileread( ar.additionalBib );
        
        fid = fopen( [savePath '/lib.bib'], 'a' );
        fprintf(fid, txt);
        fclose(fid);
    catch
        warning( 'Failed to read additional bib file' );
    end
end

if opts.bibs
    fid = fopen( [savePath '/lib.bib'], 'a' );
	for a = 1 : length( opts.bibs_args )
        txt = strrep(strrep(fileread(opts.bibs_args{a}), '\', '\\'), '%', '%%');
        fprintf( fid, txt );
	end
    fclose(fid);
end

% latex packages
copyfile(which('assurechemist.sty'), [savePath '/assurechemist.sty']);
copyfile(which('chemist.sty'), [savePath '/chemist.sty']);

% latex file
fname = 'MiniReport.tex';
fnamebib = 'MiniReport.aux';
fprintf('writing latex, file %s...', fname);
fid = fopen([savePath '/' fname], 'w');

%% Head
lp(fid, '\\nonstopmode');
lp(fid, '\\documentclass[10pt, oneside, fleqn]{article}');

lp(fid, '\\usepackage[utf8]{inputenc}');
lp(fid, '\\usepackage{amsmath, amsthm, amssymb, amstext}');
lp(fid, '\\usepackage{listings} ');
lp(fid, '\\usepackage{epsfig} ');
lp(fid, '\\usepackage{graphicx} ');
lp(fid, '\\usepackage{rotating} ');
lp(fid, '\\usepackage{lscape}');
lp(fid, '\\usepackage{color}');
lp(fid, '\\usepackage{natbib}');
lp(fid, '\\usepackage{lmodern} ');
lp(fid, '\\usepackage{xcolor} ');
lp(fid, '\\usepackage{chemist} ');
lp(fid, '\\usepackage{calc} ');
lp(fid, '\\usepackage{tabularx} ');
lp(fid, '\\usepackage{placeins} ');
lp(fid, '\\usepackage{pgfplots}  ');
lp(fid, '\\usepackage{varwidth} ');
lp(fid, '\\usepackage{float} ');
lp(fid, '\\usepackage[normal,footnotesize,bf,singlelinecheck=false]{caption}');
lp(fid, '\\usepackage{subfig}');
lp(fid, '\\usepackage{sidecap}');
lp(fid, '\\usepackage{flafter}' );
lp(fid, '\\usepackage{multirow}' );
lp(fid, '\\usepackage{colortbl}' );
lp(fid, '\\newcommand{\\altrowcol}{\\rowcolor[gray]{0.925}}');
lp(fid, '\\newcommand{\\titlerowcol}{\\rowcolor[gray]{0.825}}');
%lp(fid, '\\captionsetup[subfloat]{position=top}');
lp(fid, '\\allowdisplaybreaks');
lp(fid, '\\setlength{\\textheight}{22 cm}');
lp(fid, '\\setlength{\\textwidth}{16 cm}');
lp(fid, '\\setlength{\\topmargin}{-1.5 cm}');
lp(fid, '\\setlength{\\hoffset}{-2 cm}');
%lp(fid, '\\renewcommand*\\rmdefault{cmss} ');
lp(fid, '\\usepackage[colorlinks=true, linkcolor=blue, citecolor=blue, filecolor=blue, urlcolor=blue]{hyperref}\n');

lp(fid, '\n\\definecolor{mygray}{gray}{.5}');
lp(fid, '\n\\definecolor{red}{rgb}{.8,0,0}');
farben = lines(7);
for j=1:length(farben(:,1))
    lp(fid, '\\definecolor{line%i}{rgb}{%f,%f,%f}', j, farben(j,1), farben(j,2), farben(j,3));
end

lp(fid, '\\newcommand{\\toprule}{}'); %\\hline\\hline
lp(fid, '\\newcommand{\\midrule}{}'); %\\hline
lp(fid, '\\newcommand{\\botrule}{}'); %\\hline\\hline
lp(fid, '\\newcommand{\\dobegincenter}{\\begin{center}}');
lp(fid, '\\newcommand{\\doendcenter}{\\end{center}}');
lp(fid, '\\newcommand{\\mycaption}[3]{\\caption{\\textbf{#1}\\\\ #3 \\label{#2}}}');
lp(fid, '\\newcommand{\\mycaptionof}[3]{\\captionof{table}{\\textbf{#1}\\\\ #3 \\label{#2}}}');
lp(fid, '\\newenvironment{statictable}{}{}\n');

%lp(fid, '\\newenvironment{nobreak}');
%lp(fid, '\t{\\par\\nobreak\\vfil\\penalty0\\vfilneg')
%lp(fid, '\t\\vtop\\bgroup}')
%lp(fid, '\t{\\par\\xdef\\tpd{\\the\\prevdepth}\\egroup')
%lp(fid, '\t\\prevdepth=\\tpd}')  

lp(fid, '\\pgfplotsset{compat=newest} ');
lp(fid, '\\pgfplotsset{plot coordinates/math parser=false} ');
lp(fid, '\\newlength\\figureheight ');
lp(fid, '\\newlength\\figurewidth ');

if ( opts.alternatefont )
    lp(fid, '\\renewcommand*\\rmdefault{cmss} ');
end
lp(fid, '\\newcommand{\\crule}[1]{\\multispan{#1}{\\hspace*{\\tabcolsep}\\hrulefill\\hspace*{\\tabcolsep}}}');

lp(fid, '\n\\begin{document}\n');


lp(fid, '\\newcolumntype{K}{>{\\centering\\arraybackslash}X}');

%% Header
% lp(fid, ['\\title{Data-2-Dynamics Software' ...
%     '\\footnote{Website: \\href{http://data2dynamics.org}' ...
%     '{\\url{http://data2dynamics.org}} \\\\ ', ...
%     'Reference: \\citet{Raue:2012zt}}'...
%     ' \\\\ Modeling Report}']);

lp(fid, '\\title{%s}', project_name);
lp(fid, '\\maketitle\n');
lp(fid, '\\tableofcontents');

N = 10;

for jm=1:length(ar.model)
    lp(fid, '\\clearpage\n');
    lp(fid, '\\section{Model: %s}\n', arNameTrafo(ar.model(jm).name));
    
    targetConditions = [];
    preEq = ones( length(ar.model(jm).x), 1 );
    usesPreEquilibration = 0;
    if ( isfield(ar.model(jm), 'ss_condition') )
        usesPreEquilibration = 1;
        eqType = 'some';
        
        for ss = 1 : length( ar.model(jm).ss_condition )
            targetConditions = union( targetConditions, ar.model(jm).ss_condition(ss).ssLink );
            preEq = ar.model(jm).ss_condition(ss).ssStates.' & preEq;
        end
        if ( length( targetConditions ) == length( ar.model(jm).condition ) )
            % All conditions are pre-equilibrated
            usesPreEquilibration = 2;
            eqType = 'all';
        end
    end
    if ( usesPreEquilibration > 0 )
        if ( sum( preEq ) == length( ar.model(jm).x ) )
            eqStateCaption = sprintf( 'Note that all dynamic variables are simulated to steady state prior to perturbation in %s of the simulation experiments.', eqType );
        else
            eqStateCaption = sprintf( 'Note that initial conditions denoted by an asterisk are simulated to steady state prior to perturbation in %s of the simulation experiments.', eqType );
        end
        eqText = sprintf( 'Prior to stimulation, the cells are assumed to be in a resting (i.e. unstimulated) state. To achieve such an initial steady state, dynamic variables in the model (see Table \\ref{variables}) are pre-equilibrated prior to each simulation experiment. Model pre-equilibration was performed by simulating until the right hand side of the differential equations dropped below the equilibration threshold (%.0s).', ar.config.eq_tol );
    else
        eqText = '';
        eqStateCaption = '';
    end
    
    %% descriptions
    if(~isempty(ar.model(jm).description))
        for jd=1:length(ar.model(jm).description)
            lp(fid, '%s\\\\', strrep(strrep(ar.model(jm).description{jd}, '%', '\%'), '_', '\_'));
        end
    end
    
    %% Introduction
    lp(fid, '\\subsection{Model description}');
    lp(fid, 'The model used in this study is based on a system of Ordinary Differential Equations (ODE). These ordinary differential equations are derived by means of the law of mass-action. ' );
    lp(fid, 'The time evolution of the biochemical compounds is computed by numerically integrating these differential equations. The model contains parameters which are estimated by calibrating the model to data using a Maximum-Likelihood estimation approach. %s', eqText );
    lp(fid, 'All analyses were performed using the \\emph{Data 2 Dynamics} software package \\cite{raue2015data2dynamics}, which is available from \\href{http://data2dynamics.org}{\\url{http://data2dynamics.org}}. ');
    if ( length( ar.model(jm).c ) > 2 )
        lp( fid, 'The model considers multiple compartments (see Table \\ref{compartments}), which are taken into account by considering the relevant compartment volumes in the flux expressions. ' );
    else
        if ( length( ar.model(jm).c ) == 2 )
            lp( fid, 'The model considers two compartments namely %s (vol = %s) and %s (vol = %s), which are taken into account by considering the relevant compartment volumes in the flux expressions. ', ar.model(jm).c{1}, getVolume(jm, 1), ar.model(jm).c{2}, getVolume(jm, 2) );
        end
    end
    lp(fid, 'The %d dynamic variables used in the model are summarized in Table \\ref{variables}', length(ar.model(jm).x) );
    
    if ( length( ar.model(jm).u ) > 0 ) %#ok
        lp(fid, ', while the %i external inputs variables are summarized in Table \\ref{inputs}.', sum(~strcmp(ar.model(jm).fu, '0')));
    else
        lp(fid, '.', length(ar.model(jm).u));
    end
    
    %% compartments
    if ( length( ar.model(jm).c ) > 2 )       
        lp(fid, '\\begin{table}[h]');
        lp(fid, '\\centering');
        lp(fid, '\\begin{tabular}{@{} *2l @{}}\\toprule');
        lp(fid, '\\titlerowcol \\textbf{Compartment} & \\textbf{Volume}\\tabularnewline\\midrule' );
        for jc = 1 : length( ar.model(jm).c )
            volP = getVolume(jm, jc);
            lp(fid, '%s%s & %s\\\\', alternate(jc), ar.model(jm).c{jc}, volP );
        end
        lp(fid, '\\botrule\\end{tabular}');
        lp(fid, '\\caption[align=left,justification=justified,singlelinecheck=false]{Model compartments}\\label{compartments}');
        lp(fid, '\\end{table}');
    else
        if ( length( ar.model(jm).c ) == 1 )
            lp( fid, 'Concentrations were based on a compartment volume of %f.', getVolume(jm, 1) );
        end
    end
    
    %% species
    % Did the user fill in better descriptions for the states?
    if ~(min(strcmp(ar.model(jm).xNames, ar.model(jm).x)))
        headerChunk = '& \textbf{Description}';
        nC = 5;
    else
        headerChunk = '';
        nC = 4;
    end
    
    if(~isempty(ar.model(jm).x))
        lp(fid, '', length(ar.model(jm).x) );
        lp(fid, '\\begin{table}[h]');
        lp(fid, '\\begin{tabular}{@{} *%dl @{}}\\toprule', nC);
        lp(fid, '\\titlerowcol \\textbf{Variable} & \\textbf{Unit} & \\textbf{Compartment} & \\textbf{Initial Condition} %s\\tabularnewline\\midrule', headerChunk );
        for jx = 1 : length( ar.model(jm).x )
            if ( ( usesPreEquilibration > 0 ) && preEq(jx) && ( sum( preEq ) ~= length( ar.model(jm).x ) ) )
                asterisk = '*';
            else
                asterisk = '';
            end
            
            x0 = ar.model(jm).px0{jx};
            x0 = char(arSym(ar.model(jm).fp{strcmp(ar.model(jm).p, x0)}));
            if ~strcmp( headerChunk, '' )
                lp(fid, '%s%s & %s [%s] & %s & %s%s & %s\\\\', alternate(jx), strrep(ar.model(jm).x{jx}, '_', '\_'), ar.model(jm).xUnits{jx,3}, ar.model(jm).xUnits{jx,2}, ar.model(jm).c{ar.model(jm).cLink(jx)}, strrep(x0, '_', '\_'), asterisk, ar.model(jm).xNames{jx} );
            else
                if isempty(ar.model(jm).c)
                    lp(fid, '%s%s & %s [%s] & %s & %s%s\\\\', alternate(jx), strrep(ar.model(jm).x{jx}, '_', '\_'), ar.model(jm).xUnits{jx,3}, ar.model(jm).xUnits{jx,2}, '', strrep(x0, '_', '\_'), asterisk );
                else
                    lp(fid, '%s%s & %s [%s] & %s & %s%s\\\\', alternate(jx), strrep(ar.model(jm).x{jx}, '_', '\_'), ar.model(jm).xUnits{jx,3}, ar.model(jm).xUnits{jx,2}, ar.model(jm).c{ar.model(jm).cLink(jx)}, strrep(x0, '_', '\_'), asterisk );
                end
            end
        end
        lp(fid, '\\end{tabular}');
        lp(fid, '\\caption[align=left,justification=justified,singlelinecheck=false]{Dynamic variables used in the model. %s}\\label{variables}', eqStateCaption);
        lp(fid, '\\end{table}');
    end
    
    inputs = 1;
    graph = 0;
    obs = 1;
    
    %% inputs
    if(~isempty(ar.model(jm).u)&&inputs)
        
        % Did the user fill in better descriptions for the states?
        if (max( cellfun(@length,ar.model(jm).uNames) ) > 0)
            headerChunk = '& \textbf{Description}';
            nC = 4;
        else
            headerChunk = '';
            nC = 3;
        end        
        
        % Find the non-empty inputs
        us = find(~strcmp(ar.model(jm).fu,'0'));
        
        if ~isempty(us)
            lp(fid, '\\begin{table}');
            startFlexbox(fid, 'ainputbox');
            lp(fid, '\\begin{tabular}{@{} *%dl @{}}\\toprule', nC);
            lp(fid, '\\titlerowcol \\textbf{Input} & \\textbf{Unit} & \\textbf{Default equation} %s\\tabularnewline\\midrule', headerChunk );
            if ( ~strcmp( headerChunk, '' ) )
                for ku = 1 : length( us )
                    ju = us(ku);
                    inp = repFunc( ar.model(jm).fu{ju}, 'monospline10' );
                    lp(fid, '%s%s & %s [%s] & $%s$ & %s\\\\', alternate(ku), strrep( ar.model(jm).u{ju}, '_', '\_'), ar.model(jm).uUnits{ju,3}, ar.model(jm).uUnits{ju,2}, myFormulas(inp, jm), ar.model(jm).uNames{ju} );
                end
            else
                for ku = 1 : length( us )
                    ju = us(ku);
                    inp = repFunc( ar.model(jm).fu{ju}, 'monospline10' );
                    lp(fid, '%s%s & %s [%s] & $%s$\\\\', alternate(ku), strrep( ar.model(jm).u{ju}, '_', '\_'), ar.model(jm).uUnits{ju,3}, ar.model(jm).uUnits{ju,2}, myFormulas(inp, jm) );
                end
            end
            lp(fid, '\\botrule\\end{tabular}');
            endFlexbox(fid, 'ainputbox');
            lp(fid, '\\caption[align=left,justification=justified,singlelinecheck=false]{Inputs used in the model}\\label{inputs}');
            lp(fid, '\\end{table}');
        end
    end
    
    %% reactions
    %vs = sym('v', [1, size(ar.model(jm).N,2)]);
    %cs = sym(strcat('vol_', ar.model(jm).c));
    lp(fid, '\\newpage\\noindent The model consists of %d differential equations, which are given by the following equations:', length(ar.model(jm).x) );
    lp(fid, '{\\footnotesize');
    lp(fid, '\\begin{align}');
    for jx=1:size(ar.model(jm).N, 1) % for every species jx
        strtmp = '';
        if(~isempty(ar.model(jm).c))
            qinfluxwitheducts = ar.model(jm).N(jx,:) > 0 & sum(ar.model(jm).N < 0,1) > 0;
            eductcompartment = zeros(size(qinfluxwitheducts));
            for jj=find(qinfluxwitheducts)
                eductcompartment(jj) = unique(ar.model(jm).cLink(ar.model(jm).N(:,jj)<0)); %R2013a compatible
            end
        end
        f = 0;
        for jv = find(ar.model(jm).N(jx,:))
            if ( opts.fluxnames )
                fname = ar.model(jm).v{jv};
            else
                fname = sprintf( '%i', jv );
            end
            
            if(abs(ar.model(jm).N(jx,jv))~=1)
                
                strtmp = [strtmp sprintf(' %+i \\cdot v_{%s}', ar.model(jm).N(jx,jv), fname)]; %#ok<*AGROW>
                f = 1;
            elseif(ar.model(jm).N(jx,jv)==1)
                if ( f == 0 )
                    strtmp = [strtmp sprintf(' v_{%s}', fname)];
                else
                    strtmp = [strtmp sprintf(' + v_{%s}', fname)];
                end
                f = 1;
            elseif(ar.model(jm).N(jx,jv)==-1)
                
                strtmp = [strtmp sprintf(' - v_{%s}', fname)];
                f = 1;
            end
            if(~isempty(ar.model(jm).c) && qinfluxwitheducts(jv) && eductcompartment(jv)~=ar.model(jm).cLink(jx))
                strtmp = [strtmp sprintf(' \\cdot \\frac{%s}{%s}', myFormulas(ar.model(jm).pc{eductcompartment(jv)}, jm), ...
                    myFormulas(ar.model(jm).pc{ar.model(jm).cLink(jx)}, jm))];
            end

        end
        if(jx==size(ar.model(jm).N, 1) || mod(jx,N)==0)
            lp(fid, '\t\\mathrm{d}%s/\\mathrm{dt} &= %s \\label{%s}', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_ode%i', ar.model(jm).name, jx));
        else
            lp(fid, '\t\\mathrm{d}%s/\\mathrm{dt} &= %s \\label{%s} \\\\', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_ode%i', ar.model(jm).name, jx));
        end

        if(mod(jx,N)==0 && jx<size(ar.model(jm).N, 1))
            lp(fid, '\\end{align}\n');
            lp(fid, '\\begin{align}');
        end
    end
    lp(fid, '\\end{align}}\n\n');

    %% derived
    if(~isempty(ar.model(jm).z)&&1)
        if ( length( ar.model(jm).z ) > 1 )
            lp(fid, '\\noindent Based on numerical integration of these equations, %i derived quantities are computed:', length(ar.model(jm).z));
        else
            lp(fid, '\\noindent Based on numerical integration of these equations, the following derived quantity is computed:' );
        end
        lp(fid, '{\\footnotesize');
        lp(fid, '\\begin{align}');    
        for jz = 1:length(ar.model(jm).z)            
            lp(fid, '%s(%s) &= %s \\label{%s}\\\\', ...
                myFormulas(ar.model(jm).z{jz}, jm), ...
                myFormulas(ar.model(jm).t, jm), ...
                myFormulas(ar.model(jm).fz{jz}, jm), ...
                sprintf('%s_derived%i', ar.model(jm).name, jz));
        end
        lp(fid, '\\end{align}}');
    end    
    
    %% flux expressions
    lp(fid, '\\noindent The flux expressions corresponding to these equations are provided in Table \\ref{fluxes}. The system of ODEs was integrated using the CVODES algorithm from the SUNDIALS suite of solvers \\cite{Hindmarsh2005fb}. ');
    if ( ar.config.useSensis );
        lp(fid, 'First order derivatives were computed using the sensitivity equations and used for numerical optimization. ');
    end
    lp(fid, 'Relative and absolute tolerances were set to %.0e and %.0e respectively.\\\\\n\n', ar.config.rtol, ar.config.atol);
       
    lp(fid, '\\begin{table}[h]');
    startFlexbox(fid, 'fluxbox');
    lp(fid, '\\begin{tabular}{@{} *5l @{}}\\toprule');
    lp(fid, '\\titlerowcol \\textbf{Flux} & \\textbf{Equation} %s\\tabularnewline\\midrule', headerChunk );
    for jv = 1 : length( ar.model(jm).fv )
        if ( opts.fluxnames )
            fname = ar.model(jm).v{jv};
        else
            fname = sprintf( '%i', jv );
        end
        
        if ~strcmp( headerChunk, '' )
            lp(fid, '%s$v_{%s}$ & $%s$ & %s\\\\', alternate(jv), fname, myFormulas(ar.model(jm).fv{jv}, jm), ar.model(jm).v{jv} );
        else
            lp(fid, '%s$v_{%s}$ & $%s$\\\\', alternate(jv), fname, myFormulas(ar.model(jm).fv{jv}, jm) );
        end
    end
    lp(fid, '\\end{tabular}');
    endFlexbox(fid, 'fluxbox');
    lp(fid, '\\caption[align=left,justification=justified,singlelinecheck=false]{Model flux expressions}\\label{fluxes}');
    lp(fid, '\\end{table}');    
    
    
    %% Structure
    if(isfield(ar.model(jm), 'savePath_Graph') && ~isempty(ar.model(jm).savePath_Graph) && graph)
        savePath_Graph = [arSave '/Figures/Network' sprintf('/%s.pdf', ar.model(jm).name)];
        if(exist(savePath_Graph,'file'))
            lp(fid, '\\subsection{Model structure}');
            lp(fid, 'The model structure is depicted in Figure \\ref{%s}.', [ar.model(jm).name '_graph']);
            
            copyfile(savePath_Graph, [savePath '/' ar.model(jm).name '_graph.pdf'])
            captiontext = sprintf('\\textbf{%s network representation.} ', arNameTrafo(ar.model(jm).name));
            if(~isempty(ar.model(jm).u))
                captiontext = [captiontext 'Diamond shaped nodes correspond to model inputs, see Equation '];
                if(length(ar.model(jm).u)==1)
                    captiontext = [captiontext '\ref{' sprintf('%s_input%i', ar.model(jm).name, 1) '}. '];
                else
                    captiontext = [captiontext '\ref{' sprintf('%s_input%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_input%i', ar.model(jm).name, length(ar.model(jm).u)) '}. '];
                end
            end
            captiontext = [captiontext 'Ellipsoid shaped nodes correspond to dynamical variables described by the ODE system, see Equation '];
            captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
            captiontext = [captiontext 'Black arrows and box shaped nodes indicate reactions with corresponding rate equations given in Equation '];
            captiontext = [captiontext '\ref{' sprintf('%s_flux%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_flux%i', ar.model(jm).name, length(ar.model(jm).fv)) '}. '];
            captiontext = [captiontext 'Red T-shaped arrows indicated inhibitorial influence on a reaction and blue O-shaped arrows indicate catalysing influence on a reaction. '];
            lpfigure(fid, 0.5, [ar.model(jm).name '_graph.pdf'], captiontext, [ar.model(jm).name '_graph']);
        end
    end
        
   
    %% standard observations and error model
    if(isfield(ar.model(jm), 'y')&&obs)
        lp(fid, '\\subsection{Observables}\n');
        
        lp(fid, '\\noindent The model contains %i standard observables listed in table \\ref{defaultObservables}. Certain experiments may contain experiment specific observation functions, which are detailed in the experiment section (see \\ref{ExperimentsSection}).\\\\\n', length(ar.model(jm).y));
        
        lp(fid, '\\begin{statictable}');
        lp(fid, '\t\\centering');
        startFlexbox(fid, 'obstable' );
        lp(fid, '\\begin{tabular}{@{} lll @{}}\\toprule');
        lp(fid, '\\titlerowcol \\textbf{Observable} & & \\textbf{Equations} \\tabularnewline\\midrule' );                        

        for jy = 1 : length( ar.model(jm).y )
            strtmp = myFormulas(ar.model(jm).fy{jy}, jm);
            if(ar.model(jm).logfitting(jy))
                strtmp = ['\mathrm{log}_{10}(' strtmp ')'];
            end            
            if ( strcmp( ar.model(jm).yNames{jy}, '' ) )
                lp(fid, ' %s & y & $%s$ \\tabularnewline %s \\multirow{-2}{*}{%s [%s]} & $\\sigma$ & $%s$ \\tabularnewline\n', alternate(jy), strtmp, alternate(jy), strrep(ar.model(jm).y{jy}, '_', '\_\-'), arNameTrafo( ar.model(jm).yUnits{jy,2} ), myFormulas(ar.model(jm).fystd{jy}, jm) );
            else
                lp(fid, ' %s & y & $%s$ \\tabularnewline %s \\multirow{-2}{*}{%s [%s]} & $\\sigma$ & $%s$ \\tabularnewline\n', alternate(jy), strtmp, alternate(jy), ar.model(jm).yNames{jy}, arNameTrafo( ar.model(jm).yUnits{jy,2} ), myFormulas(ar.model(jm).fystd{jy}, jm) );
            end
        end
        lp(fid, '\\botrule\\end{tabular}');
        endFlexbox(fid, 'obstable' );
        lp(fid, '\\mycaptionof{Model observables and error models}{defaultObservables}{}' );
        lp(fid, '\\end{statictable}');       
   
    end
    
    %% Transforms
    ccount = 1;
    
    lp(fid, '\\subsection{Default parameter transformations}');
    lp(fid, '\\label{defaultpartransforms}');
    lp(fid, 'Some of the bare model parameters are transformed before each simulation experiment. ');
    lp(fid, 'This is done both to facilitate fitting (e.g. decoupling scale from dynamics), but also ')
    lp(fid, 'to allow condition specific parameter changes.');
    lp(fid, '\\noindent The ODE system is modified by the following parameter transformations:');
    
    lp(fid, '{\\footnotesize');
    lp(fid, '\\begin{align}');
    for jp=1:length(ar.model(jm).fp)
        if(~strcmp(ar.model(jm).p{jp}, ar.model(jm).fp{jp}))
            % Which parameters actually make it to the outside?
            skip = 0;
            if ( usesPreEquilibration == 2 )
                % This parameter is an initial condition and equilibrated
                if ( max( strcmp( ar.model(jm).p{jp}, strcat( 'init_', ar.model(jm).x(preEq) ) ) ) )
                    % It's zero, so not relevant.
                    if (str2num(ar.model(jm).fp{jp}) == 0) %#ok
                        skip = 1;
                    end
                end
            end
            
            if ( strcmp( ar.model(jm).p{jp}, strrep(strrep(ar.model(jm).fp{jp},'(',''),')','') ) )
                skip = 1;
            end            
            
            if ( ~skip )
                if(ccount==length(ar.model(jm).fp) || mod(ccount,N)==0)
                    lp(fid, '\t%s & \\rightarrow %s ', myFormulas(PTI(jm, ar.model(jm).p{jp},pti), jm), ...
                        myFormulas(PTI(jm, ar.model(jm).fp{jp},pti), jm));
                else
                    lp(fid, '\t%s & \\rightarrow %s \\\\', myFormulas(PTI(jm, ar.model(jm).p{jp},pti), jm), ...
                        myFormulas(PTI(jm, ar.model(jm).fp{jp},pti), jm));
                end
                if(mod(ccount,N)==0 && ccount<length(ar.model(jm).fp))
                    lp(fid, '\\end{align}\n');
                    lp(fid, '\\begin{align}');
                end
                ccount = ccount + 1;
            end
        end
        if(ccount>1 && jp==length(ar.model(jm).fp))
            lp(fid, '\\end{align}}\n\n');
        end
    end
    
    
    lp(fid, '\\subsection{Dynamic parameters}');
    if(isfield(ar,'ndata') && isfield(ar,'chi2fit') && isfield(ar,'chi2'))
    if(ar.config.fiterrors == 1)
        if ( opts.omitlikelihood )
            llterm = sprintf( ' based on a total of %i data points', ar.ndata );
        else
            llterm = sprintf( ', yielding a value of the objective function $-2 \\log(L) = %g$ for a total of %i data points', 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit, ar.ndata );
        end
        lp(fid, 'In total %i parameters are estimated from the experimental data%s.', sum(ar.qFit==1), llterm);
    else
        if ( opts.omitlikelihood )
            llterm = sprintf( ' based on a total of %i data points', ar.ndata );
        else
            llterm = sprintf( ', yielding a value of the objective function $\\chi^2 = %g$ for a total of %i data points', ar.chi2, ar.ndata );
        end        
        lp(fid, 'In total %i parameters are estimated from the experimental data%s.', sum(ar.qFit==1), llterm);
    end

    lp(fid, 'The model parameters were estimated by maximum likelihood estimation applying the MATLAB lsqnonlin algorithm.');
    end
    
    N = 51;    
    ntables = ceil(length(ar.p)/N);

    lp(fid, 'The model parameters which influence system dynamics are listed in Table \\ref{dynpars}.');
    lp(fid, 'Parameters highlighted in red color indicate parameter values close to their bounds.');
    lp(fid, 'An extensive list of all the estimated parameters (including all observational and error model parameters) is given in section \\ref{estimatedparameters}.' )
    if(ntables>1)
        lp(fid, 'These can be found in Table \\ref{paratable1} -- \\ref{paratable%i}.', ntables);
    else
        lp(fid, 'These can be found in Table \\ref{paratable1}.');
    end

    pTrans = ar.p;
    pTrans(ar.qLog10==1) = 10.^pTrans(ar.qLog10==1);

    dynPars = find(ar.qDynamic==1);
    lp(fid, '\t\\begin{table}[ht!]');
    lp(fid, '\t\\dobegincenter');
    startFlexbox(fid, 'dyntable' );
    lp(fid, '\t{\\footnotesize');
    lp(fid, '\t\t\\begin{tabular}{llllllll}');
    lp(fid, '\t\t\t\\toprule');
    lp(fid, '\\titlerowcol \t\t\t & name & $\\theta_{min}$ & $\\hat\\theta$ & $\\theta_{max}$ & log & non-log $\\hat \\theta$ & fitted \\tabularnewline');
    lp(fid, '\t\t\t\\midrule');
    count = 1;
    for dp=1:length(dynPars)   
        j = dynPars(dp);
        if ( ~masktest( ar.pLabel{j}, opts.excludedynpars_args ) )
            count = count + 1;
            if(ar.qFit(j)==1)
                if(ar.p(j) - ar.lb(j) < 0.1 || ar.ub(j) - ar.p(j) < 0.1)
                    lp(fid, '%s\t\t\t\\color{red}{%i} & \\color{red}{%s} & \\color{red}{%+8.4g} & \\color{red}{%+8.4f} & \\color{red}{%+8.4g} & \\color{red}{%i} & \\color{red}{%s} & \\color{red}{%i} \\tabularnewline', ...
                        alternate(dp), j, strrep(PTI(jm, ar.pLabel{j},pti),'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
                        ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
                else
                    lp(fid, '%s\t\t\t%i & %s & {%+8.4g} & {%+8.4f} & {%+8.4g} & %i & %s & %i \\tabularnewline', ...
                        alternate(dp), j, strrep(PTI(jm, ar.pLabel{j},pti),'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
                        ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
                end
            else
                lp(fid, '%s\t\t\t\\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%+8.4g} & \\color{mygray}{%+8.4f} & \\color{mygray}{%+8.4g} & \\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%i} \\tabularnewline', ...
                    alternate(dp), j, strrep(PTI(jm, ar.pLabel{j},pti),'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
                    ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
            end
        end
    end
    lp(fid, '\t\t\t\\botrule');
    lp(fid, '\t\t\\end{tabular}}');
    endFlexbox(fid, 'dyntable');
    lp(fid, '\t\t\\mycaption{Estimated dynamic parameter values}{dynpars}', count);
    lp(fid, '{$\\hat \\theta$ indicates the estimated value of the parameters.');
    lp(fid, '$\\theta_{min}$ and $\\theta_{max}$ indicate the upper and lower bounds for the parameters.');
    lp(fid, 'The log-column indicates if the value of a parameter was log-transformed.');
    lp(fid, 'If log = 1 the non-log-column indicates the non-logarithmic value of the estimate.');
    lp(fid, 'The fitted-column indicates if the parameter value was estimated (1), was temporarily fixed (0) or if its value was fixed to a constant value (2).}');
    lp(fid, '\t\\doendcenter');
    lp(fid, '\t\\end{table}');
    
    %% Additional description
    if isfield( ar, 'additionalTex' )
        lp(fid, '\\FloatBarrier');
        lp(fid, '\\clearpage\n');
        lp(fid, '\\section{Additional description}');
        lp(fid, strrep(fileread(ar.additionalTex), '\', '\\'));
        lp(fid, '\n');
    end
    
    %% Additional description via args
    if opts.texpages
        lp(fid, '\\FloatBarrier');
        for a = 1 : length( opts.texpages_args )
            lp(fid, strrep(strrep(fileread(opts.texpages_args{a}), '\', '\\'), '%', '%%'));
        end
    end  
    
    %% Experiments
    if ( opts.noexperiments == 0 )
        lp(fid, '\\FloatBarrier');
        lp(fid, '\\clearpage\n');
        lp(fid, '\\section{Experiments}\n');
        lp(fid, '\\label{ExperimentsSection}')
        lp(fid, 'The following pages detail the different simulation experiments that the model was calibrated and/or validated with.');
        lp(fid, 'Each page contains a brief description of the experiment, a plot with the model simulation, the raw data and additional ');
        lp(fid, 'information pertaining to the simulation protocol. Different simulation protocols are encoded by applying condition specific ');
        lp(fid, 'parameter and input transformations, which are listed in each of the experiment descriptions. Unlisted parameters either ');
        lp(fid, 'default to the transformation listed in the main model description \\ref{defaultpartransforms}, or to the raw parameter value ');
        lp(fid, 'as specified in the full table of model parameters (see \\ref{estimatedparameters}).');
        lp(fid, 'To get a global overview of the data involved in model parameterization, please refer to Table \\ref{stateData},');
        lp(fid, 'for a description of which experiments refer to which states.');
        
        if(~isempty(ar.model(jm).x))
            lp(fid, '', length(ar.model(jm).x) );
            lp(fid, '\\begin{table}[h]');
            lp(fid, '\\begin{tabularx}{\\textwidth}{@{} lp{.8\\textwidth} @{}}\\toprule');
            lp(fid, '\\titlerowcol \\textbf{Variable} & \\textbf{Datasets}\\tabularnewline\\midrule' );
            asterisks = 0;
            for jx = 1 : length( ar.model(jm).x )
                dS = arFindData( ar, jm, 'state', ar.model(jm).x{jx} );
                
                str = '';
                for jplot=1:length(ar.model(jm).plot)
                    if (ar.model(jm).qPlotYs(jplot)||~opts.omitnonplotted)
                        relevantData = intersect( dS, ar.model(jm).plot(jplot).dLink );
                        if ( length( relevantData ) > 0 ) %#ok
                            for s = 1 : length( relevantData )
                                q(s) = max( ar.model(jm).data( relevantData(s) ).qFit );
                            end
                            if ( max(q(s)) == 0 )
                                ast = '$^{*}$';
                                asterisks = 1;
                            else
                                ast = '';
                            end
                            str = sprintf( '%s \\ref{M%dExp%d}%s', str, jm, jplot, ast );
                        end
                    end
                end
                % Only include lines with data
                if ( strcmp(str, '') == 0 )
                    lp(fid, '%s%s & %s\\tabularnewline', alternate(jx), strrep(ar.model(jm).x{jx}, '_', '\_'), str );
                end
            end
            lp(fid, '\\botrule\\end{tabularx}');
            if ( asterisks )
                lp(fid, '\\caption[align=left,justification=justified,singlelinecheck=false]{Global overview of the used measurement data. Each line corresponds to a list of experiments that a particular model variable was observed in. Asterisks denote datasets that were not used for parameter optimization.}\\label{stateData}');
            else
                lp(fid, '\\caption[align=left,justification=justified,singlelinecheck=false]{Global overview of the used measurement data. Each line corresponds to a list of experiments that a particular model variable was observed in.}\\label{stateData}');
            end
            lp(fid, '\\end{table}');
        end
        
        % Cache the target replacement, since syms are expensive!
        for a = 1 : length( ar.model(jm).fp )
            fpSymString{a} = char(arSym(ar.model(jm).fp{a}));
        end

        %% do we override some of the plots
         if opts.overrideplots
            k = dir( [opts.overrideplots_args '/*.pdf'] );
            if ( ~isempty(k) )
                copyfile( [opts.overrideplots_args '/*.pdf'], [savePath '/' opts.overrideplots_args] );
                for a = 1 : length( k )
                    if ( ~strcmp( k(a).name, '.' ) && ~strcmp( k(a).name, '..' ) )
                        pdfCrop( [savePath '/' opts.overrideplots_args '/' k(a).name(1:end-4) ] )
                    end
                end
            end
            k = dir( [opts.overrideplots_args '/*.tex'] );
            if ( ~isempty(k) )
                copyfile( [opts.overrideplots_args '/*.tex'], [savePath '/' opts.overrideplots_args] );
            end
            
            for jplot=1:length(ar.model(jm).plot)
                if ( exist( [ opts.overrideplots_args '/' ar.model(jm).plot(jplot).name '.pdf' ], 'file' ) || exist( [ opts.overrideplots_args '/' ar.model(jm).plot(jplot).name '_Report.tex' ], 'file' ) )
                    ar.model(jm).plot(jplot).savePath_FigY = [opts.overrideplots_args '/' ar.model(jm).plot(jplot).name];
                    ar.model(jm).plot(jplot).override = 1;
                    fprintf( 'Overridden %s ...', ar.model(jm).plot(jplot).name );
                end
            end
        end
        for jplot=1:length(ar.model(jm).plot)
            if (ar.model(jm).qPlotYs(jplot)||~opts.omitnonplotted)
                lp(fid, '\\FloatBarrier');
                lp(fid, '\\clearpage\n');
                jd = ar.model(jm).plot(jplot).dLink(1);
                if(isfield(ar.model(jm), 'data'))
                    lp(fid, '\\subsection{Experiment: %s}\n', strrep(arNameTrafo(ar.model(jm).plot(jplot).name), '\_', ' '));
                    lp(fid, '\\label{M%dExp%d}\n', jm, jplot );
                    %% descriptions
                    if(~isempty(ar.model(jm).data(jd).description))
                        if ~strcmp( ar.model(jm).data(jd).description, 'data .def file template' )
                            lp(fid, '\\subsubsection{Comments}');
                            for jdes=1:length(ar.model(jm).data(jd).description)
                                lp(fid, '%s\\\\', strrep(strrep(ar.model(jm).data(jd).description{jdes}, '%', '\%'), '_', '\_'));
                            end
                        end
                    end

                    %% fit
                    if( ~opts.omitnonfitted || isFitted(jm, jplot) )
                        lp(fid, '\\subsubsection{Model fit and plots}');
                        
                        %% plots
                        if(isfield(ar.model(jm).plot(jplot), 'savePath_FigY') && ~isempty(ar.model(jm).plot(jplot).savePath_FigY))
                            lp(fid, 'The model observables and the experimental data is shown in Figure \\ref{%s}.', [ar.model(jm).plot(jplot).name '_y']);
                            captiontext = sprintf('\\textbf{%s observables and experimental data for the experiment.} ', arNameTrafo(ar.model(jm).plot(jplot).name));
                            if ( ~isfield( ar.model(jm).plot(jplot), 'override' ) )
                                captiontext = [captiontext 'The observables are displayed as solid lines. '];
                                captiontext = [captiontext 'The error model that describes the measurement noise ' ...
                                'is indicated by shades.'];
                            else
                                if ( isfield( ar.model(jm).plot(jplot), 'caption' ) )
                                    captiontext = [captiontext ' ' ar.model(jm).plot(jplot).caption ];
                                end
                            end
                            if(opts.texplots&&exist([ar.model(jm).plot(jplot).savePath_FigY '_Report.tex'],'file')==2)
                                copyfile([ar.model(jm).plot(jplot).savePath_FigY '_Report.tex'], ...
                                [savePath '/' ar.model(jm).plot(jplot).name '_y.tex']);
                                lpfigurePGF(fid, [ar.model(jm).plot(jplot).name '_y.tex'], captiontext, [ar.model(jm).plot(jplot).name '_y']);
                            else
                                copyfile([ar.model(jm).plot(jplot).savePath_FigY '.pdf'], ...
                                [savePath '/' ar.model(jm).plot(jplot).name '_y.pdf']);
                                if ( isfield( ar.model(jm).plot(jplot), 'nCols' ) )
                                    if ( ar.model(jm).plot(jplot).nCols == 1 )
                                        figSize = .4;
                                    elseif ( ar.model(jm).plot(jplot).nCols == 1.5 )
                                        figSize = .7;
                                    else
                                        figSize = .9;
                                    end
                                else
                                    figSize = .9;
                                end

                                lpfigure(fid, figSize, [ar.model(jm).plot(jplot).name '_y.pdf'], captiontext, [ar.model(jm).plot(jplot).name '_y']);
                            end
                        end
                    end

                    % input
                    qmod = ~strcmp(ar.model(jm).fu, ar.model(jm).data(jd).fu);
                    if ( sum(qmod) > 1 )                
                        input = sprintf('This experiment requires $%d$ custom inputs which are defined in Table \\ref{%s_input}. ', sum(qmod), ar.model(jm).plot(jplot).name );
                    elseif ( sum(qmod) == 1 )
                        input = sprintf('This experiment requires a custom input function which is defined in Table \\ref{%s_input}. ', ar.model(jm).plot(jplot).name );
                    else
                        input = '';
                    end

                    % trafo
                    [condTrans, names] = conditionSpecificParameters( jplot, jm, fpSymString, opts );
                    nUniqueCondTrafo = length(fieldnames(condTrans));
                    if ( nUniqueCondTrafo > 1 )
                        trafo = sprintf('The $%d$ necessary parameter transformations are listed in Table \\ref{%s_conditiontrafo}. Note that transformations with only one entry are the same for all experimental conditions corresponding to this experiment. ', nUniqueCondTrafo, ar.model(jm).plot(jplot).name );
                    elseif ( nUniqueCondTrafo == 1 )
                        trafo = sprintf('The necessary parameter transformation is listed in Table \\ref{%s_conditiontrafo}. ', ar.model(jm).plot(jplot).name );
                    else
                        trafo = '';
                    end

                    if ( isFitted(jm, jplot) )
                        lp(fid, '\\noindent %s%sThe experimental data is given in Table \\ref{%s_data}. ', input, trafo, ar.model(jm).plot(jplot).name);

                        % Was this plot ever simulated?
                        if(~isfield(ar.model(jm).plot(jplot), 'ndata'))
                            setchi2fields(jm,jplot);
                        end

                        
                        if(ar.config.fiterrors == 1)
                            if ( opts.omitlikelihood )
                                lp( fid, sprintf( 'This dataset contains %i data points.\\\\\\\\ \\n', ar.model(jm).plot(jplot).ndata ) );
                            else
                                lp( fid, sprintf( 'The model yields a value of $-2 \\\\log(L) = %g$ for %i data points in this data set.\\\\\\\\ \\n', 2*ar.model(jm).plot(jplot).ndata*log(sqrt(2*pi)) + ar.model(jm).plot(jplot).chi2, ar.model(jm).plot(jplot).ndata ) );
                            end
                        else
                            if ( opts.omitlikelihood )
                                lp( fid, sprintf( 'This dataset contains %i data points.\\\\\\\\ \\n', ar.model(jm).plot(jplot).ndata ) );
                            else
                                lp( fid, sprintf( 'The model yields a value of $\\\\chi^2 = %g$ for %i data points in this data set.\\\\\\\\ \\n', ar.model(jm).plot(jplot).chi2, ar.model(jm).plot(jplot).ndata ) );
                            end        
                        end
                    end

                    %% trajectories
                    if(isfield(ar.model(jm).plot(jplot), 'savePath_FigX') && ~isempty(ar.model(jm).plot(jplot).savePath_FigX))
                        lp(fid, ['The trajectories of the input, dynamic and derived variables that ' ...
                            'correspond to the experimental conditions in this experiment are shown in Figure \\ref{%s}.'], ...
                            [ar.model(jm).plot(jplot).name '_x']);
                        copyfile([ar.model(jm).plot(jplot).savePath_FigX '.pdf'], ...
                            [savePath '/' ar.model(jm).plot(jplot).name '_x.pdf'])

                        captiontext = sprintf('\\textbf{%s trajectories of the input, dynamic and derived variables.} ', ....
                            arNameTrafo(ar.model(jm).plot(jplot).name));
                        captiontext = [captiontext 'The dynamical behaviour is determined by the ODE system, see Equation '];
                        captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \ref{' ...
                            sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
                        lpfigure(fid, .9, [ar.model(jm).plot(jplot).name '_x.pdf'], captiontext, [ar.model(jm).plot(jplot).name '_x']);
                    end
                    if(isfield(ar.model(jm).plot(jplot), 'savePath_FigV') && ~isempty(ar.model(jm).plot(jplot).savePath_FigV))
                        lp(fid, 'The reaction fluxes that correspond to the experimental conditions in this experiment are shown in Figure \\ref{%s}.', ...
                            [ar.model(jm).plot(jplot).name '_v']);
                        copyfile([ar.model(jm).plot(jplot).savePath_FigV '.pdf'], ...
                            [savePath '/' ar.model(jm).plot(jplot).name '_v.pdf'])

                        captiontext = sprintf('\\textbf{%s reaction fluxes.} ', arNameTrafo(ar.model(jm).plot(jplot).name));
                        captiontext = [captiontext 'The dynamical behaviour is determined by the ODE system, see Equation '];
                        captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) ...
                            '} -- \ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
                        lpfigure(fid, .9, [ar.model(jm).plot(jplot).name '_v.pdf'], captiontext, [ar.model(jm).plot(jplot).name '_v']);
                    end            

                    %% inputs
                    if(~isempty(ar.model(jm).u))
                        % Which ones changed
                        qmod = ~strcmp(ar.model(jm).fu, ar.model(jm).data(jd).fu);

                        % Did the user fill in better descriptions for the states?
                        if (isfield(ar.model(jm).data(jd),'uNames'))&&(~isempty(ar.model(jm).data(jd).uNames))&&(max( cellfun(@length,ar.model(jm).data(jd).uNames) ) > 0)
                            headerChunk = '& \textbf{Description}';
                            nC = 4;
                        else
                            headerChunk = '';
                            nC = 3;
                        end      

                        if ( sum(qmod) > 0 )
                            %lp(fid, '\\begin{table}');
                            %lp(fid, '\t\\centering');
                            lp(fid, '\\begin{statictable}');
                            %lp(fid, '\t\\centering');
                            lp(fid, '\\begin{centering}');
                            lp(fid, '\\begin{tabularx}{\\textwidth}{@{} *%dlX @{}}\\toprule', nC-1);
                            sk = find(qmod);
                            lp(fid, '\\titlerowcol \\textbf{Input} & \\textbf{Unit} & \\textbf{Modified equation} %s\\tabularnewline\\midrule', headerChunk );
                            if ( strcmp( headerChunk, '' ) )
                                for ku = 1 : length( sk )
                                    ju = sk(ku);
                                    inp = repFunc( ar.model(jm).data(jd).fu{ju}, 'monospline10' );
                                    newl = sprintf( '$\\tabularnewline %s & & $', alternate(ku) );
                                    lp(fid, '%s%s & %s [%s] & $%s$\\tabularnewline', alternate(ku), strrep(ar.model(jm).u{ju}, '_', '\_'), ar.model(jm).uUnits{ju,3}, ar.model(jm).uUnits{ju,2}, funcWrap( myFormulas(inp, jm), opts.funcwrap, newl ) );
                                end                            
                            else
                                for ku = 1 : length( sk )
                                    ju = sk(ku);
                                    inp = repFunc( ar.model(jm).data(jd).fu{ju}, 'monospline10' );
                                    newl = sprintf( '$\\tabularnewline %s & & & $', alternate(ku) );
                                    lp(fid, '%s%s & %s [%s] & $%s$ & %s\\tabularnewline', alternate(ku), strrep(ar.model(jm).u{ju}, '_', '\_'), ar.model(jm).uUnits{ju,3}, ar.model(jm).uUnits{ju,2}, funcWrap( myFormulas(inp, jm), opts.funcwrap, newl ), ar.model(jm).data(jd).uNames{ku} );
                                end
                            end
                            lp(fid, '\\botrule\\end{tabularx}');
                            %lp(fid, '\\mycaption{Inputs modified for experiment %s}{%s_input}{}', arNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name );
                            %lp(fid, '\\end{table}');
                            lp(fid, '\\mycaptionof{Inputs modified for experiment %s}{%s_input}{}', arNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name );
                            lp(fid, '\\end{centering}');
                            lp(fid, '\\end{statictable}');       
                        end
                    end

                    %% observations and error model
                    if(isfield(ar.model(jm), 'y'))
                        qadd = ~ismember(ar.model(jm).data(jd).y, ar.model(jm).y);
                    else
                        qadd = 0;
                    end

                    if(isfield(ar.model(jm), 'y'))
                        fytmp = strrep(ar.model(jm).fy, '_filename', ['_' ar.model(jm).data(jd).name]);
                        fystdtmp = strrep(ar.model(jm).fystd, '_filename', ['_' ar.model(jm).data(jd).name]);
                        qmod = [];
                        for j=1:length(ar.model(jm).data(jd).fy)
                            if(~qadd(j))
                                iy = find(strcmp(ar.model(jm).y, ar.model(jm).data(jd).y{j}));
                                qmod(j) = ~strcmp(fytmp{iy}, ar.model(jm).data(jd).fy{j}) || ...
                                    ~strcmp(fystdtmp{iy}, ar.model(jm).data(jd).fystd{j});
                            end
                        end
                    else
                        qmod = false;
                    end

                    if(sum(qmod)>0 || sum(qadd)>0)
                        lp(fid, '\\subsubsection{Observables}\n');

                        % modified observables
                        if(sum(qmod)>0)
                            lp(fid, '\\noindent The following observables are modified in this data set:\\\\\n');

                            lp(fid, '\\begin{statictable}');
                            lp(fid, '\\begin{centering}');
                            box = sprintf('%sobsmod',latexIdentifier(jplot)) ;
                            startFlexbox(fid,  box );
                            lp(fid, '\\begin{tabular}{@{} lll @{}}\\toprule');
                            lp(fid, '\\titlerowcol \\textbf{Observable} & & \\textbf{Equations} \\tabularnewline\\midrule' );                        

                            sk = find(qmod);
                            for ky = 1 : length( sk )
                                jy = sk(ky);
                                strtmp = myFormulas(PTI(jm, ar.model(jm).data(jd).fy{jy}, pti), jm);
                                if(isfield(ar.model(jm).data(jd),'logfitting') && ar.model(jm).data(jd).logfitting(jy))
                                    strtmp = ['\mathrm{log}_{10}(' strtmp ')'];
                                end            
                                if ( strcmp( ar.model(jm).data(jd).yNames{jy}, '' ) )
                                    lp(fid, ' %s & y & $%s$ \\tabularnewline %s \\multirow{-2}{*}{%s [%s]} & $\\sigma$ & $%s$ \\tabularnewline', alternate(ky), strtmp, alternate(jy), strrep(ar.model(jm).data(jd).y{jy}, '_', '\_\-'), arNameTrafo( ar.model(jm).yUnits{jy,2} ), myFormulas(ar.model(jm).data(jd).fystd{jy}, jm) );
                                else
                                    lp(fid, ' %s & y & $%s$ \\tabularnewline %s \\multirow{-2}{*}{%s [%s]} & $\\sigma$ & $%s$ \\tabularnewline', alternate(ky), strtmp, alternate(jy), strrep(ar.model(jm).data(jd).yNames{jy}, '_', '\_\-'), arNameTrafo( ar.model(jm).yUnits{jy,2} ), myFormulas(ar.model(jm).data(jd).fystd{jy}, jm) );
                                end
                            end
                            lp(fid, '\\botrule\\end{tabular}');
                            endFlexbox(fid, box);
                            lp(fid, '\\mycaptionof{Observables modified for experiment %s}{%s_input}{}', arNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name );
                            lp(fid, '\t\\end{centering}');
                            lp(fid, '\\end{statictable}\n');                          
                        end

                        % added observables
                        if(sum(qadd)>0)
                            lp(fid, '\\noindent The following observables are added in this data set:\\\\\n');

                            lp(fid, '\\begin{statictable}');
                            lp(fid, '\\begin{centering}');
                            box = sprintf('%sobsadd',latexIdentifier(jplot)) ;
                            startFlexbox(fid,  box);
                            lp(fid, '\\begin{tabular}{@{} lll @{}}\\toprule');
                            lp(fid, '\\titlerowcol \\textbf{Observable} & & \\textbf{Equations} \\tabularnewline\\midrule' );                        

                            sk = find(qadd);
                            for ky = 1 : length( sk )
                                jy = sk(ky);
                                strtmp = myFormulas(PTI(jm, ar.model(jm).data(jd).fy{jy}, pti), jm);
                                if(isfield(ar.model(jm).data(jd),'logfitting') && ar.model(jm).data(jd).logfitting(jy))
                                    strtmp = ['\mathrm{log}_{10}(' strtmp ')'];
                                end

                                if ( strcmp( ar.model(jm).data(jd).yNames{jy}, '' ) )
                                    lp(fid, ' %s & y & $%s$ \\tabularnewline %s \\multirow{-2}{*}{%s [%s]} & $\\sigma$ & $%s$ \\tabularnewline', alternate(ky), strtmp, alternate(jy), strrep(ar.model(jm).data(jd).y{jy}, '_', '\_\-'), arNameTrafo( ar.model(jm).data(jd).yUnits{jy,2} ), myFormulas(ar.model(jm).data(jd).fystd{jy}, jm) );
                                else
                                    lp(fid, ' %s & y & $%s$ \\tabularnewline %s \\multirow{-2}{*}{%s [%s]} & $\\sigma$ & $%s$ \\tabularnewline', alternate(ky), strtmp, alternate(jy), strrep(ar.model(jm).data(jd).yNames{jy}, '_', '\_\-'), arNameTrafo( ar.model(jm).data(jd).yUnits{jy,2} ), myFormulas(ar.model(jm).data(jd).fystd{jy}, jm) );
                                end
                            end
                            lp(fid, '\\botrule\\end{tabular}');
                            endFlexbox(fid, box);
                            lp(fid, '\\mycaptionof{Observables added for experiment %s}{%s_input}{}', arNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name );
                            lp(fid, '\t\\end{centering}');
                            lp(fid, '\\end{statictable}\n');                          
                        end                    
                    end

                    %% conditions
                    vars        = fieldnames(condTrans);
                    str         = cell(0);
                    row         = cell(0);
                    var         = cell(0);
                    if ( length( vars ) > 0 ) %#ok

                        cols = '';
                        for q = 1 : length(condTrans.(vars{1}))
                            cols = [ cols 'K'];
                        end

                        q = 1;
                        for jv = 1 : length( vars )
                            str{q} = '';
                            row{q} = '';
                            var{q} = '';
                            if (length( unique(condTrans.(vars{jv})) ) == 1 )
                                variableName = strrep(PTI(jm, names.(vars{jv}), pti),'_','\_');
                                PTI(jm, names.(vars{jv}),pti)
                                str{q} = variableName;
                                formula = myFormulas(PTI(jm, condTrans.(vars{jv}){1}, pti), jm);
                                str{q} = sprintf('%s & \\multicolumn{%d}{c}{$%s$}', str{q}, length(condTrans.(vars{1})), formula );
                                var{q} = variableName;
                                row{q} = strcat( row{q}, formula, 'Q Q ' );
                                q = q + 1;
                            end
                        end
                        for jv = 1 : length( vars )
                            str{q} = '';
                            row{q} = '';    
                            var{q} = '';
                            if ~(length( unique(condTrans.(vars{jv})) ) == 1 )
                                variableName = strrep(PTI(jm, names.(vars{jv}),pti),'_','\_');
                                PTI(jm, names.(vars{jv}),pti)
                                str{q} = variableName;
                                var{q} = variableName;
                                for jdls = 1 : length(ar.model(jm).plot(jplot).dLink)
                                    formula = myFormulas(PTI(jm, condTrans.(vars{jv}){jdls}, pti), jm);
                                    str{q} = sprintf( '%s & $%s$', str{q}, formula );
                                    row{q} = strcat( row{q}, formula, 'Q Q ' );
                                end
                                q = q + 1;
                            end
                        end                           
                        row{q} = 'Q Q Condition values';
                        var{q} = 'Q Q Parameter';
                        lp(fid, '\\subsubsection{Condition dependent parameter changes}\n The following model parameters were changed to simulate these experimental conditions:\\\\\n');
                        box = sprintf('%schnk',latexIdentifier(jplot));
                        lp(fid, '\\begin{statictable}\\');
                        lp(fid, '\t\\centering');
                        startFlexbox(fid, box);
                        
                        % Have latex figure out the longest strings in each column
                        maxV = sprintf('\\widthof{$%s$Q Q}', var{1});
                        maxR = sprintf('\\widthof{$%s$}', row{1});
                        for a = 2 : length( var )
                            maxV = sprintf( '\\maxof{\n\t\\widthof{$%s$}}{%s}',var{a}, maxV );
                            maxR = sprintf( '\\maxof{\n\t\\widthof{$%s$}}{%s}',row{a}, maxR );
                        end

                        % TabularX is nicer but will fail with more than 19 columns
                        if ( length( cols ) < 20 )
                            head = sprintf( '{\\footnotesize \\begin{tabularx}{%s+%s}{c| %s} \\toprule', maxV, maxR, cols ); %@{.}
                            tail = sprintf( '\\botrule\\end{tabularx} }' ); 
                        else
                            head = sprintf( '{\\footnotesize \\begin{tabular}{c| %s} \\toprule', strrep(cols, 'K', 'c') ); %@{.}
                            tail = sprintf( '\\botrule\\end{tabular} }' );                         
                        end

                        % Alright, we assembled all the data, start printing
                        % this thing
                        lp(fid, '\t%s', head);
                        lp(fid, '\\titlerowcol \\textbf{Parameter} & \\multicolumn{%d}{c}{\\textbf{Condition values}} \\tabularnewline\\midrule', length(condTrans.(vars{1})) );
                        for a = 1 : length( str )
                            if ( ~isempty(str{a}) )
                                lp(fid, '\t\t%s\\tabularnewline', str{a});
                            end
                        end
                        lp(fid, '\t%s', tail);
                        endFlexbox(fid, box);

                        lp(fid, '\t\\mycaptionof{Model parameters modified for experiment %s. Different rows indicate different conditions.}{%s_conditiontrafo}{}', strrep(arNameTrafo(ar.model(jm).plot(jplot).name), '\_', ' '), ar.model(jm).plot(jplot).name );
                        lp(fid, '\\end{statictable}\\');
                    end

                    if( ~opts.omitnonfitted || isFitted(jm, jplot) )
                        %% experimental data
                        headstr = '\titlerowcol ';
                        headtab = '';
                        unitstr = '\titlerowcol ';
                        % time
                        unitstr = [unitstr sprintf(' %s [%s] ', arNameTrafo(ar.model(jm).data(jd).tUnits{3}), ...
                            arNameTrafo(ar.model(jm).data(jd).tUnits{2}))];
                        headstr = [headstr ' '];
                        headtab = [headtab 'r'];
                        % conditions
                        if(~isempty(ar.model(jm).data(jd).condition))
                            for jp = 1:length(ar.model(jm).data(jd).condition)
                                headstr = [headstr sprintf('& %s ', arNameTrafo(ar.model(jm).data(jd).condition(jp).parameter))];
                                unitstr = [unitstr '& '];
                                headtab = [headtab 'r'];
                            end
                        end
                        % y & ystd headers
                        for jy=1:length(ar.model(jm).data(jd).y)
                            headstr = [headstr sprintf('& %s ', arNameTrafo(ar.model(jm).data(jd).y{jy}))];
                            unitstr = [unitstr sprintf('& %s [%s] ', arNameTrafo(ar.model(jm).data(jd).yUnits{jy,3}), ...
                                arNameTrafo(ar.model(jm).data(jd).yUnits{jy,2}))];
                            headtab = [headtab 'r'];
                            if(ar.config.fiterrors == -1)
                                headstr = [headstr sprintf('& %s\\_std ', arNameTrafo(ar.model(jm).data(jd).y{jy}))];
                                unitstr = [unitstr sprintf('& %s [%s] ', arNameTrafo(ar.model(jm).data(jd).yUnits{jy,3}), ...
                                    arNameTrafo(ar.model(jm).data(jd).yUnits{jy,2}))];
                                headtab = [headtab 'r'];
                            end
                        end
                           
                        Nd = 55;
                        L = 0;
                        
                        % Avoid ending up with a last table with less than 4 entries, since
                        % this looks silly.
                        overhang = Nd-round((length(ar.p)/Nd))*Nd;
                        if ( overhang < 4 )
                            Nd = 50;
                        end
                        
                        openDataTable( fid, jm, jd, L, headtab, headstr, unitstr );

                        s = 0;

                        NAs = false;
                        for jd2 = ar.model(jm).plot(jplot).dLink
                            NAs = NAs || sum(isnan(ar.model(jm).data(jd2).yExpStd(:)))>0;
                            for j=1:length(ar.model(jm).data(jd2).tExp)
                                fprintf(fid, '%s%s ', alternate(s), sprintf('%g', ar.model(jm).data(jd2).tExp(j)));

                                % conditions
                                if(~isempty(ar.model(jm).data(jd2).condition))
                                    for jp = 1:length(ar.model(jm).data(jd2).condition)
                                        fprintf(fid, '& %s ', strrep(ar.model(jm).data(jd2).condition(jp).value,'_','\_'));
                                    end
                                end

                                % y data
                                for jj=1:size(ar.model(jm).data(jd2).yExp,2)
                                    if(ar.model(jm).data(jd2).logfitting(jj))
                                        fprintnumtab(fid, 10^ar.model(jm).data(jd2).yExp(j,jj));
                                    else
                                        fprintnumtab(fid, ar.model(jm).data(jd2).yExp(j,jj));
                                    end
                                end

                                % ystd data
                                if(ar.config.fiterrors == -1)
                                    for jj=1:size(ar.model(jm).data(jd2).yExp,2)
                                        fprintnumtab(fid, ar.model(jm).data(jd2).yExpStd(j,jj));
                                    end
                                end
                                fprintf(fid, '\\\\\n'); 
                                                                
                                % Are we going over the page limit?
                                s = s + 1;
                                if ( s > Nd )
                                    closeDataTable( fid, jplot, jm, jd, L, headtab, NAs );
                                    L = L + 1;
                                    s = 0;
                                    openDataTable( fid, jm, jd, L, headtab, headstr, unitstr );
                                end
                            end
                        end

                        closeDataTable( fid, jplot, jm, jd, L, headtab, NAs );
                    end

             %       if(ccount>1 && jp==length(ar.model(jm).data(jd).fp))
             %           lp(fid, '\\end{array}');
             %           lp(fid, '\\end{displaymath}\n}\n\n');
             %       end
                end
            end
        end
    end
end

%% Parameters
lp(fid, '\\FloatBarrier');
lp(fid, '\\clearpage\n');
lp(fid, '\\section{Estimated model parameters} \\label{estimatedparameters}\n');
lp(fid, 'This section lists all the model parameters used in the model.');
lp(fid, 'Parameters highlighted in red color indicate parameter values close to their bounds.');
if ( ~isempty( cell2mat( strfind( ar.pLabel, 'init_' ) ) ) ) 
    lp(fid, 'The parameter name prefix init\\_ indicates the initial value of a dynamic variable.');
end
if ( ~isempty( cell2mat( strfind( ar.pLabel, 'offset_' ) ) ) ) 
    lp(fid, 'The parameter name prefix offset\\_ indicates a offset of the experimental data.');
end
if ( ~isempty( cell2mat( strfind( ar.pLabel, 'scale_' ) ) ) ) 
    lp(fid, 'The parameter name prefix scale\\_ indicates a scaling factor of the experimental data.');
end
if ( ~isempty( cell2mat( strfind( ar.pLabel, 'sd_' ) ) ) ) 
    lp(fid, 'The parameter name prefix sd\\_ indicates the magnitude of the measurement noise for a specific measurement.\\\\');
end

N = 50;
% Avoid ending up with a last table with less than 4 entries, since
% this looks silly.
overhang = N-round((length(ar.p)/N))*N;
if ( overhang < 4 )
    N = 45;
end

lp(fid, '\t\\begin{table}');
lp(fid, '\t\\dobegincenter');
startFlexbox(fid, sprintf('partable%s', latexIdentifier(0)) );
lp(fid, '\t{\\footnotesize');
lp(fid, '\t\t\\begin{tabular}{llllllll}');
lp(fid, '\t\t\t\\toprule');
lp(fid, '\t\t\t \\titlerowcol & name & $\\theta_{min}$ & $\\hat\\theta$ & $\\theta_{max}$ & log & non-log $\\hat \\theta$ & fitted \\\\');
lp(fid, '\t\t\t\\midrule');
count = 1;
for j=1:length(ar.p)
    if(mod(j,N)==0)
        lp(fid, '\t\t\t\\botrule');
        lp(fid, '\t\t\\end{tabular}}');
        endFlexbox(fid, sprintf('partable%s',latexIdentifier(count-1)));
        lp(fid, '\t\t\\mycaption{Estimated parameter values}{paratable%i}', count);
        lp(fid, '{$\\hat \\theta$ indicates the estimated value of the parameters.');
        lp(fid, '$\\theta_{min}$ and $\\theta_{max}$ indicate the upper and lower bounds for the parameters.');
        lp(fid, 'The log-column indicates if the value of a parameter was log-transformed.');
        lp(fid, 'If log = 1 the non-log-column indicates the non-logarithmic value of the estimate.');
        lp(fid, 'The fitted-column indicates if the parameter value was estimated (1), was temporarily fixed (0) or if its value was fixed to a constant value (2).}');
        lp(fid, '\t\\doendcenter');
        lp(fid, '\t\\end{table}');
        lp(fid, '\t\\begin{table}');
        lp(fid, '\t\\dobegincenter');
        startFlexbox(fid, sprintf('partable%s', latexIdentifier(count)) );
        lp(fid, '\t{\\footnotesize');
        lp(fid, '\t\t\\begin{tabular}{llllllll}');
        lp(fid, '\t\t\t\\toprule');
        lp(fid, '\t\t\t \\titlerowcol & name & $\\theta_{min}$ & $\\hat\\theta$ & $\\theta_{max}$ & log & non-log $\\hat \\theta$ & fitted \\\\');
        lp(fid, '\t\t\t\\midrule');
        
        count = count + 1;
    end
    
    if(ar.qFit(j)==1)
        if(ar.p(j) - ar.lb(j) < 0.1 || ar.ub(j) - ar.p(j) < 0.1)
            lp(fid, '%s\t\t\t\\color{red}{%i} & \\color{red}{%s} & \\color{red}{%+8.4g} & \\color{red}{%+8.4f} & \\color{red}{%+8.4g} & \\color{red}{%i} & \\color{red}{%s} & \\color{red}{%i} \\\\', ...
                alternate(j), j, strrep(PTI(jm, ar.pLabel{j},pti),'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
                ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
        else
            lp(fid, '%s\t\t\t%i & %s & {%+8.4g} & {%+8.4f} & {%+8.4g} & %i & %s & %i \\\\', ...
                alternate(j), j, strrep(PTI(jm, ar.pLabel{j},pti),'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
                ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
        end
    else
        lp(fid, '%s\t\t\t\\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%+8.4g} & \\color{mygray}{%+8.4f} & \\color{mygray}{%+8.4g} & \\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%i} \\\\', ...
            alternate(j), j, strrep(PTI(jm, ar.pLabel{j},pti),'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
            ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
    end
end
lp(fid, '\t\t\t\\botrule');
lp(fid, '\t\t\\end{tabular}}');
endFlexbox(fid, sprintf('partable%s', latexIdentifier(count-1)) );

lp(fid, '\t\t\\mycaption{Estimated parameter values}{paratable%i}', count);
lp(fid, '{$\\hat \\theta$ indicates the estimated value of the parameters.');
lp(fid, '$\\theta_{min}$ and $\\theta_{max}$ indicate the upper and lower bounds for the parameters.');
lp(fid, 'The log-column indicates if the value of a parameter was log-transformed.');
lp(fid, 'If log = 1 the non-log-column indicates the non-logarithmic value of the estimate.');
lp(fid, 'The fitted-column indicates if the parameter value was estimated (1), was temporarily fixed (0) or if its value was fixed to a constant value (2).}');
lp(fid, '\t\\doendcenter');
lp(fid, '\t\\end{table}');

% Disable PLE inclusion
if 0

%% PLE
plePath = [arSave '/PLE'];
if(exist(plePath,'dir'))
    lp(fid, '\\clearpage\n');
    lp(fid, '\\section{Profile likelihood of model parameters}\n');
    
    S = load([plePath '/results.mat']);
    
    lp(fid, 'In order to evaluate the identifiability of the model parameters and to assess confidence intervals, ');
    lp(fid, 'the profile likelihood \\cite{Raue:2009ec} was calculated.');
    lp(fid, 'The mean calculation time of the profile likelihood per parameter was %s $\\pm$ %s.', ...
        secToHMS(mean(S.ar.ple.timing(S.ar.ple.q_fit(size(S.ar.ple.timing))))), ...
        secToHMS(std(S.ar.ple.timing(S.ar.ple.q_fit(size(S.ar.ple.timing))))));
    
    % Multiplot
    if(isfield(S.ar.ple, 'fighandel_multi'))
        lp(fid, 'An overview is displayed in Figure \\ref{multi_plot}.');
        
        sourcestr = [plePath '/multi_plot.eps'];
        targetstr = [savePath '/multi_plot.pdf'];
        if(ispc)
            print('-dpdf', targetstr);
        elseif(ismac)
            system(['/usr/local/bin/ps2pdf  -dEPSCrop ' sourcestr ' ' targetstr]);
        else
            system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' sourcestr ' ' targetstr]);
        end
        
        captiontext = '\textbf{Overview of the profile likelihood of the model parameters}\\';
        captiontext = [captiontext 'The solid lines indicate the profile likelihood. '];
        captiontext = [captiontext 'The broken lines indicate the threshold to assess confidence intervals. '];
        captiontext = [captiontext 'The asterisks indicate the optimal parameter values. '];
        lpfigure(fid, 1, 'multi_plot.pdf', captiontext, 'multi_plot');
    end
    
    % Singleplots
    if(isfield(S.ar.ple, 'figPath'))
        count = 0;
        for j=1:length(S.ar.ple.figPath)
            if(~isempty(S.ar.ple.chi2s{j}))
                count = count + 1;
            end
        end
        
        lp(fid, 'In Figure \\ref{ple1} -- \\ref{ple%i} the profile likelihood of each parameter is shown in more detail.', count);
        lp(fid, 'Also the functional relations to the remaining parameter are displayed.');
        
        count = 0;
        for j=1:length(S.ar.ple.figPath)
            if(~isempty(S.ar.ple.chi2s{j}))
                lp(fid, '\\clearpage\n');
                count = count + 1;
                targetstr = [savePath '/' S.ar.ple.p_labels{j} '.pdf'];
                if(ispc)
                    print('-dpdf', targetstr);
                elseif(ismac)
                    system(['/usr/local/bin/ps2pdf  -dEPSCrop ' S.ar.ple.figPath{j} '.eps ' targetstr]);
                else
                    system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' S.ar.ple.figPath{j} '.eps ' targetstr]);
                end
                
                captiontext = sprintf('\\textbf{Profile likelihood of parameter %s}\\\\', strrep(S.ar.ple.p_labels{j},'_','\_'));
                captiontext = [captiontext 'Upper panel: The solide line indicates the profile likelihood. '];
                captiontext = [captiontext 'The broken line indicates the threshold to assess confidence intervals. '];
                captiontext = [captiontext 'The asterisk indicate the optimal parameter values. '];
                captiontext = [captiontext sprintf('Lower panel: The functional relations to the other parameters along the profile likelihood of %s are displayed. ', strrep(S.ar.ple.p_labels{j},'_','\_'))];
                captiontext = [captiontext 'In the legend the top five parameters showing the strongest variations are given. '];
                captiontext = [captiontext sprintf('The calculation time was %s.', secToHMS(S.ar.ple.timing(j)))];
                lpfigure(fid, 1, [S.ar.ple.p_labels{j} '.pdf'], captiontext, sprintf('ple%i',count));
                
                lp(fid, '\n');
            end
        end
    end
    
    %% Confidence Intervals
    lp(fid, '\\clearpage\n');
    lp(fid, '\\section{Confidence intervals for the model parameters}\n');
    
    N = 30;
    ntables = 0;
    for j=1:length(S.ar.ple.p_labels)
        if(S.ar.ple.q_fit(j))
            if(j<=length(S.ar.ple.chi2s) && ~isempty(S.ar.ple.chi2s{j}))
                ntables = ntables + 1;
            end
        end
    end
    ntables = ceil(ntables/N);
    
    if(ntables>1)
        lp(fid, 'In Table \\ref{conftable1} -- \\ref{conftable%i}, %2i\\%% confidence intervals for the estimated parameter values derived by the profile likelihood \\cite{Raue:2009ec} are given.', ntables, (1-S.ar.ple.alpha_level)*100);
    else
        lp(fid, 'In Table \\ref{conftable1}, %2i\\%% confidence intervals for the estimated parameter values derived by the profile likelihood \\cite{Raue:2009ec} are given.', (1-S.ar.ple.alpha_level)*100);
    end
    
    headstr = '\t\t\t & name & $\\hat\\theta$';
    if(S.ar.ple.plot_point && S.ar.ple.plot_simu)
        headstr = [headstr '& $\\sigma^{-}_{ptw}$ & $\\sigma^{+}_{ptw}$'];
        headstr = [headstr '& $\\sigma^{-}_{sim}$ & $\\sigma^{+}_{sim}$'];
    else
        headstr = [headstr '& $\\sigma^{-}$ & $\\sigma^{+}$'];
    end
    headstr = [headstr ' \\\\'];
    
    lp(fid, '\t\\begin{table}');
    lp(fid, '\t\\dobegincenter');
    lp(fid, '\t{\\footnotesize');
    lp(fid, '\t\t\\begin{tabular}{lllllll}');
    lp(fid, '\t\t\t\\toprule');
    lp(fid, headstr);
    lp(fid, '\t\t\t\\midrule');
    
    count = 0;
    counttab = 1;
    for j=1:length(S.ar.ple.p_labels)
        if(S.ar.ple.q_fit(j))
            if(j<=length(S.ar.ple.chi2s) && ~isempty(S.ar.ple.chi2s{j}))
                count = count + 1;
                
                if(mod(count,N)==0)
                    lp(fid, '\t\t\t\\botrule');
                    lp(fid, '\t\t\\end{tabular}}');
                    lp(fid, '\t\t\\mycaption{Confidence intervals for the estimated parameter values derived by the profile likelihood}{conftable%i}', counttab);
                    lp(fid, '\t\t{$\\hat\\theta$ indicates the estimated optimal parameter value.');
                    if(S.ar.ple.plot_point && S.ar.ple.plot_simu)
                        lp(fid, '\t\t$\\sigma^{-}_{ptw}$ and $\\sigma^{+}_{ptw}$ indicate %i\\%% point-wise confidence intervals.', (1-S.ar.ple.alpha_level)*100);
                        lp(fid, '\t\t$\\sigma^{-}_{sim}$ and $\\sigma^{+}_{sim}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.ar.ple.alpha_level)*100);
                    elseif(S.ar.ple.plot_point && ~S.ar.ple.plot_simu)
                        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% point-wise confidence intervals.', (1-S.ar.ple.alpha_level)*100);
                    elseif(~S.ar.ple.plot_point && S.ar.ple.plot_simu)
                        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.ar.ple.alpha_level)*100);
                    end
                    lp(fid, '}');
                    lp(fid, '\t\\doendcenter');
                    lp(fid, '\t\\end{table}');
                    counttab = counttab + 1;
                    
                    lp(fid, '\\clearpage\\n');
                    lp(fid, '\t\\begin{table}');
                    lp(fid, '\t\\dobegincenter');
                    lp(fid, '\t{\\footnotesize');
                    lp(fid, '\t\t\\begin{tabular}{lllllll}');
                    lp(fid, '\t\t\t\\toprule');
                    lp(fid, headstr);
                    lp(fid, '\t\t\t\\midrule');
                end
                
                lp(fid, '\t\t\t%i & %s & %+8.3f & ', j, ...
                    strrep(S.ar.ple.p_labels{j},'_','\_'), S.ar.ple.p(j));
                if(S.ar.ple.plot_point)
                    lp(fid, '%+8.3f & %+8.3f &', S.ar.ple.conf_lb_point(j), S.ar.ple.conf_ub_point(j));
                end
                if(S.ar.ple.plot_simu)
                    lp(fid, '%+8.3f & %+8.3f', S.ar.ple.conf_lb(j), S.ar.ple.conf_ub(j));
                end
                lp(fid, ' \\\\');
            end
        end
    end
    lp(fid, '\t\t\t\\botrule');
    lp(fid, '\t\t\\end{tabular}}');
    lp(fid, '\t\t\\mycaption{Confidence intervals for the estimated parameter values derived by the profile likelihood}{conftable%i}', counttab);
    lp(fid, '\t\t{$\\hat\\theta$ indicates the estimated optimal parameter value.');
    if(S.ar.ple.plot_point && S.ar.ple.plot_simu)
        lp(fid, '\t\t$\\sigma^{-}_{ptw}$ and $\\sigma^{+}_{ptw}$ indicate %i\\%% point-wise confidence intervals.', (1-S.ar.ple.alpha_level)*100);
        lp(fid, '\t\t$\\sigma^{-}_{sim}$ and $\\sigma^{+}_{sim}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.ar.ple.alpha_level)*100);
    elseif(S.ar.ple.plot_point && ~S.ar.ple.plot_simu)
        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% point-wise confidence intervals.', (1-S.ar.ple.alpha_level)*100);
    elseif(~S.ar.ple.plot_point && S.ar.ple.plot_simu)
        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.ar.ple.alpha_level)*100);
    end
    lp(fid, '}');
    lp(fid, '\t\\doendcenter');
    lp(fid, '\t\\end{table}');
    
    %% Confidence intervals of model trajectories
    %     \subsection{Confidence intervals of the predicted model dynamics} \label{obsanalysis}
    % TODO
end


end

lp(fid, '\\bibliographystyle{plain}');
lp(fid, '\\bibliography{lib}');

lp(fid, '\\end{document}');

fclose(fid);
fprintf('done\n');

%% pdflatex
if(isunix && ~ismac)
    fprintf('pdflatex, file %s...', fname);
    cd(savePath);
    eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
    eval(['!bibtex ' fnamebib ' > log_bibtex.txt']);
    eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
    eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
    cd('../../..');
    try
        copyfile([savePath '/' 'MiniReport.pdf'], [savePath '/' sprintf('report_%s.pdf', datestr(now,30))])
        fprintf('done\n');
    catch
        fprintf('MiniReport.pdf was not written correctly\n');
    end
elseif(ismac)
    fprintf('pdflatex, file %s...', fname);
    cd(savePath);
    eval(['!/usr/texbin/pdflatex ' fname ' > log_pdflatex.txt']);
    eval(['!/usr/texbin/bibtex ' fnamebib ' > log_bibtex.txt']);
    eval(['!/usr/texbin/pdflatex ' fname ' > log_pdflatex.txt']);
    eval(['!/usr/texbin/pdflatex ' fname ' > log_pdflatex.txt']);
    cd('../../..');
    try
        copyfile([savePath '/' 'MiniReport.pdf'], [savePath '/' sprintf('report_%s.pdf', datestr(now,30))])
        fprintf('done\n');
    catch
        fprintf('MiniReport.pdf was not written correctly\n');
    end
end
setenv('LD_LIBRARY_PATH', library_path);

function lp(varargin)
if(nargin>2)
    fprintf(varargin{1}, sprintf('%s\n', varargin{2}), varargin{3:end});
else
    fprintf(varargin{1}, sprintf('%s\n', varargin{2}));
end

function pdfCrop(fn)
    system(sprintf('pdfcrop -hires %s.pdf %s_crop.pdf', fn, fn));
    delete('tmp-pdfcrop-*.tex');

function lpfigure(fid, textwidth, figpath, figcaption, figlabel)

global ar;
origFile = figpath(1:end-4);
fpath = [ ar.config.savepath '/Latex/' ];
system(sprintf('pdfcrop -hires %s%s.pdf %s%s_crop.pdf', fpath, origFile, fpath, origFile));
figpath = [origFile '_crop.pdf'];

lp(fid, '\\begin{figure}[H]');
%lp(fid, '\\\\noindent \\begin{minipage}{\\textwidth}');
lp(fid, '\\begin{center}');
lp(fid, '\\includegraphics[width=%f\\textwidth]{%s} \\caption{%s} \\label{%s}', textwidth, figpath, figcaption, figlabel);
lp(fid, '\\end{center}');
%lp(fid, '\\end{minipage}');
lp(fid, '\\end{figure}');

function lpfigurePGF(fid, figpath, figcaption, figlabel)
lp(fid, '\\begin{figure}[!hb]');
lp(fid, '\\centering');
lp(fid, '\\input{%s} \\caption{%s} \\label{%s}', figpath, figcaption, figlabel);
lp(fid, '\\end{figure}')

function hmstimestr = secToHMS(seconds)
hours = floor(seconds/3600);
seconds = seconds - hours*3600;
minutes = floor(seconds/60);
seconds = seconds - minutes*60;
hmstimestr = sprintf('%02i:%02i:%05.2f', hours, minutes, seconds);


function str = myFormulas(str, jm)
global ar

varlist = symvar(str)';
svarlist = sym(varlist);
shortlist = {};
for j=1:length(varlist)
    shortlist{j} = sprintf('vj%ijv', j);
end
sshortlist = sym(shortlist);

str = replaceFunctions(str, ar.config.specialFunc);
strsym = arSym(str);
sstrsym = mysubs(strsym, svarlist, sshortlist);

str = latex(sstrsym);

for j=1:length(shortlist)
    str = strrep(str, shortlist{j}, varlist{j});
end

for jx = 1:length(ar.model(jm).x)
    str = strrep(str, sprintf('\\mathrm{%s}', ar.model(jm).x{jx}), ...
        sprintf('\\mathrm{[%s]}', ar.model(jm).x{jx}));
end
for ju = 1:length(ar.model(jm).u)
    str = strrep(str, sprintf('\\mathrm{%s}', ar.model(jm).u{ju}), ...
        sprintf('\\mathrm{[%s]}', ar.model(jm).u{ju}));
end
for jz = 1:length(ar.model(jm).z)
    str = strrep(str, sprintf('\\mathrm{%s}', ar.model(jm).z{jz}), ...
        sprintf('\\mathrm{[%s]}', ar.model(jm).z{jz}));
end
for jc = 1:length(ar.model(jm).pc)
    str = strrep(str, sprintf('\\mathrm{%s}', ar.model(jm).pc{jc}), ...
        sprintf('\\mathrm{V}\\raisebox{-.4ex}{\\tiny %s}', ar.model(jm).c{jc}));
end

str = strrep(str, '_', '\_');
str = strrep(str, '\,', ' \cdot ');

function openDataTable( fid, jm, jd, L, headtab, headstr, unitstr )   
	lp(fid, '\t\\begin{table}\n\t\\dobegincenter');
	if ( headtab > 5 )
        box = sprintf('datamod%s',latexIdentifier(jd+1000*jm+10000*L)) ;
        startFlexbox(fid,  box );
    end                        
	lp(fid, '\t{\\footnotesize');
	lp(fid, '\t\t\\begin{tabular}{%s}', headtab);
	lp(fid, '\t\t\t\\toprule');
	lp(fid, '\t\t\t %s \\\\', headstr);
	lp(fid, '\t\t\t %s \\\\', unitstr);
	lp(fid, '\t\t\t\\midrule');

function closeDataTable( fid, jplot, jm, jd, L, headtab, NAs )
    global ar;

    lp(fid, '\t\t\t\\botrule');
    lp(fid, '\t\t\\end{tabular}}');
    if ( headtab > 5 )
        box = sprintf('datamod%s',latexIdentifier(jd+1000*jm+10000*L)) ;
        endFlexbox( fid,  box );
    end
    if NAs % table contains NaNs:
        lp(fid, '\t\t\\mycaption{Experimental data for the experiment %s. NA indicates data points that are not available, i.e. not measured.}{%s_data}{}', arNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name);
    else
        lp(fid, '\t\t\\mycaption{Experimental data for the experiment %s. }{%s_data}{}', arNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name);
    end
	lp(fid, '\t\\doendcenter');
	lp(fid, '\t\\end{table}');



function str = replaceFunctions(str, funcTypes, checkValidity)

    if (nargin < 3)
        checkValidity = 0;
    end

    str  = char(str);
    replaced = 0;
    for a = 1 : length( funcTypes )
        funcs = findFunc( str, funcTypes{a}{1} );
        argLayout = funcTypes{a}{3};
        
        for b = 1 : length( funcs )
            if ( length( funcs(b).args ) ~= max(argLayout) )
                msg = { 'Invalid number of function argument for function "', ...
                        funcTypes{a}{1}, '" expected ', num2str(max(argLayout)), ...
                        ' got ', num2str( length( funcs(b).args ) ) };
                error( '%s', msg{:} );
            else
                % Determine what the function should be replaced with;
                % feed the appropriate function arguments and replace it
                % Also making sure to use extra brackets for safety (e.g.
                % 5*(a+b) != 5*a+b)
                try
                    to = sprintf( ['(' funcTypes{a}{2} ')'], funcs(b).args{funcTypes{a}{3}} );
                    str = strrep( str, funcs(b).func, to );
                    replaced = replaced + 1;
                catch
                    msg = { 'Failed to replace function ', funcTypes{a}{1}, ...
                        ' in:', funcs(b).func, 'Please expression check for error.' };
                    error( '%s\n', msg{:} );
                end
            end
        end
    end
    
    % Determine whether we got a valid symbolic expression and optionally
    % simplify it
    try
        if (checkValidity)
            str = char( sym( str ) );
        end
        % Enable for input function debug purposes
        % if ( replaced > 0 )
        %     disp(sprintf( '%s =>\n\t\t%s', stro, str ));
        % end
    catch
        msg = { 'Failed to obtain valid expression from: ', ...
                str, 'Please expression check for error.' };
        error( '%s\n', msg{:} );
    end

% Function to scan for specific function name and extract its arguments
function [f] = findFunc( st, funcName )
    loc     = strfind( st, [funcName '('] );
    if ( length(loc) > 0 ) %#ok
        for a = 1 : length( loc )
            f(a) = fetchArgs( st(loc(a):end) );
            f(a).fin = f(a).fin + loc(a)-1;
        end
    else
        f = [];
    end    
    
% Function to fetch function arguments
function f = fetchArgs( st )
    commas  	= [];
    cur         = 0;
    brackets    = 0;
    while( brackets == 0 )
        cur = cur + 1;
        if ( cur > length( st ) )
            error( 'Malformed input string for argument fetcher: \n%s', st );
        end
        if ( brackets < 0 )
            error( 'Malformed input string for argument fetcher: \n%s', st );
        end
        if ( st( cur ) == '(' )
            brackets = brackets + 1;
        end
        if ( st( cur ) == ')' )
            brackets = brackets - 1;
        end        
    end
    if ( brackets < 0 )
        error( 'Malformed input string for argument fetcher: \n%s', st );
    end    
    
    f.name = strtrim( st(1:cur-1) );
    stPos = cur;
    
    while( brackets > 0 )
        cur = cur + 1;
        if ( cur > length( st ) )
            error( 'Malformed input string for argument fetcher: \n%s', st );
        end            
        if ( st( cur ) == '(' )
            brackets = brackets + 1;
        end
        if ( st( cur ) == ')' )
            brackets = brackets - 1;
        end
        if ( ( st( cur ) == ',' ) && ( brackets == 1 ) )
            commas(end+1) = cur;
        end
    end
    
    f.fin    = cur;
        
    list = [stPos, commas, f.fin];
    for b = 1 : length( list ) - 1
        f.args{b} = strtrim( st(list(b)+1:list(b+1)-1) );
    end
    
    f.func = st(1:cur);
    
    
    
function fprintnumtab(fid, num)
if isnan( num )
    fprintf(fid, '& NA ');
else
    fprintf(fid, '& %s ', sprintf('%g', num));
end

% better subs
function out = mysubs(in, old, new)
global ar
if(~isnumeric(in) && ~isempty(old) && ~isempty(char(symvar(in))))
    if(ar.config.matlabVersion>=8.1)
        out = subs(in, old(:), new(:));
    else
        out = subs(in, old(:), new(:), 0);
    end
else
    out = in;
end

function mod = getModifierStr(jm,jv,useNeg,str,sources,targets) %#ok

global ar

field1 = ['qdvd' str '_negative'];
field2 = ['qdvd' str '_nonzero'];

if(useNeg)
    qneg = ar.model(jm).(field1)(jv,:);
else
    qneg = ~ar.model(jm).(field1)(jv,:);
end

isfirst = true;
mod = '';
for jx = find(qneg & ar.model(jm).(field2)(jv,:))
    if(sum(ismember(sources, strrep(ar.model(jm).(str){jx}, '_', '\_'))) ...
            + sum(ismember(targets, strrep(ar.model(jm).(str){jx}, '_', '\_'))) == 0)
        if(~isfirst)
            mod = [mod ', '];
        end
        mod = [mod ...
            sprintf('%s', strrep(ar.model(jm).(str){jx}, '_', '\_'))];
        if(isfirst)
            isfirst = false;
        end
    end
end

% Functions for making flexible tables that scale down when they are bigger
% than the page
function startFlexbox( fid, name )
    lp(fid, '\\newsavebox{\\%s}\\savebox{\\%s}{', name, name);


function endFlexbox( fid, name )
    lp(fid, '}\\newlength{\\q%sw} \\settowidth{\\q%sw}{\\usebox\\%s}\\ifnum\\q%sw>\\linewidth\\setlength{\\q%sw}{\\linewidth}\\fi\\resizebox{\\q%sw}{!}{\\usebox\\%s}', name, name, name, name, name, name, name);


function hash = hashString(str)
    algs = {'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
    checksum = java.security.MessageDigest.getInstance(algs{2});

    checksum.update(uint8(str(:)));
    h = typecast(checksum.digest,'uint8');
    hash = dec2hex(h)';
    hash = ['F', hash(:)'];
    clear checksum
    
function formula = cleanFormula(str) %#ok

    formula = strrep(str, '\mathrm{', '');
    formula = strrep(formula, '\', '');
    formula = strrep(formula, '}', '' );
    formula = strrep(formula, '{', '' );
    
% Latex doesn't like numbers in identifiers, so we have this
function str = latexIdentifier( j )
    str = char(num2str(j)+97-48);
    
% Find condition specific parameter transforms and remove the observation parameters
% which are set to zero; and transforms which are the same as in the main
% model. fpSymString contains a list of the model variables transformed to
% syms and then transformed back to strings again, in order to be able to
% compare with data specific trafo's that underwent the same procedure
function [condTrans, names] = conditionSpecificParameters( jplot, jm, fpSymString, opts )

    global ar;
    
    % Grab all the condition specific conditions
    condTrans = struct();
    names = struct();
    numDataLinks = length(ar.model(jm).plot(jplot).dLink);
    for jdl = 1 : numDataLinks
        jd2 = ar.model(jm).plot(jplot).dLink(jdl);

        % Setup list of observational
        % parameters (could be optimized)
        obsParameters = strrep(ar.model(jm).data(jd2).py, 'filename', stripRandom( jm, jd2 ) );
        obsParameters = [ obsParameters ; strrep(ar.model(jm).data(jd2).pystd, 'filename', stripRandom( jm, jd2 ) ) ];                              

        for jpp = 1 : length( ar.model(jm).data(jd2).pold )
            skip = 0;
            parameterName = ar.model(jm).data(jd2).pold{jpp};

            % Hash it, to avoid going over the
            % struct field name length
            hash            = hashString(parameterName);

            % Check whether it doesn't transform
            % anything at all! Note that RANDOMS explicitly ignored
            tFrom = ar.model(jm).data(jd2).pold{jpp};
            if ( ~(opts.keeprandoms) )
                for jra = 1 : length( ar.model(jm).data(jd2).prand )
                    tFrom = strrep( tFrom, ar.model(jm).data(jd2).prand{jra}, sprintf('%s%s', ar.model(jm).data(jd2).prand{jra}, ar.model(jm).data(jd2).fprand(jra) ) );
                end
            end
            if strcmp( tFrom, ar.model(jm).data(jd2).fp{jpp} )
                skip = 1;
            end

            % Check whether it is an observational parameter
            if (~skip) && ( max( strcmp( parameterName, obsParameters ) ) )
                % If so, check whether it is zero
                % and if so => skip it
                if strcmp( ar.model(jm).data(jd2).fp{jpp}, '0' )
                    skip = 1;
                end
                
                % Check whether it exists in the base model,
                % and whether it is identical. If so => skip it
                if ( ~skip )
                    replacementID = find(strcmp( ar.model(jm).p, strrep( parameterName, ar.model(jm).data(jd2).name, 'filename' ) ) );
                    if ( ~isempty(replacementID) )
                        % It does exist! Is it the same? Look at
                        % the precached model sym strings
                        mainModelReplacement = fpSymString{replacementID};

                        expression = char( arSym( strrep( ar.model(jm).data(jd2).fp{jpp}, ar.model(jm).data(jd2).name, 'filename' ) ) );
                        if ( strcmp( mainModelReplacement, expression ) )
                            skip = 1;
                        end
                    end                            
                end
            end

            if ( ~skip )
                if ( ~isfield( condTrans, hash ) )
                    condTrans.(hash) = cell( numDataLinks, 1 );
                end
                condTrans.(hash){jdl}   = ar.model(jm).data(jd2).fp{jpp};
                names.(hash)            = parameterName;                            
            end
        end
    end

    %% make empty cells a cell of empty strings 
    vars = fieldnames(condTrans);
    for jv=1:length(vars)
        for j=1:length(condTrans.(vars{jv}))
            if isempty(condTrans.(vars{jv}){j})
                condTrans.(vars{jv}){j} = '';
            end
        end
    end
    
    % Remove the transforms that are the same as in
    % the main model
    for jv = 1 : length( vars )
        localReplacement = unique(condTrans.(vars{jv}));
        if ( length(localReplacement) == 1 )
            replacementID = find(strcmp( ar.model(jm).p, names.(vars{jv}) ));
            if ( ~isempty(replacementID) )
                mainModelReplacement = fpSymString{replacementID};
                if ( strcmp( mainModelReplacement, char(arSym(localReplacement)) ) )
                    condTrans = rmfield( condTrans, vars{jv} );
                end
            end
        end
    end
    
% Performs regexprep which transforms func(t,args) => func(t)
function str = repFunc( str, funcName )

    % Pattern that matches func(args)
    pattern = sprintf('%s[\\(](\\w)([\\[\\]\\-\\.\\s,\\w]*)[\\)]', funcName); %\w*
    
    % Compute the mask for the printf
    % Performs regexprep which transforms func(t,args) => func(t)
    mask = regexprep(str, pattern, sprintf('%s($1)', funcName) );
    
    str = sprintf( mask );
    
function str = stripRandom( jm, jd2 )
    global ar;
    str = ar.model(jm).data(jd2).name;
    if ( isfield( ar.model(jm).data(jd2), 'fprand' ) && isfield( ar.model(jm).data(jd2), 'prand' ) )
        for a = 1 : length( ar.model(jm).data(jd2).prand )
            str = strrep(str, [ '_', ar.model(jm).data(jd2).prand{a}, ar.model(jm).data(jd2).fprand(a) ], '' );
        end
    end
    
% Function which makes a hash of all the file identifiers
function str = PTI(jm,  str, pti )
    global ar;
    
    if ( pti == 1 )
        %fileList    = {ar.model.plot.name};
        fileList = cell(0);
        for jm=1:length(ar.model)
            fileList    = unique(union(fileList,unique({ar.model(jm).data.name})));
        end
        for a = 1 : length( fileList )
            str = strrep( str, fileList{a}, num2str(a) );
        end
    end
    
    
function str = arNameTrafo(str)
    str = strrep(str, '_', '\_\-');
    str = strrep(str, '%', '\%');
    
function str = alternate(i)
    if mod(i,2) == 1
        str = '\altrowcol ';
    else
        str = '';
    end
    
    
function [opts] = targSwitch( switches, extraArgs, description, varargin )

    for a = 1 : length(switches)
        if ( extraArgs(a) == 0 )
            opts.(lower(switches{a})) = 0;
        else
            opts.(lower(switches{a})) = {};
        end
    end
    
    a = 1;
    if ~isempty( varargin{1} )
        while( a <= length( varargin{1} ) )
            if ( max( strcmpi( switches, varargin{1}{a} ) ) == 0 )
                error( 'Invalid switch argument was provided %s', varargin{1}{a} );
            end
            if ( extraArgs( strcmpi( switches, varargin{1}{a} ) ) == 0 )
                opts.(lower(varargin{1}{a})) = 1;
            else
                try
                    opts.(lower(varargin{1}{a})) = varargin{1}{a+1};
                    a = a + 1;
                catch
                    error( 'Did not provide arguments for flag %s', varargin{1}{a} );
                end
            end
            a = a + 1;
        end
    end
    for a = 1 : length( switches )
        if ( extraArgs(a) == 0 )
            fprintf( '%s\n', description{a}{2-opts.(lower(switches{a}))} );
        else
            fprintf( '%s\n', description{a}{1+isempty(opts.(lower(switches{a})))} );
        end
    end
    
    
function setchi2fields(jm, jplot)
    global ar;
    
    % chi^2, ndata and dr_times
    chi2 = zeros(1,ar.model(jm).plot(jplot).ny);
    ndata = zeros(1,ar.model(jm).plot(jplot).ny);
    %dr_times = [];
    for jd = ar.model(jm).plot(jplot).dLink
        if(isfield(ar.model(jm),'data'))              
            ny = length(ar.model(jm).data(jd).y);
            for jy = 1:ny
                % chi^2 & ndata
                if(ar.model(jm).data(jd).qFit(jy)==1)
                    chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2(jy);
                    ndata(jy) = ndata(jy) + ar.model(jm).data(jd).ndata(jy);
                    if(ar.config.fiterrors==1)
                        chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2err(jy);
                    end
                end
            end
        end
    end
        
    ar.model(jm).plot(jplot).chi2  = sum(chi2);
    ar.model(jm).plot(jplot).ndata = sum(ndata);
    
function qF = isFitted(jm, jp)
    global ar;

    qF = 0;
    for jd = 1 : length( ar.model(jm).plot(jp).dLink )
        qF = max( union( qF, ar.model(jm).data(ar.model(jm).plot(jp).dLink(jd)).qFit ) );
    end
    
function match = masktest( str, masks )
    match = 0;
    if ( iscell(masks) )
        for a = 1 : length( masks )
            match = max( [ match, regexp( str, masks{a} ) > 0 ] );
        end
    else
        match = max( [ match, regexp( str, masks ) > 0 ] );
    end
    
function wrapped = funcWrap( str, maxLen, newl )
    c = 1;
    cChunk = 0;
    operators = '+-*/';
    cb = 0;
    rb = 0;
    
    wrapped = '';
    while ( c <= length( str ) )
        if ( str(c) == '(' );
            rb = rb + 1;
        end
        if ( str(c) == ')' )
            rb = rb - 1;
        end
        if ( str(c) == '{' );
            cb = cb + 1;
        end
        if ( str(c) == '}' )
            cb = cb - 1;
        end      
        if ( ismember( str(c), operators ) )
            % Potential split
            if ( ( cb == 0 ) && ( rb == 0 ) )
                if ( cChunk > maxLen )
                    % Next line!
                    wrapped = sprintf('%s%s', wrapped, newl);
                    cChunk = 0;
                end
            end
        end
        wrapped = [wrapped str(c)];
        cChunk = cChunk + 1;
        c = c + 1;
    end
    
function volP = getVolume(jm, jc)
    global ar;
    
    volP = ar.model(jm).pc{jc};
    if ( isempty( str2num( volP ) ) ) %#ok
        volP = num2str( str2num( ar.model(jm).fp{ find( strcmp( ar.model(jm).p, volP ) ) } ) ); %#ok
	end