% create report

function arReport(project_name)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

savePath = [arSave '/Latex'];

if(~exist([cd '/' savePath], 'dir'))
    mkdir([cd '/' savePath])
end

% empty LD_LIBRARY_PATH (MATLAB-shipped libraries conflict with libs pdflatex is linked to)
library_path = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH', '');

% copy lib.bib
copyfile(which('lib.bib'), [savePath '/lib.bib']);

% latex packages
copyfile(which('assurechemist.sty'), [savePath '/assurechemist.sty']);
copyfile(which('chemist.sty'), [savePath '/chemist.sty']);

% latex file
fname = 'report.tex';
fnamebib = 'report.aux';
fprintf('writing latex, file %s...', fname);
fid = fopen([savePath '/' fname], 'w');

%% Head
lp(fid, '\\nonstopmode');
lp(fid, '\\documentclass[10pt, oneside, fleqn]{article}');

lp(fid, '\\usepackage[utf8]{inputenc}');
lp(fid, '\\usepackage{amsmath, amsthm, amssymb}');
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
lp(fid, '\\usepackage[normal,footnotesize,bf]{caption}');
lp(fid, '\\usepackage{subfig}');
lp(fid, '\\usepackage{sidecap}');
lp(fid, '\\captionsetup[subfloat]{position=top}');

lp(fid, '\\setlength{\\textheight}{22 cm}');
lp(fid, '\\setlength{\\textwidth}{16 cm}');
lp(fid, '\\setlength{\\topmargin}{-1.5 cm}');
lp(fid, '\\setlength{\\hoffset}{-2 cm}');

lp(fid, '\\usepackage[colorlinks=true, linkcolor=blue, citecolor=blue, filecolor=blue, urlcolor=blue]{hyperref}\n');

lp(fid, '\n\\definecolor{mygray}{gray}{.5}');
lp(fid, '\n\\definecolor{red}{rgb}{.8,0,0}');
farben = lines(7);
for j=1:length(farben(:,1))
    lp(fid, '\\definecolor{line%i}{rgb}{%f,%f,%f}', j, farben(j,1), farben(j,2), farben(j,3));
end

lp(fid, '\\newcommand{\\toprule}{\\hline\\hline}');
lp(fid, '\\newcommand{\\midrule}{\\hline}');
lp(fid, '\\newcommand{\\botrule}{\\hline\\hline}');
lp(fid, '\\newcommand{\\dobegincenter}{\\begin{center}}');
lp(fid, '\\newcommand{\\doendcenter}{\\end{center}}');

lp(fid, '\\newcommand{\\mycaption}[3]{\\caption{\\textbf{#1}\\\\ #3 \\label{#2}}}');

lp(fid, '\n\\begin{document}\n');

%% Header
% lp(fid, ['\\title{Data-2-Dynamics Software' ...
%     '\\footnote{Website: \\href{https://bitbucket.org/d2d-development/d2d-software}' ...
%     '{\\url{https://bitbucket.org/d2d-development/d2d-software}} \\\\ ', ...
%     'Reference: \\citet{Raue:2012zt}}'...
%     ' \\\\ Modeling Report}']);

if(nargin==0)
    project_name = 'Data 2 Dynamics Software -- Modeling Report';
end

lp(fid, '\\title{%s}', project_name);
lp(fid, '\\author{%s}', ar.config.username);
lp(fid, '\\date{%s}', datestr(now));
lp(fid, '\\maketitle\n');

if(nargin>0)
    lp(fid, '\\noindent {\\bf Data 2 Dynamics Software -- Modeling Report}\\\\\\\\');
end

lp(fid, ['\\noindent {\\bf Website:} \\href{https://bitbucket.org/d2d-development/d2d-software}' ...
    '{\\url{https://bitbucket.org/d2d-development/d2d-software}} \\\\\\\\']);
lp(fid, ['{\\bf Reference:} A~Raue, M~Schilling, J~Bachmann, A~Matteson, M~Schelker, ' ...
    'D~Kaschek, S~Hug, C~Kreutz, BD~Harms, F~Theis, U~Klingm{\\"u}ller, and J~Timmer.', ...
    ' Lessons learned from quantitative dynamical modeling in systems biology.', ...
    ' {\\em PLOS ONE}, 8(9):e74335, 2013.']);

lp(fid, '\\tableofcontents');

N = 10;

for jm=1:length(ar.model)
    lp(fid, '\\clearpage\n');
    lp(fid, '\\section{Model: %s}\n', myNameTrafo(ar.model(jm).name));
    
    %% descriptions
    if(~isempty(ar.model(jm).description))
        lp(fid, '\\subsection{Comments}');
        for jd=1:length(ar.model(jm).description)
            lp(fid, '%s\\\\', strrep(strrep(ar.model(jm).description{jd}, '%', '\%'), '_', '\_'));
        end
    end
    
    %% species
    if(~isempty(ar.model(jm).x))
        lp(fid, '\\subsection{Dynamic variables}');
        lp(fid, 'The model contains %i dynamic variables:', length(ar.model(jm).x));
        lp(fid, '\\begin{itemize}');
        for jx = 1:length(ar.model(jm).x)
            lp(fid, '\\item {\\bf Dynamic variable %i:} %s', ...
                jx, strrep(ar.model(jm).x{jx}, '_', '\_'));
            
            lp(fid, '{\\footnotesize');
            lp(fid, '\\begin{equation}');
            lp(fid, '\\begin{aligned}');
            lp(fid, '%s(%s=0) & = & %s \\label{%s}', ...
                myFormulas(ar.model(jm).x{jx}, jm), ...
                myFormulas(ar.model(jm).t, jm), ...
                myFormulas(ar.model(jm).px0{jx}, jm), ...
                sprintf('%s_init%i', ar.model(jm).name, jx));
            lp(fid, '\\end{aligned}');
            lp(fid, '\\end{equation}}');
            
            lp(fid, 'Unit: %s [%s]', ar.model(jm).xUnits{jx,3}, ar.model(jm).xUnits{jx,2});
        end
        lp(fid, '\\end{itemize}');
    end
    
    %% inputs
    if(~isempty(ar.model(jm).u))
        lp(fid, '\\subsection{Input variables}');
        lp(fid, 'The model contains %i external inputs variables:', length(ar.model(jm).u));
        lp(fid, '\\begin{itemize}');
        for ju = 1:length(ar.model(jm).u)
            lp(fid, '\\item {\\bf Input variable %i:} %s', ...
                ju, strrep(ar.model(jm).u{ju}, '_', '\_'));
            
            if(~isempty(ar.model(jm).fu{ju}))
                lp(fid, '{\\footnotesize');
                lp(fid, '\\begin{equation}');
                lp(fid, '\\begin{aligned}');
                lp(fid, '%s(%s) & = & %s \\label{%s}', ...
                    myFormulas(ar.model(jm).u{ju}, jm), ...
                    myFormulas(ar.model(jm).t, jm), ...
                    myFormulas(ar.model(jm).fu{ju}, jm), ...
                    sprintf('%s_input%i', ar.model(jm).name, ju));
                lp(fid, '\\end{aligned}');
                lp(fid, '\\end{equation}}');
            end
            
            lp(fid, 'Unit: %s [%s]', ar.model(jm).uUnits{ju,3}, ar.model(jm).uUnits{ju,2});
        end
        lp(fid, '\\end{itemize}');
    end
    
    %% reactions
    if(~isempty(ar.model(jm).x))
        lp(fid, '\\subsection{Reactions}');
        lp(fid, 'The model contains %i reactions:', ...
            length(ar.model(jm).fv));
        lp(fid, '\\begin{itemize}');
        for jv = 1:length(ar.model(jm).fv)
            lp(fid, '\\item {\\bf Reaction %i:}', jv);
            
            % source
            isfirst = true;
            sourcestr = '';
            sources = {};
            for jx = 1:length(ar.model(jm).fx)
                if(ar.model(jm).N(jx,jv)<0)
                    sources{end+1} = strrep(ar.model(jm).x{jx}, '_', '\_');
                    if(~isfirst)
                        sourcestr = [sourcestr '+ '];
                    end
                    if(abs(ar.model(jm).N(jx,jv))>1)
                        sourcestr = [sourcestr ...
                            sprintf('%i \\cdot %s ', abs(ar.model(jm).N(jx,jv)), ...
                            strrep(ar.model(jm).x{jx}, '_', '\_'))];
                    else
                        sourcestr = [sourcestr ...
                            sprintf('%s ', strrep(ar.model(jm).x{jx}, '_', '\_'))];
                    end
                    if(isfirst)
                        isfirst = false;
                    end
                end
            end
            if(isfirst)
                sourcestr = [sourcestr '\varnothing'];
            end
                        
            % target
            isfirst = true;
            targetstr = '';
            targets = {};
            for jx = 1:length(ar.model(jm).fx)
                if(ar.model(jm).N(jx,jv)>0)
                    targets{end+1} = strrep(ar.model(jm).x{jx}, '_', '\_');
                    if(~isfirst)
                        targetstr = [targetstr '+ '];
                    end
                    if(abs(ar.model(jm).N(jx,jv))>1)
                        targetstr = [targetstr ...
                            sprintf('%i \\cdot %s ', abs(ar.model(jm).N(jx,jv)), strrep(ar.model(jm).x{jx}, '_', '\_'))];
                    else
                        targetstr = [targetstr ...
                            sprintf('%s ', strrep(ar.model(jm).x{jx}, '_', '\_'))];
                    end
                    if(isfirst)
                        isfirst = false;
                    end
                end
            end
            if(isfirst)
                targetstr = [targetstr '\varnothing'];
            end
            lp(fid, '\\\\');
            
            % modifiers
            mod_pos = getModifierStr(jm,jv,false,'x',sources,targets);
            mod_posu = getModifierStr(jm,jv,false,'u',sources,targets);
            if(~isempty(mod_posu))
                if(~isempty(mod_pos))
                    mod_pos = [mod_pos ', ' mod_posu];
                else
                    mod_pos = mod_posu;
                end
            end
            mod_neg = getModifierStr(jm,jv,true,'x',sources,targets);
            mod_negu = getModifierStr(jm,jv,true,'u',sources,targets);
            if(~isempty(mod_negu))
                if(~isempty(mod_neg))
                    mod_neg = [mod_neg ', ' mod_negu];
                else
                    mod_neg = mod_negu;
                end
            end
            
            % chem reaction
            lp(fid, '\\begin{chemmath}');
            lp(fid, '%s', sourcestr);
            if(length(mod_pos)>length(mod_neg))
                lp(fid, '\\reactrarrow{1pt}{\\widthof{%s}+0.5cm}{%s}{\\color{red}%s}', ...
                    mod_pos, mod_pos, mod_neg);
            else
                lp(fid, '\\reactrarrow{1pt}{\\widthof{%s}+0.5cm}{%s}{\\color{red}%s}', ...
                    mod_neg, mod_pos, mod_neg);
            end
            lp(fid, '%s', targetstr);
            lp(fid, '\\end{chemmath}');
            
            % rate equation
            lp(fid, '{\\footnotesize');
            lp(fid, '\\begin{equation}');
            lp(fid, '\\begin{aligned}');
            lp(fid, 'v_{%i} & = & %s \\label{%s}', jv, myFormulas(ar.model(jm).fv{jv}, jm), ...
                sprintf('%s_flux%i', ar.model(jm).name, jv));
            lp(fid, '\\end{aligned}');
            lp(fid, '\\end{equation}}');
            
            fprintf(fid, '\n');
        end
        lp(fid, '\\end{itemize}');
    end
    
    %% Structure
    if(isfield(ar.model(jm), 'savePath_Graph') && ~isempty(ar.model(jm).savePath_Graph))
        savePath_Graph = [arSave '/Figures/Network' sprintf('/%s.pdf', ar.model(jm).name)];
        if(exist(savePath_Graph,'file'))
            lp(fid, '\\subsection{Model structure}');
            lp(fid, 'The model structure is depicted in Figure \\ref{%s}.', [ar.model(jm).name '_graph']);
            
            copyfile(savePath_Graph, [savePath '/' ar.model(jm).name '_graph.pdf'])
            captiontext = sprintf('\\textbf{%s network representation.} ', myNameTrafo(ar.model(jm).name));
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
    
    %% ODE system
    if(~isempty(ar.model(jm).x))
        lp(fid, '\\subsection{ODE system}');
        
        lp(fid, '\\noindent The ODE system determining the time evolution of the dynamical variables is given by:');
        lp(fid, '{\\footnotesize');
        lp(fid, '\\begin{eqnarray}');
        for jx=1:size(ar.model(jm).N, 1) % for every species jx
            strtmp = '';
            if(~isempty(ar.model(jm).c))
                qinfluxwitheducts = ar.model(jm).N(jx,:) > 0 & sum(ar.model(jm).N < 0,1) > 0;
                eductcompartment = zeros(size(qinfluxwitheducts));
                for jj=find(qinfluxwitheducts)
                    eductcompartment(jj) = unique(ar.model(jm).cLink(ar.model(jm).N(:,jj)<0)); %R2013a compatible
                end
            end
            for jv = find(ar.model(jm).N(jx,:))
                if(abs(ar.model(jm).N(jx,jv))~=1)
                    strtmp = [strtmp sprintf(' %+i \\cdot v_{%i}', ar.model(jm).N(jx,jv), jv)]; %#ok<*AGROW>
                elseif(ar.model(jm).N(jx,jv)==1)
                    strtmp = [strtmp sprintf(' + v_{%i}', jv)];
                elseif(ar.model(jm).N(jx,jv)==-1)
                    strtmp = [strtmp sprintf(' - v_{%i}', jv)];
                end
                if(~isempty(ar.model(jm).c) && qinfluxwitheducts(jv) && eductcompartment(jv)~=ar.model(jm).cLink(jx))
                    strtmp = [strtmp sprintf(' \\cdot \\frac{%s}{%s}', myFormulas(ar.model(jm).pc{eductcompartment(jv)}, jm), ...
                        myFormulas(ar.model(jm).pc{ar.model(jm).cLink(jx)}, jm))];
                end
                
            end
            if(jx==size(ar.model(jm).N, 1) || mod(jx,N)==0)
                lp(fid, '\t\\mathrm{d}%s/\\mathrm{dt} & = & %s \\label{%s}', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_ode%i', ar.model(jm).name, jx));
            else
                lp(fid, '\t\\mathrm{d}%s/\\mathrm{dt} & = & %s \\label{%s} \\\\', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_ode%i', ar.model(jm).name, jx));
            end
            
            if(mod(jx,N)==0 && jx<size(ar.model(jm).N, 1))
                lp(fid, '\\end{eqnarray}\n');
                lp(fid, '\\begin{eqnarray}');
            end
        end
        lp(fid, '\\end{eqnarray}}\n\n');
        
        lp(fid, '\\noindent The ODE system was solved by a parallelized implementation of the CVODES algorithm \\cite{Hindmarsh:2005fb}.');
        lp(fid, 'It also supplies the parameter sensitivities utilized for parameter estimation.\\\\\n\n');
    end
    
    %% derived
    if(~isempty(ar.model(jm).z))
        lp(fid, '\\subsection{Derived variables}');
        lp(fid, 'The model contains %i derived variables:', length(ar.model(jm).z));
        lp(fid, '\\begin{itemize}');
        for jz = 1:length(ar.model(jm).z)
            lp(fid, '\\item {\\bf Derived variable %i:} %s', ...
                jz, strrep(ar.model(jm).z{jz}, '_', '\_'));
            
            lp(fid, '{\\footnotesize');
            lp(fid, '\\begin{equation}');
            lp(fid, '\\begin{aligned}');
            lp(fid, '%s(%s) & = & %s \\label{%s}', ...
                myFormulas(ar.model(jm).z{jz}, jm), ...
                myFormulas(ar.model(jm).t, jm), ...
                myFormulas(ar.model(jm).fz{jz}, jm), ...
                sprintf('%s_derived%i', ar.model(jm).name, jz));
            lp(fid, '\\end{aligned}');
            lp(fid, '\\end{equation}}');
            
            lp(fid, 'Unit: %s [%s]', ar.model(jm).zUnits{jz,3}, ar.model(jm).zUnits{jz,2});
        end
        lp(fid, '\\end{itemize}');
    end
    
    %% standard observations and error model
    if(isfield(ar.model(jm), 'y'))
        lp(fid, '\\subsection{Observables}');
        lp(fid, 'The model contains %i standard observables:', length(ar.model(jm).y));
        lp(fid, '\\begin{itemize}');
        for jy = 1:length(ar.model(jm).y)
            lp(fid, '\\item {\\bf Observable %i:} %s', ...
                jy, strrep(ar.model(jm).y{jy}, '_', '\_'));
            
            strtmp = myFormulas(ar.model(jm).fy{jy}, jm);
            if(ar.model(jm).logfitting(jy))
                strtmp = ['\log_{10}(' strtmp ')'];
            end
            
            lp(fid, '{\\footnotesize');
            lp(fid, '\\begin{equation}');
            lp(fid, '\\begin{aligned}');
            lp(fid, '%s(%s) & = & %s \\label{%s}', ...
                myFormulas(ar.model(jm).y{jy}, jm), ...
                myFormulas(ar.model(jm).t, jm), ...
                strtmp, ...
                sprintf('%s_std_obs%i', ar.model(jm).name, jy));
            lp(fid, '\\end{aligned}');
            lp(fid, '\\end{equation}}');
            
            lp(fid, 'Unit: %s [%s]; ', ar.model(jm).yUnits{jy,3}, ar.model(jm).yUnits{jy,2});
            
            lp(fid, 'With error model:');
            
            lp(fid, '{\\footnotesize');
            lp(fid, '\\begin{equation}');
            lp(fid, '\\begin{aligned}');
            lp(fid, '\\sigma\\{%s\\}(%s) & = & %s \\label{%s}', ...
                myFormulas(ar.model(jm).y{jy}, jm), ...
                myFormulas(ar.model(jm).t, jm), ...
                myFormulas(ar.model(jm).fystd{jy}, jm), ...
                sprintf('%s_std_err%i', ar.model(jm).name, jy));
            lp(fid, '\\end{aligned}');
            lp(fid, '\\end{equation}}');
        end
        lp(fid, '\\end{itemize}');
    end
    
    %% Conditions
    ccount = 1;
    for jp=1:length(ar.model(jm).fp)
        if(~strcmp(ar.model(jm).p{jp}, ar.model(jm).fp{jp}))
            if(ccount==1)
                lp(fid, '\\subsection{Parameter transformations}');
                lp(fid, '\\noindent The ODE system is modified by the following parameter transformations:');
                lp(fid, '{\\footnotesize');
                lp(fid, '\\begin{eqnarray}');
            end
            
            if(ccount==length(ar.model(jm).fp) || mod(ccount,N)==0)
                lp(fid, '\t%s & \\rightarrow & %s', myFormulas(ar.model(jm).p{jp}, jm), ...
                    myFormulas(ar.model(jm).fp{jp}, jm));
            else
                lp(fid, '\t%s & \\rightarrow & %s \\\\', myFormulas(ar.model(jm).p{jp}, jm), ...
                    myFormulas(ar.model(jm).fp{jp}, jm));
            end
            if(mod(ccount,N)==0 && ccount<length(ar.model(jm).fp))
                lp(fid, '\\end{eqnarray}\n');
                lp(fid, '\\begin{eqnarray}');
            end
            ccount = ccount + 1;
        end
        if(ccount>1 && jp==length(ar.model(jm).fp))
            lp(fid, '\\end{eqnarray}}\n\n');
        end
    end
    
    %% Experiments
    lp(fid, '\\clearpage\n');
    
    for jplot=1:length(ar.model(jm).plot)
        jd = ar.model(jm).plot(jplot).dLink(1);
        if(isfield(ar.model(jm), 'data'))
            lp(fid, '\\clearpage\n');
            lp(fid, '\\subsection{Experiment: %s}\n', myNameTrafo(ar.model(jm).plot(jplot).name));
            
            %% descriptions
            if(~isempty(ar.model(jm).data(jd).description))
                lp(fid, '\\subsubsection{Comments}');
                for jdes=1:length(ar.model(jm).data(jd).description)
                    lp(fid, '%s\\\\', strrep(strrep(ar.model(jm).data(jd).description{jdes}, '%', '\%'), '_', '\_'));
                end
            end
            
            %% fit
            lp(fid, '\\subsubsection{Model fit and plots}');
            lp(fid, '\\noindent The agreement of the model observables and the experimental data, given in Table \\ref{%s_data}, ', ar.model(jm).plot(jplot).name);
            
            if(ar.config.fiterrors == 1)
                lp(fid, 'yields a value of the objective function $-2 \\log(L) = %g$ for %i data points in this data set.', ...
                    2*ar.model(jm).plot(jplot).ndata*log(sqrt(2*pi)) + ar.model(jm).plot(jplot).chi2, ar.model(jm).plot(jplot).ndata);
            else
                lp(fid, 'yields a value of the objective function $\\chi^2 = %g$ for %i data points in this data set.', ...
                    ar.model(jm).plot(jplot).chi2, ar.model(jm).plot(jplot).ndata);
            end
            
            %% plots
            if(isfield(ar.model(jm).plot(jplot), 'savePath_FigY') && ~isempty(ar.model(jm).plot(jplot).savePath_FigY))
                copyfile([ar.model(jm).plot(jplot).savePath_FigY '.pdf'], ...
                    [savePath '/' ar.model(jm).plot(jplot).name '_y.pdf']);
                
                lp(fid, 'The model observables and the experimental data is shown in Figure \\ref{%s}.', [ar.model(jm).plot(jplot).name '_y']);
                captiontext = sprintf('\\textbf{%s observables and experimental data for the experiment.} ', myNameTrafo(ar.model(jm).plot(jplot).name));
                captiontext = [captiontext 'The observables are displayed as solid lines. '];
                captiontext = [captiontext 'The error model that describes the measurement noise ' ...
                    'is indicated by shades.'];
                lpfigure(fid, 1, [ar.model(jm).plot(jplot).name '_y.pdf'], captiontext, [ar.model(jm).plot(jplot).name '_y']);
            end
            
            %% experimental data
            headstr = '';
            headtab = '';
            unitstr = '';
            % time
            unitstr = [unitstr sprintf('%s [%s] ', myNameTrafo(ar.model(jm).data(jd).tUnits{3}), ...
                myNameTrafo(ar.model(jm).data(jd).tUnits{2}))];
            headstr = [headstr ' '];
            headtab = [headtab 'r'];
            % conditions
            if(~isempty(ar.model(jm).data(jd).condition))
                for jp = 1:length(ar.model(jm).data(jd).condition)
                    headstr = [headstr sprintf('& %s ', myNameTrafo(ar.model(jm).data(jd).condition(jp).parameter))];
                    unitstr = [unitstr '& '];
                    headtab = [headtab 'r'];
                end
            end
            % y & ystd headers
            for jy=1:length(ar.model(jm).data(jd).y)
                headstr = [headstr sprintf('& %s ', myNameTrafo(ar.model(jm).data(jd).y{jy}))];
                unitstr = [unitstr sprintf('& %s [%s] ', myNameTrafo(ar.model(jm).data(jd).yUnits{jy,3}), ...
                    myNameTrafo(ar.model(jm).data(jd).yUnits{jy,2}))];
                headtab = [headtab 'r'];
                if(ar.config.fiterrors == -1)
                    headstr = [headstr sprintf('& %s\_std ', myNameTrafo(ar.model(jm).data(jd).y{jy}))];
                    unitstr = [unitstr sprintf('& %s [%s] ', myNameTrafo(ar.model(jm).data(jd).yUnits{jy,3}), ...
                        myNameTrafo(ar.model(jm).data(jd).yUnits{jy,2}))];
                    headtab = [headtab 'r'];
                end
            end
            lp(fid, '\t\\begin{table}');
            lp(fid, '\t\\dobegincenter');
            lp(fid, '\t{\\footnotesize');
            lp(fid, '\t\t\\begin{tabular}{%s}', headtab);
            lp(fid, '\t\t\t\\toprule');
            lp(fid, '\t\t\t %s \\\\', headstr);
            lp(fid, '\t\t\t %s \\\\', unitstr);
            lp(fid, '\t\t\t\\midrule');
            
            for jd2 = ar.model(jm).plot(jplot).dLink
                for j=1:length(ar.model(jm).data(jd2).tExp)
                    fprintf(fid, '%s ', sprintf('%g', ar.model(jm).data(jd2).tExp(j)));
                    
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
                end
            end
            
            lp(fid, '\t\t\t\\botrule');
            lp(fid, '\t\t\\end{tabular}}');
            lp(fid, '\t\t\\mycaption{Experimental data for the experiment %s}{%s_data}{}', myNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name);
            lp(fid, '\t\\doendcenter');
            lp(fid, '\t\\end{table}');
            
            %% trajectories
            if(isfield(ar.model(jm).plot(jplot), 'savePath_FigX') && ~isempty(ar.model(jm).plot(jplot).savePath_FigX))
                lp(fid, ['The trajectories of the input, dynamic and derived variables that ' ...
                    'correspond to the experimental conditions in this experiment are shown in Figure \\ref{%s}.'], ...
                    [ar.model(jm).plot(jplot).name '_x']);
                copyfile([ar.model(jm).plot(jplot).savePath_FigX '.pdf'], ...
                    [savePath '/' ar.model(jm).plot(jplot).name '_x.pdf'])
                
                captiontext = sprintf('\\textbf{%s trajectories of the input, dynamic and derived variables.} ', ....
                    myNameTrafo(ar.model(jm).plot(jplot).name));
                captiontext = [captiontext 'The dynamical behaviour is determined by the ODE system, see Equation '];
                captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \ref{' ...
                    sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
                lpfigure(fid, 1, [ar.model(jm).plot(jplot).name '_x.pdf'], captiontext, [ar.model(jm).plot(jplot).name '_x']);
            end
            if(isfield(ar.model(jm).plot(jplot), 'savePath_FigV') && ~isempty(ar.model(jm).plot(jplot).savePath_FigV))
                lp(fid, 'The reaction fluxes that correspond to the experimental conditions in this experiment are shown in Figure \\ref{%s}.', ...
                    [ar.model(jm).plot(jplot).name '_v']);
                copyfile([ar.model(jm).plot(jplot).savePath_FigV '.pdf'], ...
                    [savePath '/' ar.model(jm).plot(jplot).name '_v.pdf'])
                
                captiontext = sprintf('\\textbf{%s reaction fluxes.} ', myNameTrafo(ar.model(jm).plot(jplot).name));
                captiontext = [captiontext 'The dynamical behaviour is determined by the ODE system, see Equation '];
                captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) ...
                    '} -- \ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
                lpfigure(fid, 1, [ar.model(jm).plot(jplot).name '_v.pdf'], captiontext, [ar.model(jm).plot(jplot).name '_v']);
            end            
            
            %% inputs
            if(~isempty(ar.model(jm).u))
                qmod = ~strcmp(ar.model(jm).fu, ar.model(jm).data(jd).fu);
                if(sum(qmod)>0)
                    lp(fid, '\\subsubsection{Input variables}');
                    lp(fid, 'The following inputs variables are modified in this data set:');
                    lp(fid, '\\begin{itemize}');
                    for ju = find(qmod)
                        lp(fid, '\\item {\\bf Input variable %i:} %s', ju, strrep(ar.model(jm).u{ju}, '_', '\_'));
                        
                        lp(fid, '{\\footnotesize');
                        lp(fid, '\\begin{equation}');
                        lp(fid, '\\begin{aligned}');
                        lp(fid, '%s(%s) & = & %s \\label{%s}', ...
                            myFormulas(ar.model(jm).u{ju}, jm), ...
                            myFormulas(ar.model(jm).t, jm), ...
                            myFormulas(ar.model(jm).data(jd).fu{ju}, jm), ...
                            sprintf('%s_input%i', ar.model(jm).plot(jplot).name, ju));
                        lp(fid, '\\end{aligned}');
                        lp(fid, '\\end{equation}}');
                    end
                    lp(fid, '\\end{itemize}');
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
                lp(fid, '\\subsubsection{Observables}');
                
                % modified observables
                if(sum(qmod)>0)
                    lp(fid, 'The following observables are modified in this data set:');
                    lp(fid, '\\begin{itemize}');
                    for jy = find(qmod)
                        lp(fid, '\\item {\\bf Observable:} %s', ...
                            strrep(ar.model(jm).data(jd).y{jy}, '_', '\_'));
                        
                        strtmp = myFormulas(ar.model(jm).data(jd).fy{jy}, jm);
                        if(ar.model(jm).logfitting(jy))
                            strtmp = ['\log_{10}(' strtmp ')'];
                        end
                        
                        lp(fid, '{\\footnotesize');
                        lp(fid, '\\begin{equation}');
                        lp(fid, '\\begin{aligned}');
                        lp(fid, '%s(%s) & = & %s \\label{%s}', ...
                            myFormulas(ar.model(jm).data(jd).y{jy}, jm), ...
                            myFormulas(ar.model(jm).t, jm), ...
                            strtmp, ...
                            sprintf('%s_obs%i', ar.model(jm).plot(jplot).name, jy));
                        lp(fid, '\\end{aligned}');
                        lp(fid, '\\end{equation}}');
                        
                        lp(fid, 'Unit: %s [%s]; ', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2});
                        
                        lp(fid, 'With error model:');
                        
                        lp(fid, '{\\footnotesize');
                        lp(fid, '\\begin{equation}');
                        lp(fid, '\\begin{aligned}');
                        lp(fid, '\\sigma\\{%s\\}(%s) & = & %s \\label{%s}', ...
                            myFormulas(ar.model(jm).data(jd).y{jy}, jm), ...
                            myFormulas(ar.model(jm).t, jm), ...
                            myFormulas(ar.model(jm).data(jd).fystd{jy}, jm), ...
                            sprintf('%s_err%i', ar.model(jm).plot(jplot).name, jy));
                        lp(fid, '\\end{aligned}');
                        lp(fid, '\\end{equation}}');
                    end
                    lp(fid, '\\end{itemize}');
                end
                
                % additional observables
                if(sum(qadd)>0)
                    lp(fid, 'The following observables are added in this data set:');
                    lp(fid, '\\begin{itemize}');
                    for jy = find(qadd)
                        lp(fid, '\\item {\\bf Observable:} %s', ...
                            strrep(ar.model(jm).data(jd).y{jy}, '_', '\_'));
                        
                        strtmp = myFormulas(ar.model(jm).data(jd).fy{jy}, jm);
                        if(isfield(ar.model(jm),'logfitting') && ar.model(jm).logfitting(jy))
                            strtmp = ['\log_{10}(' strtmp ')'];
                        end
                        
                        lp(fid, '{\\footnotesize');
                        lp(fid, '\\begin{equation}');
                        lp(fid, '\\begin{aligned}');
                        lp(fid, '%s(%s) & = & %s \\label{%s}', ...
                            myFormulas(ar.model(jm).data(jd).y{jy}, jm), ...
                            myFormulas(ar.model(jm).t, jm), ...
                            strtmp, ...
                            sprintf('%s_obs%i', ar.model(jm).plot(jplot).name, jy));
                        lp(fid, '\\end{aligned}');
                        lp(fid, '\\end{equation}}');
                        
                        lp(fid, 'Unit: %s [%s]; ', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2});
                        
                        lp(fid, 'With error model:');
                        
                        lp(fid, '{\\footnotesize');
                        lp(fid, '\\begin{equation}');
                        lp(fid, '\\begin{aligned}');
                        lp(fid, '\\sigma\\{%s\\}(%s) & = & %s \\label{%s}', ...
                            myFormulas(ar.model(jm).data(jd).y{jy}, jm), ...
                            myFormulas(ar.model(jm).t, jm), ...
                            myFormulas(ar.model(jm).data(jd).fystd{jy}, jm), ...
                            sprintf('%s_err%i', ar.model(jm).plot(jplot).name, jy));
                        lp(fid, '\\end{aligned}');
                        lp(fid, '\\end{equation}}');
                    end
                    lp(fid, '\\end{itemize}');
                end
            end
            
            %% conditions
            ccount = 1;
            for jp=1:length(ar.model(jm).data(jd).fp)
                % check if this observable was removed
                wasRemoved = false;
                if(sum(strcmp(ar.model(jm).data(jd).py, ar.model(jm).data(jd2).pold{jp}))>0 || ...
                        sum(strcmp(ar.model(jm).data(jd).pystd, ar.model(jm).data(jd2).pold{jp}))>0)
                    if(sum(strcmp(ar.p, ar.model(jm).data(jd2).pold{jp}))==0)
                        wasRemoved = true;
                    end
                end
                
                if(~wasRemoved)
                    % check if this is really a condition
                    qlocalcondi = false;
                    arraystr = 'll';
                    for jd2 = ar.model(jm).plot(jplot).dLink
                        qlocalcondi = qlocalcondi || ~strcmp(ar.model(jm).data(jd2).pold{jp}, ar.model(jm).data(jd2).fp{jp});
                        arraystr = [arraystr 'r'];
                        if(jd2 ~= ar.model(jm).plot(jplot).dLink(end))
                            arraystr = [arraystr '|'];
                        end
                    end
                    
                    if(qlocalcondi)
                        % check is already shown in model part
                        qdyn = ismember(ar.model(jm).p, ar.model(jm).data(jd).pold{jp}); %R2013a compatible
                        if(sum(qdyn)>0)
                            qalreadyset = true;
                            for jd2 = ar.model(jm).plot(jplot).dLink
                                qalreadyset = qalreadyset && strcmp(ar.model(jm).fp{qdyn}, ar.model(jm).data(jd2).fp{jp});
                            end
                        else
                            qalreadyset = false;
                        end
                        
                        if(~qalreadyset)
                            if(ccount==1)
                                lp(fid, '\\subsubsection{Conditions}');
                                lp(fid, ['\\noindent To evaluate the ODE system of Equation \\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \\ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}']);
                                lp(fid, 'for the conditions in this experiment, the following parameter transformations are applied:');
                                lp(fid, '{\\footnotesize');
                                lp(fid, '\\begin{displaymath}');
                                lp(fid, '\\begin{array}{%s}', arraystr);
                            else
                                if(~(mod(ccount-1,N)==0 && ccount-1<length(ar.model(jm).data(jd).fp)))
                                    lp(fid, ' \\\\');
                                end
                            end
                            
                            lp(fid, '\t%s & \\rightarrow & ', myFormulas(ar.model(jm).data(jd).pold{jp}, jm));
                            for jd2 = ar.model(jm).plot(jplot).dLink
                                try
                                    ddouble = double(sym(ar.model(jm).data(jd2).fp{jp}));
                                    if(ddouble ~= 0 && abs(log10(ddouble)) > 3)
                                        ddouble = sprintf('%.1e', ddouble);
                                        ddouble = strrep(ddouble, 'e', '\cdot 10^{');
                                        ddouble = [ddouble '}'];
                                    else
                                        ddouble = sprintf('%g', ddouble);
                                    end
                                    lp(fid, '%s', ddouble);
                                catch  %#ok<CTCH>
                                    lp(fid, '%s', myFormulas(ar.model(jm).data(jd2).fp{jp}, jm));
                                end
                                if(length(ar.model(jm).plot(jplot).dLink)>1 && jd2~=ar.model(jm).plot(jplot).dLink(end))
                                    lp(fid, '& ');
                                end
                            end
                            
                            if(mod(ccount,N)==0 && ccount<length(ar.model(jm).data(jd).fp))
                                lp(fid, '\\end{array}');
                                lp(fid, '\\end{displaymath}\n');
                                lp(fid, '\\begin{displaymath}');
                                lp(fid, '\\begin{array}{%s}', arraystr);
                            end
                            ccount = ccount + 1;
                        end
                    end
                end    
                
                if(ccount>1 && jp==length(ar.model(jm).data(jd).fp))
                	lp(fid, '\\end{array}');
                	lp(fid, '\\end{displaymath}\n}\n\n');
                end
            end
        end
    end
end

%% Parameters
lp(fid, '\\clearpage\n');
lp(fid, '\\section{Estimated model parameters} \\label{estimatedparameters}\n');

if(isfield(ar,'ndata') && isfield(ar,'chi2fit') && isfield(ar,'chi2'))
    if(ar.config.fiterrors == 1)
        lp(fid, 'In total %i parameters are estimated from the experimental data, yielding a value of the objective function $-2 \\log(L) = %g$ for a total of %i data points.', ...
            sum(ar.qFit==1), 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit, ar.ndata);
    else
        lp(fid, 'In total %i parameters are estimated from the experimental data, yielding a value of the objective function $\\chi^2 = %g$ for a total of %i data points.', ...
            sum(ar.qFit==1), ar.chi2, ar.ndata);
    end

    lp(fid, 'The model parameters were estimated by maximum likelihood estimation applying the MATLAB lsqnonlin algorithm.');
end

N = 50;
ntables = ceil(length(ar.p)/N);
if(ntables>1)
    lp(fid, 'In Table \\ref{paratable1} -- \\ref{paratable%i} the estimated parameter values are given.', ntables);
else
    lp(fid, 'In Table \\ref{paratable1} the estimated parameter values are given.');
end
lp(fid, 'Parameters highlighted in red color indicate parameter values close to their bounds.');
lp(fid, 'The parameter name prefix init\\_ indicates the initial value of a dynamic variable.');
lp(fid, 'The parameter name prefix offset\\_ indicates a offset of the experimental data.');
lp(fid, 'The parameter name prefix scale\\_ indicates a scaling factor of the experimental data.');
lp(fid, 'The parameter name prefix sd\\_ indicates the magnitude of the measurement noise for a specific measurement.\\\\');

pTrans = ar.p;
pTrans(ar.qLog10==1) = 10.^pTrans(ar.qLog10==1);

lp(fid, '\t\\begin{table}');
lp(fid, '\t\\dobegincenter');
lp(fid, '\t{\\footnotesize');
lp(fid, '\t\t\\begin{tabular}{llllllll}');
lp(fid, '\t\t\t\\toprule');
lp(fid, '\t\t\t & name & $\\theta_{min}$ & $\\hat\\theta$ & $\\theta_{max}$ & log & non-log $\\hat \\theta$ & fitted \\\\');
lp(fid, '\t\t\t\\midrule');
count = 1;
for j=1:length(ar.p)
    if(mod(j,N)==0)
        lp(fid, '\t\t\t\\botrule');
        lp(fid, '\t\t\\end{tabular}}');
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
        lp(fid, '\t{\\footnotesize');
        lp(fid, '\t\t\\begin{tabular}{llllllll}');
        lp(fid, '\t\t\t\\toprule');
        lp(fid, '\t\t\t & name & $\\theta_{min}$ & $\\hat\\theta$ & $\\theta_{max}$ & log & non-log $\\hat \\theta$ & fitted \\\\');
        lp(fid, '\t\t\t\\midrule');
        
        count = count + 1;
    end
    
    if(ar.qFit(j)==1)
        if(ar.p(j) - ar.lb(j) < 0.1 || ar.ub(j) - ar.p(j) < 0.1)
            lp(fid, '\t\t\t\\color{red}{%i} & \\color{red}{%s} & \\color{red}{%+8.0g} & \\color{red}{%+8.4f} & \\color{red}{%+8.0g} & \\color{red}{%i} & \\color{red}{%s} & \\color{red}{%i} \\\\', ...
                j, strrep(ar.pLabel{j},'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
                ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
        else
            lp(fid, '\t\t\t%i & %s & {%+8.0g} & {%+8.4f} & {%+8.0g} & %i & %s & %i \\\\', ...
                j, strrep(ar.pLabel{j},'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
                ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
        end
    else
        lp(fid, '\t\t\t\\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%+8.0g} & \\color{mygray}{%+8.4f} & \\color{mygray}{%+8.0g} & \\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%i} \\\\', ...
            j, strrep(ar.pLabel{j},'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
            ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
    end
end
lp(fid, '\t\t\t\\botrule');
lp(fid, '\t\t\\end{tabular}}');
lp(fid, '\t\t\\mycaption{Estimated parameter values}{paratable%i}', count);
lp(fid, '{$\\hat \\theta$ indicates the estimated value of the parameters.');
lp(fid, '$\\theta_{min}$ and $\\theta_{max}$ indicate the upper and lower bounds for the parameters.');
lp(fid, 'The log-column indicates if the value of a parameter was log-transformed.');
lp(fid, 'If log = 1 the non-log-column indicates the non-logarithmic value of the estimate.');
lp(fid, 'The fitted-column indicates if the parameter value was estimated (1), was temporarily fixed (0) or if its value was fixed to a constant value (2).}');
lp(fid, '\t\\doendcenter');
lp(fid, '\t\\end{table}');

%% PLE
plePath = [arSave '/PLE'];
if(exist(plePath,'dir'))
    lp(fid, '\\clearpage\n');
    lp(fid, '\\section{Profile likelihood of model parameters}\n');
    
    S = load([plePath '/results.mat']);
    
    lp(fid, 'In order to evaluate the identifiability of the model parameters and to assess confidence intervals, ');
    lp(fid, 'the profile likelihood \\cite{Raue:2009ec} was calculated.');
    lp(fid, 'The mean calculation time of the profile likelihood per parameter was %s $\\pm$ %s.', ...
        secToHMS(mean(S.pleGlobals.timing(S.pleGlobals.q_fit(size(S.pleGlobals.timing))))), ...
        secToHMS(std(S.pleGlobals.timing(S.pleGlobals.q_fit(size(S.pleGlobals.timing))))));
    
    % Multiplot
    if(isfield(S.pleGlobals, 'fighandel_multi'))
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
    if(isfield(S.pleGlobals, 'figPath'))
        count = 0;
        for j=1:length(S.pleGlobals.figPath)
            if(~isempty(S.pleGlobals.chi2s{j}))
                count = count + 1;
            end
        end
        
        lp(fid, 'In Figure \\ref{ple1} -- \\ref{ple%i} the profile likelihood of each parameter is shown in more detail.', count);
        lp(fid, 'Also the functional relations to the remaining parameter are displayed.');
        
        count = 0;
        for j=1:length(S.pleGlobals.figPath)
            if(~isempty(S.pleGlobals.chi2s{j}))
                lp(fid, '\\clearpage\n');
                count = count + 1;
                targetstr = [savePath '/' S.pleGlobals.p_labels{j} '.pdf'];
                if(ispc)
                    print('-dpdf', targetstr);
                elseif(ismac)
                    system(['/usr/local/bin/ps2pdf  -dEPSCrop ' S.pleGlobals.figPath{j} '.eps ' targetstr]);
                else
                    system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' S.pleGlobals.figPath{j} '.eps ' targetstr]);
                end
                
                captiontext = sprintf('\\textbf{Profile likelihood of parameter %s}\\\\', strrep(S.pleGlobals.p_labels{j},'_','\_'));
                captiontext = [captiontext 'Upper panel: The solide line indicates the profile likelihood. '];
                captiontext = [captiontext 'The broken line indicates the threshold to assess confidence intervals. '];
                captiontext = [captiontext 'The asterisk indicate the optimal parameter values. '];
                captiontext = [captiontext sprintf('Lower panel: The functional relations to the other parameters along the profile likelihood of %s are displayed. ', strrep(S.pleGlobals.p_labels{j},'_','\_'))];
                captiontext = [captiontext 'In the legend the top five parameters showing the strongest variations are given. '];
                captiontext = [captiontext sprintf('The calculation time was %s.', secToHMS(S.pleGlobals.timing(j)))];
                lpfigure(fid, 1, [S.pleGlobals.p_labels{j} '.pdf'], captiontext, sprintf('ple%i',count));
                
                lp(fid, '\n');
            end
        end
    end
    
    %% Confidence Intervals
    lp(fid, '\\clearpage\n');
    lp(fid, '\\section{Confidence intervals for the model parameters}\n');
    
    N = 30;
    ntables = 0;
    for j=1:length(S.pleGlobals.p_labels)
        if(S.pleGlobals.q_fit(j))
            if(j<=length(S.pleGlobals.chi2s) && ~isempty(S.pleGlobals.chi2s{j}))
                ntables = ntables + 1;
            end
        end
    end
    ntables = ceil(ntables/N);
    
    if(ntables>1)
        lp(fid, 'In Table \\ref{conftable1} -- \\ref{conftable%i}, %2i\\%% confidence intervals for the estimated parameter values derived by the profile likelihood \\cite{Raue:2009ec} are given.', ntables, (1-S.pleGlobals.alpha_level)*100);
    else
        lp(fid, 'In Table \\ref{conftable1}, %2i\\%% confidence intervals for the estimated parameter values derived by the profile likelihood \\cite{Raue:2009ec} are given.', (1-S.pleGlobals.alpha_level)*100);
    end
    
    headstr = '\t\t\t & name & $\\hat\\theta$';
    if(S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
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
    for j=1:length(S.pleGlobals.p_labels)
        if(S.pleGlobals.q_fit(j))
            if(j<=length(S.pleGlobals.chi2s) && ~isempty(S.pleGlobals.chi2s{j}))
                count = count + 1;
                
                if(mod(count,N)==0)
                    lp(fid, '\t\t\t\\botrule');
                    lp(fid, '\t\t\\end{tabular}}');
                    lp(fid, '\t\t\\mycaption{Confidence intervals for the estimated parameter values derived by the profile likelihood}{conftable%i}', counttab);
                    lp(fid, '\t\t{$\\hat\\theta$ indicates the estimated optimal parameter value.');
                    if(S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
                        lp(fid, '\t\t$\\sigma^{-}_{ptw}$ and $\\sigma^{+}_{ptw}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
                        lp(fid, '\t\t$\\sigma^{-}_{sim}$ and $\\sigma^{+}_{sim}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
                    elseif(S.pleGlobals.plot_point && ~S.pleGlobals.plot_simu)
                        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
                    elseif(~S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
                        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
                    end
                    lp(fid, '}');
                    lp(fid, '\t\\doendcenter');
                    lp(fid, '\t\\end{table}');
                    counttab = counttab + 1;
                    
                    lp(fid, '\t\\begin{table}');
                    lp(fid, '\t\\dobegincenter');
                    lp(fid, '\t{\\footnotesize');
                    lp(fid, '\t\t\\begin{tabular}{lllllll}');
                    lp(fid, '\t\t\t\\toprule');
                    lp(fid, headstr);
                    lp(fid, '\t\t\t\\midrule');
                end
                
                lp(fid, '\t\t\t%i & %s & %+8.3f & ', j, ...
                    strrep(S.pleGlobals.p_labels{j},'_','\_'), S.pleGlobals.p(j));
                if(S.pleGlobals.plot_point)
                    lp(fid, '%+8.3f & %+8.3f &', S.pleGlobals.conf_lb_point(j), S.pleGlobals.conf_ub_point(j));
                end
                if(S.pleGlobals.plot_simu)
                    lp(fid, '%+8.3f & %+8.3f', S.pleGlobals.conf_lb(j), S.pleGlobals.conf_ub(j));
                end
                lp(fid, ' \\\\');
            end
        end
    end
    lp(fid, '\t\t\t\\botrule');
    lp(fid, '\t\t\\end{tabular}}');
    lp(fid, '\t\t\\mycaption{Confidence intervals for the estimated parameter values derived by the profile likelihood}{conftable%i}', counttab);
    lp(fid, '\t\t{$\\hat\\theta$ indicates the estimated optimal parameter value.');
    if(S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
        lp(fid, '\t\t$\\sigma^{-}_{ptw}$ and $\\sigma^{+}_{ptw}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
        lp(fid, '\t\t$\\sigma^{-}_{sim}$ and $\\sigma^{+}_{sim}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
    elseif(S.pleGlobals.plot_point && ~S.pleGlobals.plot_simu)
        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
    elseif(~S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
        lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
    end
    lp(fid, '}');
    lp(fid, '\t\\doendcenter');
    lp(fid, '\t\\end{table}');
    
    %% Confidence intervals of model trajectories
    %     \subsection{Confidence intervals of the predicted model dynamics} \label{obsanalysis}
    % TODO
end

lp(fid, '\\bibliographystyle{plain}');
lp(fid, '\\bibliography{lib}');

lp(fid, '\\end{document}');



fclose(fid);
fprintf('done\n');

%% pdflatex
if(isunix)
    fprintf('pdflatex, file %s...', fname);
    cd(savePath);
    eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
    eval(['!bibtex ' fnamebib ' > log_bibtex.txt']);
    eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
    eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
    cd('../../..');
    try
        copyfile([savePath '/' 'report.pdf'], [savePath '/' sprintf('report_%s.pdf', datestr(now,30))])
        fprintf('done\n');
    catch
        fprintf('report.pdf was not written correctly\n');
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
        copyfile([savePath '/' 'report.pdf'], [savePath '/' sprintf('report_%s.pdf', datestr(now,30))])
        fprintf('done\n');
    catch
        fprintf('report.pdf was not written correctly\n');
    end
end
setenv('LD_LIBRARY_PATH', library_path);

function lp(varargin)
if(nargin>2)
    fprintf(varargin{1}, sprintf('%s\n', varargin{2}), varargin{3:end});
else
    fprintf(varargin{1}, sprintf('%s\n', varargin{2}));
end

function lpfigure(fid, textwidth, figpath, figcaption, figlabel)
lp(fid, '\\begin{figure}');
lp(fid, '\\begin{center}');
lp(fid, '\\includegraphics[width=%f\\textwidth]{%s} \\caption{%s} \\label{%s}', textwidth, figpath, figcaption, figlabel);
lp(fid, '\\end{center}');
lp(fid, '\\end{figure}');

function hmstimestr = secToHMS(seconds)
hours = floor(seconds/3600);
seconds = seconds - hours*3600;
minutes = floor(seconds/60);
seconds = seconds - minutes*60;
hmstimestr = sprintf('%02i:%02i:%05.2f', hours, minutes, seconds);

function str = myNameTrafo(str)
str = strrep(str, '_', '\_');
str = strrep(str, '%', '\%');

function str = myFormulas(str, jm)
global ar

varlist = symvar(str)';
svarlist = sym(varlist);
shortlist = {};
for j=1:length(varlist)
    shortlist{j} = sprintf('vj%ijv', j);
end
sshortlist = sym(shortlist);

strsym = sym(str);
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

str = strrep(str, '_', '\_');
str = strrep(str, '\,', ' \cdot ');





function fprintnumtab(fid, num)
fprintf(fid, '& %s ', sprintf('%g', num));

% better subs
function out = mysubs(in, old, new)
if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
    matVer = ver('MATLAB');
    if(str2double(matVer.Version)>=8.1)
        out = subs(in, old(:), new(:));
    else
        out = subs(in, old(:), new(:), 0);
    end
else
    out = in;
end

function mod = getModifierStr(jm,jv,useNeg,str,sources,targets)

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

% %% Residuals
% if(isfield(ar.report, 'residualPlotPath'))
% 	lp(fid, '\\subsection{Residual plots}');
% 	copyfile([ar.report.residualPlotPath '.pdf'], ...
% 		['./' savePath '/residuals.pdf'])
% 	lp(fid, '\\begin{center}');
% 	lp(fid, '\\includegraphics[width=\\textwidth]{residuals.pdf}');
% 	lp(fid, '\\end{center}');
% end