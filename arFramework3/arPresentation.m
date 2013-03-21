% create presentation

function arPresentation(doEquations)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end
if(~exist('doEquations','var'))
	doEquations = true;
end

savePath = [arSave '/Latex'];

if(~exist([cd '/' savePath], 'dir'))
    mkdir([cd '/' savePath])
end

% latex file
fname = 'presentation.tex';
fnamebib = 'presentation.aux';
fprintf('writing latex, file %s...', fname);
fid = fopen([savePath '/' fname], 'w');

%% Head
lp(fid, '\\nonstopmode');
lp(fid, '\\documentclass[8pt, oneside]{beamer}');
lp(fid, '\\usetheme{Boadilla}');

lp(fid, '\\usepackage[utf8]{inputenc}');
lp(fid, '\\usepackage{amsmath, amsthm, amssymb}');
lp(fid, '\\usepackage{listings} ');
lp(fid, '\\usepackage{epsfig} ');
lp(fid, '\\usepackage{graphicx} ');
lp(fid, '\\usepackage{rotating} ');
lp(fid, '\\usepackage{lscape}');
lp(fid, '\\usepackage{color}');
% lp(fid, '\\usepackage{natbib}');
lp(fid, '\\usepackage{lmodern} ');
lp(fid, '\\usepackage{xcolor} ');
lp(fid, '\\usepackage[normal,footnotesize,bf]{caption}');
% lp(fid, '\\usepackage{subfig}');
lp(fid, '\\usepackage{sidecap}');
lp(fid, '\\captionsetup[subfloat]{position=top}');

% lp(fid, '\\setlength{\\textheight}{22 cm}');
% lp(fid, '\\setlength{\\textwidth}{16 cm}');
% lp(fid, '\\setlength{\\topmargin}{-1.5 cm}');
% lp(fid, '\\setlength{\\hoffset}{-2 cm}');
% lp(fid, '\\usepackage[colorlinks=true, linkcolor=blue, citecolor=blue, filecolor=blue, urlcolor=blue]{hyperref}\n');

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
lp(fid, '\\title{Modeling Report}');
lp(fid, '\\author{%s}', ar.config.username);
lp(fid, '\\institute{%s}', 'University of Freiburg');
lp(fid, '\\date{%s}', datestr(now));
lp(fid, '\\maketitle\n');

% lp(fid, '\\tableofcontents');

N = 10;

for jm=1:length(ar.model)
    lp(fid, '\\section{Model: %s}\n', myNameTrafo(ar.model(jm).name));
    lp(fid, '\\subsection{General Definition}\n');
    
    if(doEquations)
        lp(fid, '\\begin{frame}');
        lp(fid, '\\frametitle{Model: %s}\n', myNameTrafo(ar.model(jm).name));
        if(~isempty(ar.model(jm).description))
            lp(fid, '\\begin{itemize}');
            for jdes=1:length(ar.model(jm).description)
                lp(fid, '\\item %s\\\\', strrep(strrep(ar.model(jm).description{jdes}, '%', '\%'), '_', '\_'));
            end
            lp(fid, '\\end{itemize}');
        end
        
        lp(fid, 'The model consists of %i external inputs, %i dynamical variables indicated by square brackets and %i reactions', ...
            length(ar.model(jm).u), length(ar.model(jm).x), length(ar.model(jm).fv));
        lp(fid, ' and was evaluated for %i experimental conditions.', length(ar.model(jm).condition));
        if(isfield(ar.model(jm), 'savePath_Graph') && ~isempty(ar.model(jm).savePath_Graph))
            lp(fid, 'The model structure is depicted in Figure \\ref{%s}.', [ar.model(jm).name '_graph']);
        end
        if(isfield(ar.model(jm), 'data'))
            if(ar.config.fiterrors == 1)
                lp(fid, 'In total %i parameters are estimated from the experimental data, yielding a value of the objective function $-2 \\log(L) = %g$ for a total of %i data points.', ...
                    sum(ar.qFit==1), 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit, ar.ndata);
            else
                lp(fid, 'In total %i parameters are estimated from the experimental data, yielding a value of the objective function $\\chi^2 = %g$ for a total of %i data points.', ...
                    sum(ar.qFit==1), ar.chi2, ar.ndata);
            end
        end
        % 	lp(fid, 'The estimated parameter values are given in Section \\ref{estimatedparameters}.\\\\');
        lp(fid, '\\end{frame}');
    end

    
	% Structure
    savePath_Graph = [arSave '/Figures/Network' sprintf('/%s.pdf', ar.model(jm).name)];
    if(exist(savePath_Graph,'file'))
        lp(fid, '\\begin{frame}');
        copyfile(savePath_Graph, [savePath '/' ar.model(jm).name '_graph.pdf'])
% 		captiontext = sprintf('\\textbf{Network representation of the model %s}\\\\ ', myNameTrafo(ar.model(jm).name));
% 		if(~isempty(ar.model(jm).u))
% 			captiontext = [captiontext 'Diamond shaped nodes correspond to model inputs, see Equation '];
% 			if(length(ar.model(jm).u)==1)
% 				captiontext = [captiontext '\ref{' sprintf('%s_input%i', ar.model(jm).name, 1) '}. '];
% 			else
% 				captiontext = [captiontext '\ref{' sprintf('%s_input%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_input%i', ar.model(jm).name, length(ar.model(jm).u)) '}. '];
% 			end
% 		end
% 		captiontext = [captiontext 'Ellipsoid shaped nodes correspond to dynamical variables described by the ODE system, see Equation '];
% 		captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
% 		captiontext = [captiontext 'Black arrows and box shaped nodes indicate reactions with corresponding rate equations given in Equation '];
% 		captiontext = [captiontext '\ref{' sprintf('%s_flux%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_flux%i', ar.model(jm).name, length(ar.model(jm).fv)) '}. '];
% 		captiontext = [captiontext 'Red T-shaped arrows indicated inhibitorial influence on a reaction and blue O-shaped arrows indicate catalysing influence on a reaction. '];
		lpfigure(fid, 0.9, [ar.model(jm).name '_graph.pdf'], [], [ar.model(jm).name '_graph']);
        lp(fid, '\\end{frame}');
    end
    
    % Inputs
    if(doEquations && ~isempty(ar.model(jm).u))
        lp(fid, '\\begin{frame}');
        lp(fid, '\\noindent The model dynamics depends on the external inputs:');
        lp(fid, '{\\tiny');
        lp(fid, '\\begin{eqnarray}');
        for ju=1:length(ar.model(jm).fu)
            if(~isempty(ar.model(jm).fu{ju}))
                strfu = ['& = & ' myFormulas(ar.model(jm).fu{ju}, jm)];
            else
                strfu = '';
            end
            if(ju==length(ar.model(jm).fu) || mod(ju,N)==0)
                lp(fid, '\t%s(%s) %s \\label{%s}', ...
					myFormulas(ar.model(jm).u{ju}, jm), myFormulas(ar.model(jm).t, jm), strfu, sprintf('%s_input%i', ar.model(jm).name, ju));
            else
                lp(fid, '\t%s(%s) %s \\label{%s} \\\\', ...
					myFormulas(ar.model(jm).u{ju}, jm), myFormulas(ar.model(jm).t, jm), strfu, sprintf('%s_input%i', ar.model(jm).name, ju));
            end
            if(mod(ju,N)==0 && ju<length(ar.model(jm).fu))
                lp(fid, '\\end{eqnarray}}\n');
                lp(fid, '\\end{frame}');
                lp(fid, '\\begin{frame}');
                lp(fid, '\\quad');
                lp(fid, '{\\tiny');
                lp(fid, '\\begin{eqnarray}');
            end
        end
        lp(fid, '\\end{eqnarray}}\n\n');
        lp(fid, '\\end{frame}');
    end
    
    % Fluxes
    if(~isempty(ar.model(jm).fv))
        lp(fid, '\\begin{frame}');
        if(doEquations)
            lp(fid, '\\noindent The rate equations corresponding to the reactions included in the model are give by:');
        else
            lp(fid, '\\noindent \\quad \\\\');
        end
        
%         lp(fid, '{\\tiny');
%         lp(fid, '\\begin{eqnarray}');
%         for jv = 1:length(ar.model(jm).fv);
%             strtmp = myFormulas(ar.model(jm).fv{jv}, jm);
%             if(jv == length(ar.model(jm).fv) || mod(jv,N)==0)
%                 lp(fid, '\tv_{%i} & = & %s \\label{%s}', jv, strtmp, sprintf('%s_flux%i', ar.model(jm).name, jv));
%             else
%                 lp(fid, '\tv_{%i} & = & %s \\label{%s} \\\\', jv, strtmp, sprintf('%s_flux%i', ar.model(jm).name, jv));
%             end
%             
%             if(mod(jv,N)==0 && jv<length(ar.model(jm).fv))
%                 lp(fid, '\\end{eqnarray}}\n');
%                 lp(fid, '\\end{frame}');
%                 lp(fid, '\\begin{frame}');
%                 lp(fid, '\\quad');
%                 lp(fid, '{\\tiny');
%                 lp(fid, '\\begin{eqnarray}');
%             end
%         end
%         lp(fid, '\\end{eqnarray}}\n\n');
        
        for jv = 1:length(ar.model(jm).fv);
            lp(fid, '{\\tiny');
            str = sprintf('v_{%i}: ',jv);
            ssource = find(ar.model(jm).N(:,jv)<0);
            for js = 1:length(ssource)
                if(abs(ar.model(jm).N(ssource(js),jv)) == 1)
                    if(js>1)
                        str = [str ' + ' myNameTrafo(ar.model(jm).x{ssource(js)})];
                    else
                        str = [str myNameTrafo(ar.model(jm).x{ssource(js)})];
                    end
                else
                    if(js>1)
                        str = [str ' + ' num2str(abs(ar.model(jm).N(ssource(js),jv))) ...
                            '*' myNameTrafo(ar.model(jm).x{ssource(js)})];
                    else
                        str = [str num2str(abs(ar.model(jm).N(ssource(js),jv))) ...
                            '*' myNameTrafo(ar.model(jm).x{ssource(js)})];
                    end
                end
            end
            if(isempty(ssource))
                str = [str '$\emptyset$'];
            end
            str = [str ' ~ $\rightarrow$ '];
            starget = find(ar.model(jm).N(:,jv)>0);
            for js = 1:length(starget)
                if(js>1)
                    str = [str ' + ' myNameTrafo(ar.model(jm).x{starget(js)})];
                else
                    str = [str myNameTrafo(ar.model(jm).x{starget(js)})];
                end
            end
            if(isempty(starget))
                str = [str '$\emptyset$'];
            end
            fprintf(fid, strrep(str, '\','\\'));
            fprintf(fid, '\\quad : ');
            
            
            strtmp = myFormulas(ar.model(jm).fv{jv}, jm);
            lp(fid, '$%s$}\\\\\\\n', strtmp);        
            
            if(mod(jv,N)==0 && jv<length(ar.model(jm).fv))
                lp(fid, '\\end{frame}');
                lp(fid, '\\begin{frame}');
                lp(fid, '\\noindent \\quad \\\\');
            end
        end

        lp(fid, '\\end{frame}');
    end
    
    % ODE system
    if(doEquations && ~isempty(ar.model(jm).x))
        lp(fid, '\\begin{frame}');
        lp(fid, '\\noindent The ODE system determining the time evolution of the dynamical variables is given by:');
        lp(fid, '{\\tiny');
        lp(fid, '\\begin{eqnarray}');
        for jx=1:size(ar.model(jm).N, 1) % for every species jx
            strtmp = '';
            if(~isempty(ar.model(jm).c))
                qinfluxwitheducts = ar.model(jm).N(jx,:) > 0 & sum(ar.model(jm).N < 0,1) > 0;
                eductcompartment = zeros(size(qinfluxwitheducts));
                for jj=find(qinfluxwitheducts)
                    eductcompartment(jj) = unique(ar.model(jm).cLink(ar.model(jm).N(:,jj)<0));
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
                lp(fid, '\t\\mathrm{d}%s/\\mathrm{d}t & = & %s \\label{%s}', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_ode%i', ar.model(jm).name, jx));
            else
                lp(fid, '\t\\mathrm{d}%s/\\mathrm{d}t & = & %s \\label{%s} \\\\', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_ode%i', ar.model(jm).name, jx));
            end
            
            if(mod(jx,N)==0 && jx<size(ar.model(jm).N, 1))
                lp(fid, '\\end{eqnarray}}\n');
                lp(fid, '\\end{frame}');
                lp(fid, '\\begin{frame}');
                lp(fid, '\\quad');
                lp(fid, '{\\tiny');
                lp(fid, '\\begin{eqnarray}');
            end
        end
        lp(fid, '\\end{eqnarray}}\n\n');
        
        lp(fid, '\\noindent The ODE system was solved by a parallelized implementation of the CVODES algorithm \\cite{Hindmarsh:2005fb}.');
        lp(fid, 'It also supplies the parameter sensitivities utilized for parameter estimation.\\\\\n\n');
        lp(fid, '\\end{frame}');
    end
    
    % ODE initials
    if(doEquations && ~isempty(ar.model(jm).x))
        lp(fid, '\\begin{frame}');
        lp(fid, '\\noindent The initial conditions for the ODE system are given by:');
        lp(fid, '{\\tiny');
        lp(fid, '\\begin{eqnarray}');
        for jx=1:size(ar.model(jm).N, 1) % for every species jx
            strtmp = myFormulas(ar.model(jm).px0{jx}, jm);
            if(jx==size(ar.model(jm).N, 1) || mod(jx,N)==0)
                lp(fid, '\t%s(0) & = & %s \\label{%s}', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_init%i', ar.model(jm).name, jx));
            else
                lp(fid, '\t%s(0) & = & %s \\label{%s} \\\\', myFormulas(ar.model(jm).x{jx}, jm), strtmp, sprintf('%s_init%i', ar.model(jm).name, jx));
            end
            
            if(mod(jx,N)==0 && jx<size(ar.model(jm).N, 1))
                lp(fid, '\\end{eqnarray}}\n');
                lp(fid, '\\end{frame}');
                lp(fid, '\\begin{frame}');
                lp(fid, '\\quad');
                lp(fid, '{\\tiny');
                lp(fid, '\\begin{eqnarray}');
            end
        end
        lp(fid, '\\end{eqnarray}}\n\n');
        lp(fid, '\\end{frame}');
    end
    
    % Conditions
    if(doEquations)
        ccount = 1;
        for jp=1:length(ar.model(jm).fp)
            if(~strcmp(ar.model(jm).p{jp}, ar.model(jm).fp{jp}))
                if(ccount==1)
                    lp(fid, '\\begin{frame}');
                    lp(fid, '\\noindent The ODE system is modified by the following parameter transformations:');
                    lp(fid, '{\\tiny');
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
                    lp(fid, '\\end{eqnarray}}\n');
                    lp(fid, '\\end{frame}');
                    lp(fid, '\\begin{frame}');
                    lp(fid, '\\quad');
                    lp(fid, '{\\tiny');
                    lp(fid, '\\begin{eqnarray}');
                end
                ccount = ccount + 1;
            end
            if(ccount>1 && jp==length(ar.model(jm).fp))
                lp(fid, '\\end{eqnarray}}\n\n');
                lp(fid, '\\end{frame}');
            end
        end
    end
    
%% Experiments    
    
    for jplot=1:length(ar.model(jm).plot)
        lp(fid, '\\subsection{Experiment: %s}\n', myNameTrafo(ar.model(jm).plot(jplot).name));
        jd = ar.model(jm).plot(jplot).dLink(1);

        if(isfield(ar.model(jm), 'data'))
            if(doEquations)
                lp(fid, '\\begin{frame}');
                lp(fid, '\\frametitle{Experiment: %s}', myNameTrafo(ar.model(jm).plot(jplot).name));
                if(~isempty(ar.model(jm).data(jd).description))
                    lp(fid, '\\begin{itemize}');
                    for jdes=1:length(ar.model(jm).data(jd).description)
                        lp(fid, '\\item %s\\\\', strrep(strrep(ar.model(jm).data(jd).description{jdes}, '%', '\%'), '_', '\_'));
                    end
                    lp(fid, '\\end{itemize}');
                end
                %             lp(fid, '\\noindent The agreement of the model outpxuts and the experimental data, given in Table \\ref{%s_data}, ', ar.model(jm).plot(jplot).name);
                lp(fid, '\\noindent The agreement of the model outpxuts and the experimental data, ');
                
                if(ar.config.fiterrors == 1)
                    lp(fid, 'yields a value of the objective function $-2 \\log(L) = %g$ for %i data points in this data set.', ...
                        2*ar.model(jm).plot(jplot).ndata*log(sqrt(2*pi)) + ar.model(jm).plot(jplot).chi2, ar.model(jm).plot(jplot).ndata);
                else
                    lp(fid, 'yields a value of the objective function $\\chi^2 = %g$ for %i data points in this data set.', ...
                        ar.model(jm).plot(jplot).chi2, ar.model(jm).plot(jplot).ndata);
                end
                lp(fid, '\\end{frame}');
            end
            
            % observations
            if(doEquations)
                lp(fid, '\\begin{frame}');
                lp(fid, '\\noindent The model outputs available in this data set are defined by:');
                lp(fid, '{\\tiny');
                lp(fid, '\\begin{eqnarray}');
                for jy=1:length(ar.model(jm).data(jd).fy)
                    if(ar.model(jm).data(jd).logfitting(jy))
                        strtmp = ['\log_{10}(' myFormulas(ar.model(jm).data(jd).fy{jy}, jm) ')'];
                    else
                        strtmp = myFormulas(ar.model(jm).data(jd).fy{jy}, jm);
                    end
                    if(jy==length(ar.model(jm).data(jd).fy) || mod(jy,N)==0)
                        lp(fid, '\t%s & = & %s \\label{%s}', myFormulas(ar.model(jm).data(jd).y{jy}, jm), strtmp, ...
                            sprintf('%s_obs%i', ar.model(jm).plot(jplot).name, jy));
                    else
                        lp(fid, '\t%s & = & %s \\label{%s} \\\\', myFormulas(ar.model(jm).data(jd).y{jy}, jm), strtmp, ...
                            sprintf('%s_obs%i', ar.model(jm).plot(jplot).name, jy));
                    end
                    
                    if(mod(jy,N)==0 && jy<length(ar.model(jm).data(jd).fy))
                        lp(fid, '\\end{eqnarray}}\n');
                        lp(fid, '\\end{frame}');
                        lp(fid, '\\begin{frame}');
                        lp(fid, '\\quad');
                        lp(fid, '{\\tiny');
                        lp(fid, '\\begin{eqnarray}');
                    end
                end
                lp(fid, '\\end{eqnarray}\n}\n\n');
                lp(fid, '\\end{frame}');
            end
            
            % error functions
            if(doEquations)
                lp(fid, '\\begin{frame}');
                lp(fid, '\\noindent The error model that describes the measurement noise for each model output is given by:');
                lp(fid, '{\\tiny');
                lp(fid, '\\begin{eqnarray}');
                for jy=1:length(ar.model(jm).data(jd).fystd)
                    if(jy==length(ar.model(jm).data(jd).fystd) || mod(jy,N)==0)
                        lp(fid, '\t%s & = & %s \\label{%s}', myFormulas(ar.model(jm).data(jd).y{jy}, jm), ...
                            myFormulas(ar.model(jm).data(jd).fystd{jy}, jm), sprintf('%s_err%i', ar.model(jm).plot(jplot).name, jy));
                    else
                        lp(fid, '\t%s & = & %s \\label{%s} \\\\', myFormulas(ar.model(jm).data(jd).y{jy}, jm), ...
                            myFormulas(ar.model(jm).data(jd).fystd{jy}, jm), sprintf('%s_err%i', ar.model(jm).plot(jplot).name, jy));
                    end
                    
                    if(mod(jy,N)==0 && jy<length(ar.model(jm).data(jd).fy))
                        lp(fid, '\\end{eqnarray}}\n');
                        lp(fid, '\\end{frame}');
                        lp(fid, '\\begin{frame}');
                        lp(fid, '\\quad');
                        lp(fid, '{\\tiny');
                        lp(fid, '\\begin{eqnarray}');
                    end
                end
                lp(fid, '\\end{eqnarray}\n}\n\n');
                lp(fid, '\\end{frame}');
            end
        end
        
        % Inputs
        if(doEquations && ~isempty(ar.model(jm).u))
            ucount = 1;
            for ju=1:length(ar.model(jm).data(jd).fu)
                if(~strcmp(ar.model(jm).data(jd).fu{ju}, ar.model(jm).fu{ju}))
                    if(ucount==1)
                        lp(fid, '\\begin{frame}');
                        lp(fid, ['\\noindent To evaluate the ODE system of Equation \\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \\ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}']);
                        lp(fid, ' for the conditions in this experiment, the following external inputs are given:');
                        lp(fid, '{\\tiny');
                        lp(fid, '\\begin{eqnarray}');
                    end
                    if(ucount==length(ar.model(jm).data(jd).fu) || mod(ucount,N)==0)
                        lp(fid, '\t%s(%s) & = & %s', myFormulas(ar.model(jm).u{ju}, jm), myFormulas(ar.model(jm).t, jm), ...
                            myFormulas(ar.model(jm).data(jd).fu{ju}, jm));
                    else
                        lp(fid, '\t%s(%s) & = & %s \\\\', myFormulas(ar.model(jm).u{ju}, jm), myFormulas(ar.model(jm).t, jm), ...
                            myFormulas(ar.model(jm).data(jd).fu{ju}, jm));
                    end
                    
                    if(mod(ucount,N)==0 && ucount<length(ar.model(jm).data(jd).fu))
                        lp(fid, '\\end{eqnarray}}\n');
                        lp(fid, '\\end{frame}');
                        lp(fid, '\\begin{frame}');
                        lp(fid, '\\quad');
                        lp(fid, '{\\tiny');
                        lp(fid, '\\begin{eqnarray}');
                    end
                    ucount = ucount + 1;
                end
                if(ucount>1 && ju==length(ar.model(jm).data(jd).fu))
                    lp(fid, '\\end{eqnarray}\n}\n\n');
                    lp(fid, '\\end{frame}');
                end
            end
        end
        
        if(isfield(ar.model(jm), 'data'))
            if(doEquations)
                % conditions
                ccount = 1;
                for jp=1:length(ar.model(jm).data(jd).fp)
                    
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
                        qdyn = ismember(ar.model(jm).p, ar.model(jm).data(jd).pold{jp});
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
                                lp(fid, '\\begin{frame}');
                                lp(fid, '\\noindent To evaluate the ODE system the following parameter transformations are applied:');
                                lp(fid, '{\\tiny');
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
                                lp(fid, '\\end{array}}\n');
                                lp(fid, '\\end{displaymath}\n');
                                lp(fid, '\\end{frame}');
                                lp(fid, '\\begin{frame}');
                                lp(fid, '\\quad');
                                lp(fid, '\\begin{displaymath}');
                                lp(fid, '\\begin{array}{%s}', arraystr);
                            end
                            ccount = ccount + 1;
                        end
                    end
                    
                    if(ccount>1 && jp==length(ar.model(jm).data(jd).fp))
                        lp(fid, '\\end{array}');
                        lp(fid, '\\end{displaymath}\n}\n\n');
                        lp(fid, '\\end{frame}');
                    end
                end
            end
            
            % plots
            if(isfield(ar.model(jm).plot(jplot), 'savePath_FigY') && ~isempty(ar.model(jm).plot(jplot).savePath_FigY))
                lp(fid, '\\begin{frame}');
                copyfile([ar.model(jm).plot(jplot).savePath_FigY '.pdf'], ...
                    [savePath '/' ar.model(jm).plot(jplot).name '_y.pdf']);
                
%                 lp(fid, 'The model outputs and the experimental data is show in Figure \\ref{%s}.', [ar.model(jm).plot(jplot).name '_y']);
%                 captiontext = sprintf('\\textbf{Agreement of model outputs and experimental data for the experiment %s}\\\\', myNameTrafo(ar.model(jm).plot(jplot).name));
%                 captiontext = [captiontext 'The displayed model outputs (solid lines) are defined by Equation '];
%                 captiontext = [captiontext '\ref{' sprintf('%s_obs%i', ar.model(jm).plot(jplot).name, 1) '} -- \ref{' sprintf('%s_obs%i', ar.model(jm).plot(jplot).name, length(ar.model(jm).data(jd).fy)) '}. '];
%                 captiontext = [captiontext 'The error model that describes the measurement noise for each model output is indicated by shades around the model outputs and is given by Equation '];
%                 captiontext = [captiontext '\ref{' sprintf('%s_err%i', ar.model(jm).plot(jplot).name, 1) '} -- \ref{' sprintf('%s_err%i', ar.model(jm).plot(jplot).name, length(ar.model(jm).data(jd).fy)) '}. '];
                lpfigure(fid, 0.9, [ar.model(jm).plot(jplot).name '_y.pdf'], [], [ar.model(jm).plot(jplot).name '_y']);
                lp(fid, '\\end{frame}');
            end
            
%             % experimental data
%             
%             headstr = '';
%             headtab = '';
%             unitstr = '';
%             % time
%             unitstr = [unitstr sprintf('%s [%s] ', myNameTrafo(ar.model(jm).data(jd).tUnits{3}), ...
%                 myNameTrafo(ar.model(jm).data(jd).tUnits{2}))];
%             headstr = [headstr ' '];
%             headtab = [headtab 'r'];
%             % conditions
%             if(~isempty(ar.model(jm).data(jd).condition))
%                 for jp = 1:length(ar.model(jm).data(jd).condition)
%                     headstr = [headstr sprintf('& %s ', myNameTrafo(ar.model(jm).data(jd).condition(jp).parameter))];
%                     unitstr = [unitstr '& '];
%                     headtab = [headtab 'r'];
%                 end
%             end
%             % y & ystd headers
%             for jy=1:length(ar.model(jm).data(jd).y)
%                 headstr = [headstr sprintf('& %s ', myNameTrafo(ar.model(jm).data(jd).y{jy}))];
%                 unitstr = [unitstr sprintf('& %s [%s] ', myNameTrafo(ar.model(jm).data(jd).yUnits{jy,3}), ...
%                     myNameTrafo(ar.model(jm).data(jd).yUnits{jy,2}))];
%                 headtab = [headtab 'r'];
%                 if(ar.config.fiterrors == -1)
%                     headstr = [headstr sprintf('& %s\_std ', myNameTrafo(ar.model(jm).data(jd).y{jy}))];
%                     unitstr = [unitstr sprintf('& %s [%s] ', myNameTrafo(ar.model(jm).data(jd).yUnits{jy,3}), ...
%                         myNameTrafo(ar.model(jm).data(jd).yUnits{jy,2}))];
%                     headtab = [headtab 'r'];
%                 end
%             end
%             lp(fid, '\t\\begin{table}');
%             lp(fid, '\t\\dobegincenter');
%             lp(fid, '\t{\\tiny');
%             lp(fid, '\t\t\\begin{tabular}{%s}', headtab);
%             lp(fid, '\t\t\t\\toprule');
%             lp(fid, '\t\t\t %s \\\\', headstr);
%             lp(fid, '\t\t\t %s \\\\', unitstr);
%             lp(fid, '\t\t\t\\midrule');
%             
%             for jd2 = ar.model(jm).plot(jplot).dLink
%                 for j=1:length(ar.model(jm).data(jd2).tExp)
%                     fprintf(fid, '%s ', sprintf('%f', ar.model(jm).data(jd2).tExp(j)));
%                     
%                     % conditions
%                     if(~isempty(ar.model(jm).data(jd2).condition))
%                         for jp = 1:length(ar.model(jm).data(jd2).condition)
%                             fprintf(fid, '& %s ', ar.model(jm).data(jd2).condition(jp).value);
%                         end
%                     end
%                     
%                     % y data
%                     for jj=1:size(ar.model(jm).data(jd2).yExp,2)
%                         if(ar.model(jm).data(jd2).logfitting(jj))
%                             fprintnumtab(fid, 10^ar.model(jm).data(jd2).yExp(j,jj));
%                         else
%                             fprintnumtab(fid, ar.model(jm).data(jd2).yExp(j,jj));
%                         end
%                     end
%                     
%                     % ystd data
%                     if(ar.config.fiterrors == -1)
%                         for jj=1:size(ar.model(jm).data(jd2).yExp,2)
%                             fprintnumtab(fid, ar.model(jm).data(jd2).yExpStd(j,jj));
%                         end
%                     end
%                     fprintf(fid, '\\\\\n');
%                 end
%             end
%             
%             lp(fid, '\t\t\t\\botrule');
%             lp(fid, '\t\t\\end{tabular}}');
%             lp(fid, '\t\t\\mycaption{Experimental data for the experiment %s}{%s_data}{}', myNameTrafo(ar.model(jm).plot(jplot).name), ar.model(jm).plot(jplot).name);
%             lp(fid, '\t\\doendcenter');
%             lp(fid, '\t\\end{table}');
        end
        
        % trajectories
        if(isfield(ar.model(jm).plot(jplot), 'savePath_FigX'))
            lp(fid, '\\begin{frame}');
            lp(fid, 'The trajectories of the dynamical variables and external inputs that correspond to the experimental conditions in this experiment are shown in Figure \\ref{%s}.', [ar.model(jm).plot(jplot).name '_x']);
            copyfile([ar.model(jm).plot(jplot).savePath_FigX '.pdf'], ...
                [savePath '/' ar.model(jm).plot(jplot).name '_x.pdf'])  
%             captiontext = sprintf('\\textbf{Trajectories of the dynamical variables and external inputs for the experiment %s}\\\\', myNameTrafo(ar.model(jm).plot(jplot).name));
%             captiontext = [captiontext 'The dynamical behaviour is determined by the ODE system, see Equation '];
%             captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
            lpfigure(fid, 1, [ar.model(jm).plot(jplot).name '_x.pdf'], [], [ar.model(jm).plot(jplot).name '_x']);
            lp(fid, '\\end{frame}');
        end
        if(isfield(ar.model(jm).plot(jplot), 'savePath_FigV'))
            lp(fid, '\\begin{frame}');
            lp(fid, 'The reaction fluxes that correspond to the experimental conditions in this experiment are shown in Figure \\ref{%s}.', [ar.model(jm).plot(jplot).name '_v']);
            copyfile([ar.model(jm).plot(jplot).savePath_FigV '.pdf'], ...
                [savePath '/' ar.model(jm).plot(jplot).name '_v.pdf'])
%             captiontext = sprintf('\\textbf{Reaction fluxes for the experiment %s}\\\\', myNameTrafo(ar.model(jm).plot(jplot).name));
%             captiontext = [captiontext 'The dynamical behaviour is determined by the ODE system, see Equation '];
%             captiontext = [captiontext '\ref{' sprintf('%s_ode%i', ar.model(jm).name, 1) '} -- \ref{' sprintf('%s_ode%i', ar.model(jm).name, length(ar.model(jm).x)) '}. '];
            lpfigure(fid, 1, [ar.model(jm).plot(jplot).name '_v.pdf'], [], [ar.model(jm).plot(jplot).name '_v']);
            lp(fid, '\\end{frame}');
        end
    end
end

% %% Parameters
% lp(fid, '\\clearpage\n');
% lp(fid, '\\section{Estimated model parameters} \\label{estimatedparameters}\n');
% lp(fid, 'The model parameter were estimated by maximum likelihood estimation applying the MATLAB lsqnonlin algorithm.');
% 
% N = 25;
% ntables = ceil(length(ar.p)/N);
% if(ntables>1)
%     lp(fid, 'In Table \\ref{paratable1} -- \\ref{paratable%i} the estimated parameter values are given.', ntables); 
% else
%     lp(fid, 'In Table \\ref{paratable1} the estimated parameter values are given.'); 
% end
% lp(fid, 'Parameters highlighted in red color indicate parameter values close to their bounds.');
% lp(fid, 'The parameter name prefix init\\_ indicates the initial value of a dynamic variable.');
% lp(fid, 'The parameter name prefix offset\\_ indicates a offset of the experimental data.');
% lp(fid, 'The parameter name prefix scale\\_ indicates a scaling factor of the experimental data.');
% lp(fid, 'The parameter name prefix sd\\_ indicates the magnitude of the measurement noise for a specific measurement.\\\\');
% 
% pTrans = ar.p;
% pTrans(ar.qLog10==1) = 10.^pTrans(ar.qLog10==1);
% 
% lp(fid, '\t\\begin{table}');
% lp(fid, '\t\\dobegincenter');
% lp(fid, '\t{\\tiny');
% lp(fid, '\t\t\\begin{tabular}{llllllll}');
% lp(fid, '\t\t\t\\toprule');
% lp(fid, '\t\t\t & name & $\\theta_{min}$ & $\\hat\\theta$ & $\\theta_{max}$ & log & non-log $\\hat \\theta$ & fitted \\\\');
% lp(fid, '\t\t\t\\midrule');
% count = 1;
% for j=1:length(ar.p)
%     if(mod(j,N)==0)
%         lp(fid, '\t\t\t\\botrule');
%         lp(fid, '\t\t\\end{tabular}}');
%         lp(fid, '\t\t\\mycaption{Estimated parameter values}{paratable%i}', count);
%         lp(fid, '{$\\hat \\theta$ indicates the estimated value of the parameters.');
%         lp(fid, '$\\theta_{min}$ and $\\theta_{max}$ indicate the upper and lower bounds for the parameters.');
%         lp(fid, 'The log-column indicates if the value of a parameter was log-transformed.');
%         lp(fid, 'If log = 1 the non-log-column indicates the non-logarithmic value of the estimate.');
%         lp(fid, 'The fitted-column indicates if the parameter value was estimated (1), was temporarily fixed (0) or if its value was fixed to a constant value (2).}');
%         lp(fid, '\t\\doendcenter');
%         lp(fid, '\t\\end{table}');
%         lp(fid, '\t\\begin{table}');
%         lp(fid, '\t\\dobegincenter');
%         lp(fid, '\t{\\tiny');
%         lp(fid, '\t\t\\begin{tabular}{llllllll}');
%         lp(fid, '\t\t\t\\toprule');
%         lp(fid, '\t\t\t & name & $\\theta_{min}$ & $\\hat\\theta$ & $\\theta_{max}$ & log & non-log $\\hat \\theta$ & fitted \\\\');
%         lp(fid, '\t\t\t\\midrule');
%         
%         count = count + 1;
%     end
%     
%     if(ar.qFit(j)==1)
%         if(ar.p(j) - ar.lb(j) < 0.1 || ar.ub(j) - ar.p(j) < 0.1)
%             lp(fid, '\t\t\t\\color{red}{%i} & \\color{red}{%s} & \\color{red}{%+8.0g} & \\color{red}{%+8.4f} & \\color{red}{%+8.0g} & \\color{red}{%i} & \\color{red}{%s} & \\color{red}{%i} \\\\', ...
%                 j, strrep(ar.pLabel{j},'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
%                 ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
%         else
%             lp(fid, '\t\t\t%i & %s & {%+8.0g} & {%+8.4f} & {%+8.0g} & %i & %s & %i \\\\', ...
%                 j, strrep(ar.pLabel{j},'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
%                 ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
%         end
%     else
%         lp(fid, '\t\t\t\\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%+8.0g} & \\color{mygray}{%+8.4f} & \\color{mygray}{%+8.0g} & \\color{mygray}{%i} & \\color{mygray}{%s} & \\color{mygray}{%i} \\\\', ...
%             j, strrep(ar.pLabel{j},'_','\_'), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), ...
%             ['$' strrep(sprintf('%+5.2e',pTrans(j)), 'e', '\cdot 10^{') '}$'], ar.qFit(j));
%     end
% end
% lp(fid, '\t\t\t\\botrule');
% lp(fid, '\t\t\\end{tabular}}');
% lp(fid, '\t\t\\mycaption{Estimated parameter values}{paratable%i}', count);
% lp(fid, '{$\\hat \\theta$ indicates the estimated value of the parameters.');
% lp(fid, '$\\theta_{min}$ and $\\theta_{max}$ indicate the upper and lower bounds for the parameters.');
% lp(fid, 'The log-column indicates if the value of a parameter was log-transformed.');
% lp(fid, 'If log = 1 the non-log-column indicates the non-logarithmic value of the estimate.');
% lp(fid, 'The fitted-column indicates if the parameter value was estimated (1), was temporarily fixed (0) or if its value was fixed to a constant value (2).}');
% lp(fid, '\t\\doendcenter');
% lp(fid, '\t\\end{table}');
% 
% %% PLE
% plePath = [arSave '/PLE'];
% if(exist(plePath,'dir'))
%     lp(fid, '\\clearpage\n');
% 	lp(fid, '\\section{Profile likelihood of model parameters}\n');
%     
% 	S = load([plePath '/results.mat']);
% 
%     lp(fid, 'In order to evaluate the identifiability of the model parameters and to assess confidence intervals');
%     lp(fid, 'the profile likelihood \\cite{Raue:2009ec} was calculated.');
%     lp(fid, 'The mean calculation time of the profile likelihood per parameter was %s $\\pm$ %s.', ...
%         secToHMS(mean(S.pleGlobals.timing(S.pleGlobals.q_fit(size(S.pleGlobals.timing))))), ...
%         secToHMS(std(S.pleGlobals.timing(S.pleGlobals.q_fit(size(S.pleGlobals.timing))))));
%     
% 	% Multiplot
%     if(isfield(S.pleGlobals, 'fighandel_multi'))
%         lp(fid, 'An overview is displayed in Figure \\ref{multi_plot}.');
%         
%         sourcestr = [plePath '/multi_plot.eps'];
%         targetstr = [savePath '/multi_plot.pdf'];
%         eval(['!ps2pdf  -dEPSCrop ' sourcestr ' ' targetstr]);
%                 
%         captiontext = '\textbf{Overview of the profile likelihood of the model parameters}\\';
%         captiontext = [captiontext 'The solide lines indicate the profile likelihood. '];
%         captiontext = [captiontext 'The broken lines indicate the threshold to assess confidence intervals. '];
%         captiontext = [captiontext 'The asterisk indicate the optimal parameter values. '];
%         lpfigure(fid, 1, 'multi_plot.pdf', captiontext, 'multi_plot');
%     end
%     
%     % Singleplots
%     if(isfield(S.pleGlobals, 'figPath'))
%         count = 0;
%         for j=1:length(S.pleGlobals.figPath)
%             if(~isempty(S.pleGlobals.chi2s{j}))
%                 count = count + 1;
%             end
%         end
%         
%         lp(fid, 'In Figure \\ref{ple1} -- \\ref{ple%i} the profile likelihood of each parameter is shown in more detail.', count); 
%         lp(fid, 'Also the functional relations to the remaining parameter are displayed.');
%         
%         count = 0;
%         for j=1:length(S.pleGlobals.figPath)
%             if(~isempty(S.pleGlobals.chi2s{j}))
%                 lp(fid, '\\clearpage\n');
%                 count = count + 1;
%                 targetstr = [savePath '/' S.pleGlobals.p_labels{j} '.pdf'];
%                 eval(['!ps2pdf  -dEPSCrop ' S.pleGlobals.figPath{j} '.eps ' targetstr]);
%                 
%                 captiontext = sprintf('\\textbf{Profile likelihood of parameter %s}\\\\', strrep(S.pleGlobals.p_labels{j},'_','\_'));
%                 captiontext = [captiontext 'Upper panel: The solide line indicates the profile likelihood. '];
%                 captiontext = [captiontext 'The broken line indicates the threshold to assess confidence intervals. '];
%                 captiontext = [captiontext 'The asterisk indicate the optimal parameter values. '];
%                 captiontext = [captiontext sprintf('Lower panel: The functional relations to the other parameters along the profile likelihood of %s are displayed. ', strrep(S.pleGlobals.p_labels{j},'_','\_'))];
%                 captiontext = [captiontext 'In the legend the top five parameters showing the strongest variations are given. '];
%                 captiontext = [captiontext sprintf('The calculation time was %s.', secToHMS(S.pleGlobals.timing(j)))];
%                 lpfigure(fid, 1, [S.pleGlobals.p_labels{j} '.pdf'], captiontext, sprintf('ple%i',count)); 
%                 
%                 lp(fid, '\n');
%             end
%         end
%     end
%     
%     %% Confidence Intervals
%     lp(fid, '\\clearpage\n');
%     lp(fid, '\\section{Confidence intervals for the model parameters}\n');
%     
%     N = 30;
%     ntables = 0;
%     for j=1:length(S.pleGlobals.p_labels)
%         if(S.pleGlobals.q_fit(j))
%             if(j<=length(S.pleGlobals.chi2s) && ~isempty(S.pleGlobals.chi2s{j}))
%                 ntables = ntables + 1;
%             end
%         end
%     end
%     ntables = ceil(ntables/N);
%     
%     if(ntables>1)
%         lp(fid, 'In Table \\ref{conftable1} -- \\ref{conftable%i}, %2i\\%% confidence intervals for the estimated parameter values derived by the profile likelihood \\cite{Raue:2009ec} are given.', ntables, (1-S.pleGlobals.alpha_level)*100);
%     else
%         lp(fid, 'In Table \\ref{conftable1}, %2i\\%% confidence intervals for the estimated parameter values derived by the profile likelihood \\cite{Raue:2009ec} are given.', (1-S.pleGlobals.alpha_level)*100);
%     end
%     
%     headstr = '\t\t\t & name & $\\hat\\theta$';
%     if(S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
%         headstr = [headstr '& $\\sigma^{-}_{ptw}$ & $\\sigma^{+}_{ptw}$'];
%         headstr = [headstr '& $\\sigma^{-}_{sim}$ & $\\sigma^{+}_{sim}$'];
%     else
%         headstr = [headstr '& $\\sigma^{-}$ & $\\sigma^{+}$'];
%     end
%     headstr = [headstr ' \\\\'];
%     
%     lp(fid, '\t\\begin{table}');
%     lp(fid, '\t\\dobegincenter');
%     lp(fid, '\t{\\tiny');
%     lp(fid, '\t\t\\begin{tabular}{lllllll}');
%     lp(fid, '\t\t\t\\toprule');
%     lp(fid, headstr);
%     lp(fid, '\t\t\t\\midrule');
%     
%     count = 0;
%     counttab = 1;
%     for j=1:length(S.pleGlobals.p_labels)
%         if(S.pleGlobals.q_fit(j))
%             if(j<=length(S.pleGlobals.chi2s) && ~isempty(S.pleGlobals.chi2s{j}))
%                 count = count + 1;
%                 
%                 if(mod(count,N)==0)
%                     lp(fid, '\t\t\t\\botrule');
%                     lp(fid, '\t\t\\end{tabular}}');
%                     lp(fid, '\t\t\\mycaption{Confidence intervals for the estimated parameter values derived by the profile likelihood}{conftable%i}', counttab);
%                     lp(fid, '\t\t{$\\hat\\theta$ indicates the estimated optimal parameter value.');
%                     if(S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
%                         lp(fid, '\t\t$\\sigma^{-}_{ptw}$ and $\\sigma^{+}_{ptw}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%                         lp(fid, '\t\t$\\sigma^{-}_{sim}$ and $\\sigma^{+}_{sim}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%                     elseif(S.pleGlobals.plot_point && ~S.pleGlobals.plot_simu)
%                         lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%                     elseif(~S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
%                         lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%                     end
%                     lp(fid, '}');
%                     lp(fid, '\t\\doendcenter');
%                     lp(fid, '\t\\end{table}');
%                     counttab = counttab + 1;
%                     
%                     lp(fid, '\t\\begin{table}');
%                     lp(fid, '\t\\dobegincenter');
%                     lp(fid, '\t{\\tiny');
%                     lp(fid, '\t\t\\begin{tabular}{lllllll}');
%                     lp(fid, '\t\t\t\\toprule');
%                     lp(fid, headstr);
%                     lp(fid, '\t\t\t\\midrule');
%                 end
%                 
%                 lp(fid, '\t\t\t%i & %s & %+8.3f & ', j, ...
%                     strrep(S.pleGlobals.p_labels{j},'_','\_'), S.pleGlobals.p(j));
%                 if(S.pleGlobals.plot_point)
%                     lp(fid, '%+8.3f & %+8.3f &', S.pleGlobals.conf_lb_point(j), S.pleGlobals.conf_ub_point(j));
%                 end
%                 if(S.pleGlobals.plot_simu)
%                     lp(fid, '%+8.3f & %+8.3f', S.pleGlobals.conf_lb(j), S.pleGlobals.conf_ub(j));
%                 end
%                 lp(fid, ' \\\\');
%             end
%         end
%     end
%     lp(fid, '\t\t\t\\botrule');
%     lp(fid, '\t\t\\end{tabular}}');
%     lp(fid, '\t\t\\mycaption{Confidence intervals for the estimated parameter values derived by the profile likelihood}{conftable%i}', counttab);
%     lp(fid, '\t\t{$\\hat\\theta$ indicates the estimated optimal parameter value.');
%     if(S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
%         lp(fid, '\t\t$\\sigma^{-}_{ptw}$ and $\\sigma^{+}_{ptw}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%         lp(fid, '\t\t$\\sigma^{-}_{sim}$ and $\\sigma^{+}_{sim}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%     elseif(S.pleGlobals.plot_point && ~S.pleGlobals.plot_simu)
%         lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% point-wise confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%     elseif(~S.pleGlobals.plot_point && S.pleGlobals.plot_simu)
%         lp(fid, '\t\t$\\sigma^{-}$ and $\\sigma^{+}$ indicate %i\\%% simultaneous confidence intervals.', (1-S.pleGlobals.alpha_level)*100);
%     end
%     lp(fid, '}');
%     lp(fid, '\t\\doendcenter');
%     lp(fid, '\t\\end{table}');
%     
%     %% Confidence intervals of model trajectories
% %     \subsection{Confidence intervals of the predicted model dynamics} \label{obsanalysis}
%     % TODO
% end

lp(fid, '\\bibliographystyle{plain}');
lp(fid, '\\bibliography{/Users/araue/Sites/pub/bibtex/Library}');

lp(fid, '\\end{document}');



fclose(fid);
fprintf('done\n');

%% pdflatex
fprintf('pdflatex, file %s...', fname);
cd(savePath);
eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
eval(['!bibtex ' fnamebib ' > log_bibtex.txt']);
eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
eval(['!pdflatex ' fname ' > log_pdflatex.txt']);
cd('../../..');
copyfile([savePath '/' 'presentation.pdf'], [savePath '/' sprintf('presentation_%s.pdf', datestr(now,30))])
fprintf('done\n');


function lp(varargin)
if(nargin>2)
    fprintf(varargin{1}, sprintf('%s\n', varargin{2}), varargin{3:end});
else
    fprintf(varargin{1}, sprintf('%s\n', varargin{2}));
end

function lpfigure(fid, textwidth, figpath, figcaption, figlabel)
lp(fid, '\\begin{figure}');
lp(fid, '\\begin{center}');
if(~isempty(figcaption))
    lp(fid, '\\includegraphics[height=%f\\textheight]{%s} \\caption{%s} \\label{%s}', textwidth, figpath, figcaption, figlabel);
else
    lp(fid, '\\includegraphics[height=%f\\textheight]{%s} \\label{%s}', textwidth, figpath, figlabel);
end
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

str = latex(sym(str));

if(verLessThan('symbolic', '4.9')) % maple
	str = strrep(str, '\,', ' \cdot ');

	xtmp = strrep(ar.model(jm).x, '_', '\_');
	utmp = strrep(ar.model(jm).u, '_', '\_');
	ptmp = strrep(ar.model(jm).p, '_', '\_');

	for jx = 1:length(ar.model(jm).x)
		str = strrep(str, sprintf('{\\it %s}', xtmp{jx}), sprintf('\\mathrm{[%s]}', xtmp{jx}));
	end
	for ju = 1:length(ar.model(jm).u)
		str = strrep(str, sprintf('{\\it %s}', utmp{ju}), sprintf('\\mathrm{[%s]}', utmp{ju}));
	end
	for jp = 1:length(ar.model(jm).p)
		str = strrep(str, sprintf('{\\it %s}', ptmp{jp}), sprintf('\\mathrm{%s}', ptmp{jp}));
	end
else
	str = strrep(str, '_', '\_');
	str = strrep(str, '}\_{', '\_');
	str = strrep(str, '\,', ' \cdot ');
	str = strrep(str, '}\mathrm{', '');

	xtmp = strrep(ar.model(jm).x, '_', '\_');
	utmp = strrep(ar.model(jm).u, '_', '\_');

	for jx = 1:length(ar.model(jm).x)
		str = strrep(str, sprintf('\\mathrm{%s}', xtmp{jx}), sprintf('\\mathrm{[%s]}', xtmp{jx}));
	end
	for ju = 1:length(ar.model(jm).u)
		str = strrep(str, sprintf('\\mathrm{%s}', utmp{ju}), sprintf('\\mathrm{[%s]}', utmp{ju}));
	end
end
if(iscell(str))
    for j=1:length(str)
        str{j} = [str{j} ' \nonumber'];
    end
else
    str = [str ' \nonumber'];
end


function fprintnumtab(fid, num)
fprintf(fid, '& %s ', sprintf('%g', num));



% %% Residuals
% if(isfield(ar.report, 'residualPlotPath'))
% 	lp(fid, '\\subsection{Residual plots}');
% 	copyfile([ar.report.residualPlotPath '.pdf'], ...
% 		['./' savePath '/residuals.pdf'])
% 	lp(fid, '\\begin{center}');
% 	lp(fid, '\\includegraphics[width=\\textwidth]{residuals.pdf}');
% 	lp(fid, '\\end{center}');
% end
%
%
%
