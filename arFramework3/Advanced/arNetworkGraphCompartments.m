% arNetworkGraphCompartments([m])
%
% writes network graph with respect to the defined compartments as a .dot 
% graphics file and compile to pdf 
%
% m:	model index (default: all models within the current loaded project)
%
% Note: graphviz package hase to be installed. Homepage: https://graphviz.gitlab.io
% Under windows, dot has to be excecuded by hand to translate the source file to a pdf.
%
% See also: arNetworkGraph

function arNetworkGraphCompartments(m)

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('m', 'var'))
	for j=1:length(ar.model)
		arNetworkGraphCompartments(j);
	end
	return
end

if(~ar.model(m).isReactionBased)
    fprintf('arPlotV: model %i is not reaction based, plotting omitted\n', m);
end

savePath = [arSave '/Figures/NetworkCompartments'];

if(~exist(savePath, 'dir'))
    mkdir(savePath)
end

savePath = [savePath sprintf('/%s', ar.model(m).name) '_Compartments'];

fid = fopen([savePath '.dot'], 'w');


fprintf(fid, 'digraph G {\n');

fprintf(fid, '\tgraph [overlap=false];\n');

% measured?
if(isfield(ar.model(m), 'data'))
    qu = false(size(ar.model(m).us));
    qx = false(size(ar.model(m).xs));
    if(isfield(ar.model(1).data(1), 'qu_measured'))
        for jd=1:length(ar.model(m).data)
            if(~isempty(ar.model(m).us) && ~isempty(ar.model(m).data(jd).qu_measured))
                qu = qu | ar.model(m).data(jd).qu_measured;
            end
            if(~isempty(ar.model(m).data(jd).qx_measured))
                qx = qx | ar.model(m).data(jd).qx_measured;
            end
        end
    end
else
    qu = true(size(ar.model(m).us));
    qx = true(size(ar.model(m).xs));
end

% measured inputs
if(sum(qu)>0)
    for j=find(qu)
        fprintf(fid, '\t{node [shape=diamond]; %s;}\n', ar.model(m).u{j});
    end
end

% not measured inputs
if(sum(~qu)>0)
    for j=find(~qu)
        fprintf(fid, '\t{node [style=dashed]; %s;}\n', ar.model(m).u{j});
    end
end

for jc = 1:length(ar.model(m).c)
    qc = ar.model(m).cLink == jc;
    
    fprintf(fid, '\tsubgraph cluster_%i {\n', jc);
    
    % measured species
    if(sum(qx & qc)>0)
        for jj=find(qx & qc)
            fprintf(fid, '\t\t{%s;}\n', ar.model(m).x{jj});
        end
    end
    
    % not measured species
    if(sum(~qx & qc)>0)
        for jj=find(~qx & qc)
            fprintf(fid, '\t\t{node [style=dashed]; %s ;}\n', ar.model(m).x{jj});
        end
    end
    
    % fluxes inside compartments
    for jf=1:size(ar.model(m).N,2)
        qsource = ar.model(m).N(:,jf) < 0;
        csource = ar.model(m).cLink(qsource);
        csource = unique(csource);
        qtarget = ar.model(m).N(:,jf) > 0;
        ctarget = ar.model(m).cLink(qtarget);
        ctarget = unique(ctarget);
        if(length(csource)>1)
            error('different compartments for educts');
        end
        if(isempty(csource))
            if(ctarget == jc)
                fprintf(fid, '\t\t{node [label="" shape=point]; source%i;}\n', jf);
            end
        elseif(isempty(ctarget))
            if(csource == jc)
                fprintf(fid, '\t\t{node [label="" shape=point]; target%i;}\n', jf);
            end
        end
        if((isempty(csource) && ctarget == jc) || (isempty(ctarget) && csource == jc) || (csource == jc && ctarget == jc))
            fprintf(fid, '\t\t{node [shape=box fontsize=10 width=0.1 height=0.1]; v%i ;}\n', jf);
        end
    end
    
	fprintf(fid, '\t\tlabel="%s"\n', ar.model(m).c{jc});
    
    fprintf(fid, '\t}\n');
end

% fluxes outside compartments
for jf=1:size(ar.model(m).N,2)
    qsource = ar.model(m).N(:,jf) < 0;
    csource = ar.model(m).cLink(qsource);
    csource = unique(csource);
    qtarget = ar.model(m).N(:,jf) > 0;
    ctarget = ar.model(m).cLink(qtarget);
    ctarget = unique(ctarget);
    if(length(csource)>1)
        error('different compartments for educts');
    end
    if(~isempty(csource) && ~isempty(ctarget) && csource ~= ctarget)
        fprintf(fid, '\t\t{node [shape=box fontsize=10 width=0.1 height=0.1]; v%i ;}\n', jf);
    end
end

% edges
normweight = 10;
for jj = 1:size(ar.model(m).N,2)
	source = find(ar.model(m).N(:,jj) < 0);
	target = find(ar.model(m).N(:,jj) > 0);
	
	if(~isempty(source))
		for j=1:length(source)
			fprintf(fid, '\t\t%s -> v%i [arrowhead=none weight=%i];\n', ...
				ar.model(m).x{source(j)}, jj, normweight);
		end
	else
% 		fprintf(fid, '\t\t{node [label="" shape=point]; source%i;}\n', jj);
		fprintf(fid, '\t\tsource%i -> v%i [arrowhead=none weight=%i];\n', ...
			jj, jj, normweight);
	end

	if(~isempty(target))
		for j=1:length(target)
			fprintf(fid, '\t\tv%i -> %s [weight=%i];\n', ...
				jj, ar.model(m).x{target(j)}, normweight);
		end
	else
% 		fprintf(fid, '\t\t{node [label="" shape=point]; target%i;}\n', jj);
		fprintf(fid, '\t\tv%i -> target%i [weight=%i];\n', jj, jj, normweight);
	end
end

% effect edges species
Nspecies = transpose(ar.model(m).qdvdx_nonzero);
Nspecies = Nspecies & (ar.model(m).N == 0);
for jj = 1:size(Nspecies,2)
	ij = find(Nspecies(:,jj)==1);
	for j = 1:length(ij);
		if(ar.model(m).qdvdx_negative(jj,ij(j)) == 1)
			fprintf(fid, '\t\t%s -> v%i [arrowhead=tee color=red];\n', ar.model(m).x{ij(j)}, jj);
		else
			fprintf(fid, '\t\t%s -> v%i [arrowhead=odot color=blue];\n', ar.model(m).x{ij(j)}, jj);
		end
	end
end

% effect edges inputs
Ninput = transpose(ar.model(m).qdvdu_nonzero);
for jj = 1:size(Ninput,2)
	ij = find(Ninput(:,jj)==1);
	for j = 1:length(ij);
		if(ar.model(m).qdvdu_negative(jj,ij(j)) == 1)
			fprintf(fid, '\t\t%s -> v%i [arrowhead=tee color=red];\n', ar.model(m).u{ij(j)}, jj);
		else
			fprintf(fid, '\t\t%s -> v%i [arrowhead=odot color=blue];\n', ar.model(m).u{ij(j)}, jj);
		end
	end
end

fprintf(fid, '}\n');
fclose(fid);

%% Render and save
if(ismac)
    eval(['!/usr/local/bin/dot -Tpdf ' savePath '.dot -o' savePath '.pdf']);
    % eval(['!neato -Tpdf ' savePath '.dot -o' savePath '.pdf']);
    % eval(['!fdp -Tpdf ' savePath '.dot -o' savePath '.pdf']);
elseif(isunix)
      eval(['!dot -Tpdf ' savePath '.dot -o' savePath '.pdf']);
    % eval(['!neato -Tpdf ' savePath '.dot -o' savePath '.pdf']);
    % eval(['!fdp -Tpdf ' savePath '.dot -o' savePath '.pdf']);
else
    fprintf('%s has been saved.\n',[savePath '.dot']);
    warning('Under windows, dot has to be excecuded by hand to translate the source file to a pdf.')
end