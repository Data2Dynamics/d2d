% write .dot graphics file and compile to pdf
%
% arNetworkGraph(m)
% m:	model index

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

savePath = [savePath sprintf('/%s', ar.model(m).name)];

fid = fopen([savePath '.dot'], 'w');


fprintf(fid, 'digraph G {\n');

fprintf(fid, '\tgraph [overlap=false];\n');

% measured?
qu = false(size(ar.model(m).us));
qx = false(size(ar.model(m).xs));
if(isfield(ar.model(m), 'data') && isfield(ar.model(1).data(1), 'qu_measured'))
    for jd=1:length(ar.model(m).data)
        if(~isempty(ar.model(m).us) && ~isempty(ar.model(m).data(jd).qu_measured))
            qu = qu | ar.model(m).data(jd).qu_measured;
        end
        if(~isempty(ar.model(m).data(jd).qx_measured))
            qx = qx | ar.model(m).data(jd).qx_measured;
        end
    end
end

% measured inputs
if(sum(qu)>0)
    fprintf(fid, '\t{node [shape=diamond];');
    for j=find(qu)
        fprintf(fid, '%s ', ar.model(m).u{j});
    end
    fprintf(fid, ';}\n');
end

% not measured inputs
if(sum(~qu)>0)
    fprintf(fid, '\t{node [style=dashed];');
    for j=find(~qu)
        fprintf(fid, '%s ', ar.model(m).u{j});
    end
    fprintf(fid, ';}\n');
end

for jc = 1:length(ar.model(m).c)
    qc = ar.model(m).cLink == jc;
    
    fprintf(fid, '\tsubgraph cluster_%i {\n', jc);
    
    % measured species
    if(sum(qx & qc)>0)
        fprintf(fid, '\t\t{');
        for jj=find(qx & qc)
            fprintf(fid, '%s ', ar.model(m).x{jj});
        end
        fprintf(fid, ';}\n');
    end
    
    % not measured species
    if(sum(~qx & qc)>0)
        fprintf(fid, '\t\t{node [style=dashed]; ');
        for jj=find(~qx & qc)
            fprintf(fid, '%s ', ar.model(m).x{jj});
        end
        fprintf(fid, ';}\n');
    end
    
    % fluxes
    fprintf(fid, '\t\t{node [shape=box fontsize=10 width=0.1 height=0.1]; ');
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
            csource = ctarget;
            if(csource == jc)
                fprintf(fid, '{node [label="" shape=point]; source%i;} ', jf);
            end
        end
        if(isempty(ctarget))
            if(csource == jc)
                fprintf(fid, '{node [label="" shape=point]; target%i;}', jf);
            end
        end
        if(csource == jc)
            fprintf(fid, 'v%i ', jf);
        end
    end
    fprintf(fid, ';}\n');
    
	fprintf(fid, '\t\tlabel="%s"\n', ar.model(m).c{jc});
end
for jc = 1:length(ar.model(m).c)
	fprintf(fid, '\t}\n');
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
		fprintf(fid, '\t\t{node [label="" shape=point]; source%i;}\n', jj);
		fprintf(fid, '\t\tsource%i -> v%i [arrowhead=none weight=%i];\n', ...
			jj, jj, normweight);
	end

	if(~isempty(target))
		for j=1:length(target)
			fprintf(fid, '\t\tv%i -> %s [weight=%i];\n', ...
				jj, ar.model(m).x{target(j)}, normweight);
		end
	else
		fprintf(fid, '\t\t{node [label="" shape=point]; target%i;}\n', jj);
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
eval(['!dot -Tpdf ' savePath '.dot -o' savePath '.pdf']);
% eval(['!neato -Tpdf ' savePath '.dot -o' savePath '.pdf']);
% eval(['!fdp -Tpdf ' savePath '.dot -o' savePath '.pdf']);

