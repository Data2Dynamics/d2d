% write .dot graphics file and compile to pdf
%
% arNetworkGraph(m)
% m:	model index

function arNetworkGraph(m)

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('m', 'var'))
	for j=1:length(ar.model)
		arNetworkGraph(j);
	end
	return
end

if(~ar.model(m).isReactionBased)
    fprintf('arPlotV: model %i is not reaction based, plotting omitted\n', m);
end

if(isempty(ar.model(m).x))
    return;
end

savePath = [arSave '/Figures/Network'];

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
    fprintf(fid, '\t{node [shape=diamond style=bold style=filled fillcolor=gray90];');
    for j=find(qu)
        fprintf(fid, '%s ', ar.model(m).u{j});
    end
    fprintf(fid, ';}\n');
end

% not measured inputs
if(sum(~qu)>0)
    fprintf(fid, '\t{node [shape=diamond style=filled fillcolor=gray90];');
    for j=find(~qu)
        fprintf(fid, '"%s" ', ar.model(m).u{j});
    end
    fprintf(fid, ';}\n');
end

% measured species
if(sum(qx)>0)
	fprintf(fid, '\t\t{node [style=bold];');
	for jj=find(qx)
		fprintf(fid, '"%s" ', ar.model(m).xNames{jj});
	end
	fprintf(fid, ';}\n');
end

% not measured species
if(sum(~qx)>0)
    fprintf(fid, '\t\t{');
	for jj=find(~qx)
		fprintf(fid, '"%s" ', ar.model(m).xNames{jj});
	end
	fprintf(fid, ';}\n');
end

% fluxes
fprintf(fid, '\t\t{node [shape=box fontsize=10 width=0.1 height=0.1]; ');
for jf=1:size(ar.model(m).N,2)
	fprintf(fid, 'v%i ', jf);
end
fprintf(fid, ';}\n');

% edges
for jj = 1:size(ar.model(m).N,2)
	source = find(ar.model(m).N(:,jj) < 0);
	target = find(ar.model(m).N(:,jj) > 0);
	
	if(~isempty(source))
		for j=1:length(source)
			fprintf(fid, '\t\t"%s" -> v%i [arrowhead=none weight=10];\n', ...
				ar.model(m).xNames{source(j)}, jj);
		end
	else
		fprintf(fid, '\t\t{node [label="" shape=point]; source%i;}\n', jj);
		fprintf(fid, '\t\tsource%i -> v%i [arrowhead=none weight=10];\n', jj, jj);
	end

	if(~isempty(target))
		for j=1:length(target)
			fprintf(fid, '\t\tv%i -> "%s" [weight=10];\n', ...
				jj, ar.model(m).xNames{target(j)});
		end
	else
		fprintf(fid, '\t\t{node [label="" shape=point]; target%i;}\n', jj);
		fprintf(fid, '\t\tv%i -> target%i [weight=10];\n', jj, jj);
	end
end

% effect edges species
Nspecies = transpose(ar.model(m).qdvdx_nonzero);
for c=2:length(ar.model(m).condition)
    Nspecies = Nspecies | transpose(ar.model(m).qdvdx_nonzero);
end
Nspecies = Nspecies & (ar.model(m).N == 0);
for jj = 1:size(Nspecies,2)
	ij = find(Nspecies(:,jj)==1);
	for j = 1:length(ij);
		if(ar.model(m).qdvdx_negative(jj,ij(j)) == 1)
			fprintf(fid, '\t\t"%s" -> v%i [arrowhead=tee color=red weight=1];\n', ar.model(m).xNames{ij(j)}, jj);
		else
			fprintf(fid, '\t\t"%s" -> v%i [arrowhead=odot color=blue weight=5];\n', ar.model(m).xNames{ij(j)}, jj);
		end
	end
end

% effect edges inputs
Ninput = transpose(ar.model(m).qdvdu_nonzero);
for c=2:length(ar.model(m).condition)
    Ninput = Ninput | transpose(ar.model(m).qdvdu_nonzero);
end
for jj = 1:size(Ninput,2)
	ij = find(Ninput(:,jj)==1);
	for j = 1:length(ij);
		if(ar.model(m).qdvdu_negative(jj,ij(j)) == 1)
			fprintf(fid, '\t\t%s -> v%i [arrowhead=tee color=red weight=1];\n', ar.model(m).u{ij(j)}, jj);
		else
			fprintf(fid, '\t\t%s -> v%i [arrowhead=odot color=blue weight=5];\n', ar.model(m).u{ij(j)}, jj);
		end
	end
end

fprintf(fid, '}\n');
fclose(fid);

%% Render and save
if(~ispc)
    eval(['!dot -Tpdf ' savePath '.dot -o' savePath '.pdf']);
    % eval(['!neato -Tpdf ' savePath '.dot -o' savePath '.pdf']);
    % eval(['!fdp -Tpdf ' savePath '.dot -o' savePath '.pdf']);
end
