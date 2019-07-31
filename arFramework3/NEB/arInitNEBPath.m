function arInitNEBPath
% Initializes NEB path nodes and their parameters (not initial, e.g. direct path)

global ar

ar.merger.neb.ps_index = [];
% get NEB parameter indicies

ar.qFit(ar.qFit==1) = 0;

% for basic model step 0
for j = 1:length(ar.merger.neb.qFitlable)
    	ar.merger.neb.ps_index(1,j) = arGetParsIndex( [char(ar.merger.neb.qFitlable(j))] );
        ar.qFit(ar.merger.neb.ps_index(1,j)) = 0;
end

% for NEB nodes
for i = 1:ar.merger.neb.steps
    for j = 1:length(ar.merger.neb.qFitlable)
    	ar.merger.neb.ps_index(i+1,j) = arGetParsIndex( [char(ar.merger.neb.qFitlable(j)) '_NEB_' sprintf('%03d', i)] );
        ar.qFit(ar.merger.neb.ps_index(i+1,j)) = 1;
        
        ar.lb(ar.merger.neb.ps_index(i+1,j)) = ar.lb(ar.merger.neb.ps_index(1,j));
        ar.up(ar.merger.neb.ps_index(i+1,j)) = ar.ub(ar.merger.neb.ps_index(1,j));
        
    end
end

% for basic model step at end
for j = 1:length(ar.merger.neb.qFitlable)
    	ar.merger.neb.ps_index(ar.merger.neb.steps+2,j) = arGetParsIndex( [char(ar.merger.neb.qFitlable(j)) '_NEB_end'] );
        ar.qFit(ar.merger.neb.ps_index(ar.merger.neb.steps+2,j)) = 0;
end
