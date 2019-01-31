% l1_PrintDynPars([m],[cell_string],[omit])
% 
% This function prints possible transformations for replacements for
% cell-specific fold-factors to the matlab cmd line,
% 
% m             [1] is default
%               Model where parameters have to be substituted
% cell_string   ['is_other'] is default
%               identifier of e.g. different cell-lines
% omit          [{}] is default    
% 
% Substitute all parameters that are fitted to L1-conform parameters
% Substitution also holds for parameter transformations formulated in the
% CONDITIONS section of the .def file
% 
% To obtain reltos of dynamic parameters between e.g. two cell lines, the
% output of l1_PrintDynParscan can be parsed into the CONDITIONS section of
% the model-def file.

function l1_PrintDynPars(m,cell_string,omit)
global ar
if(~exist('m','var') || isempty(m))
    warning('no model specified, will take m=1')
    m = 1;
end
if(~exist('cell_string','var') || isempty(cell_string))
    warning('no cell string specified, will take is_other to identify changes')
    cell_string = 'is_other';
end
if(~exist('omit','var') || isempty(omit))
    omit = {};
end
warning('off','all')
L1_pars = ar.pLabel(ar.qDynamic==1 & ar.qFit==1);
del_L1 = NaN(1,length(L1_pars));
for i = 1:length(L1_pars)
    is_modelPar = find(strcmp(ar.model(m).p,L1_pars{i}));
    if(~isempty(is_modelPar) && any(strcmp(ar.model(m).p{is_modelPar},omit)))
        del_L1(i) = i;  
    end
end
L1_pars(~isnan(del_L1)) = [];


[found_SS,found_cells] = setdiff(ar.model(m).p,ar.model(m).fp);
del_SS = [];
for iSS = 1:length(found_cells)
   tmp_regexp = fprintf('^(%s\)$',ar.model(m).p{found_cells(iSS)});
   if(contains(ar.model(m).fp{found_cells(iSS)},['relto_' ar.model(m).p{found_cells(iSS)}]) || ~isempty(regexp(ar.model.fp{found_cells(iSS)},'^(\d+\)$','ONCE')) || ~isempty(regexp(ar.model.fp{found_cells(iSS)},tmp_regexp,'ONCE')))
       del_SS(end+1) = iSS;
   end
end

L1_pars(ismember(L1_pars,found_SS(del_SS))) = [];
L1_pars_new = L1_pars;
found_cells(del_SS) = [];

saved_L1 = cell(length(found_cells),2);
for i = 1:length(L1_pars)
    new_var = sprintf('( %s *(1+%s * (relto_%s -1 )))', L1_pars{i}, cell_string, L1_pars{i});
    L1_pars_new{i} = strrep(L1_pars_new{i},L1_pars{i}, new_var);
    is_modelPar = find(strcmp(ar.model(m).p,L1_pars{i}));
    if(~isempty(is_modelPar))
        fprintf('%s  "%s" \n',L1_pars{i},L1_pars_new{i});
        saved_L1{i,1} = L1_pars{i};
        saved_L1{i,2} = L1_pars_new{i};
    end
end


for i = 1:length(found_cells)
    L1_tmp = sym(ar.model(m).fp{found_cells(i)});
    L1_tmp = char(subs(L1_tmp,L1_pars,L1_pars_new)); 
    fprintf('%s  "%s" \n',ar.model.p{found_cells(i)},L1_tmp); 
    saved_L1{i,1} = ar.model.p{found_cells(i)};
    saved_L1{i,2} = L1_tmp;
end
