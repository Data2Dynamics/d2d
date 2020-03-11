% Tcond = arLoadCondPEtab(expcondfilename, Tarr)
%
% This function can be used to process experimental condition files in the format of PEtab.
%
%   expcondfilename    name of file.
%   Tarr               cell array of meas and obs tables
%
% In this data format, there is one .tsv-file that contains the info about
% the experimental conditions. This should be called before arCompile.
% This data format shall allow easier transitions between modeling
% tools.
%
% See also arLoadDataPEtab
%
% References
%   - https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md
%

function Tcond = arLoadCondPEtab(expcondfilename, Tarr)

global ar;

if ~contains(expcondfilename,'.tsv')
    if ~contains(expcondfilename,'.')
        expcondfilename = [expcondfilename '.tsv'];
    else
        error('this file type is not supported!')
    end
end

%% Read in tsv file
T = tdfread(expcondfilename); % all data of the model
fns = fieldnames(T);
for i = 1:length(fns)
    if ischar(T.(fns{i}))
        T.(fns{i}) = regexprep(string(T.(fns{i})),' ','');
    end
end

if ~isfield(T,'conditionId')
    T.conditionId = T.conditionID;
end

% load pre-equilibration condition
m = 1;
Tobs = Tarr{m,1};
Tmeas = Tarr{m,2};
noDataPreEqConds = setdiff(cellstr(T.conditionId), {ar.model.data.name});
for iPreEq = 1:numel(noDataPreEqConds)
    condId = find(T.conditionId == noDataPreEqConds{iPreEq});
    
    Ttable = struct2table(T);
    pold = Ttable.Properties.VariableNames(2:end);
    fp = cellfun(@num2str,table2cell(Ttable(condId, 2:end)),...
        'UniformOutput',false);
        
    args = cell(0);
    Sd2d.name = char(Ttable.conditionId(condId));
    fns2 = fieldnames(Sd2d);
    for i = 1:length(fns2)
        args(end+1) = fns2(i);
        args{end+1} = Sd2d.(fns2{i});
    end

    D = arCreateDataStruct(m,pold,fp, args{:});
    arAddDataStruct(D,m)
end

for m = 1:length(ar.model)
    for i = 1:length(T.conditionId)
        for j = 1:length(ar.model(m).data)
            if strcmp(ar.model(m).data(j).name,T.conditionId(i))
                for k = 2:(length(fns))
                    % Check if parameter is replaced by condition
                    if sum(contains(ar.model(m).data(j).fp,fns{k}))>0    % changed p to fp to catch cases in which initial value was renamed from a0 to init_A_state
                        ar.model(m).data(j).fp{ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k}))} = ...
                            num2str(T.(fns{k})(i));
                        
%                         if isnumeric(T.(fns{k})) % if condition parameter is replaced by a parameter instead of number
%                             ar.model(m).data(j).fp{ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k}))} = ...
%                                 num2str(T.(fns{k})(i));
%                         else
%                             ar.model(m).data(j).fp{ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k}))} = ...
%                                 T.(fns{k}){i};
%                         end
                        % Check if initial state is set by condition
                    elseif any(strcmp(ar.model(m).xNames,fns{k}))
                        InitialsSet = strcmp(ar.model(m).xNames,fns{k});
                        if sum(InitialsSet) > 1
                            error('Problem in finding right initial state set in condition table!')
                        end
                        initialVariable = ar.model(m).px0{InitialsSet};
                        InitialDataVariableIndex = strcmp(initialVariable, ar.model(m).data(j).p);
                        ar.model(m).data(j).fp{InitialDataVariableIndex} = num2str(T.(fns{k})(i));
                        
                     % Check if compartment parameter is set condition
                     % specific
                    elseif any(strcmp(ar.model(m).c,fns{k}))
                        CompartmentSet = strcmp(ar.model(m).c,fns{k});
                        CompartmentVariable = ar.model(m).pc{CompartmentSet};
                        CompartmentDataVariableIndex = strcmp(CompartmentVariable, ar.model(m).data(j).p);
                        ar.model(m).data(j).fp{CompartmentDataVariableIndex} = num2str(T.(fns{k})(i));
                    end
                end
            end
        end
    end
end

Tcond = struct2table(T);
end