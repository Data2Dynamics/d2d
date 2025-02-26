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
    T.conditionId = T.conditionID; % old version of arExportPEtab used 'conditionID'
end

% load pre-equilibration condition
m = 1;
Tobs = Tarr{m,1};
Tmeas = Tarr{m,2};
noDataPreEqConds = setdiff(cellstr(T.conditionId), cellstr(Tobs.simulationConditionId));
for iPreEq = 1:numel(noDataPreEqConds)
    condId = find(T.conditionId == noDataPreEqConds{iPreEq});
    
    Ttable = struct2table(T);
    idx = Ttable.Properties.VariableNames(~strcmpi(Ttable.Properties.VariableNames,'conditionId') & ~strcmpi(Ttable.Properties.VariableNames,'conditionName'));
    pold = Ttable.Properties.VariableNames(idx);
    fp = cellfun(@num2str,table2cell(Ttable(condId, idx)),...
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
                    %% NN: here strcmp is not sufficient. fp has to evaluated by arSym.
                    % _0002 has the problem that fp='(ao)' but p='a0'
                    if any(ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k})))    % changed p to fp to catch cases in which initial value was renamed from a0 to init_A_state
                    
                        % condition parameter is replaced by a number
                        if isnumeric(T.(fns{k})(i))
                            if sum(ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k})))==length(T.(fns{k})(i)) % if not ismember, can't replace
                                ar.model(m).data(j).fp{ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k}))} = ...
                                    num2str(T.(fns{k})(i));
                            end
                        % condition parameter is replaced by string of a number
                        % has to be treated separatly for technical reasons
                        elseif ~isnan(str2double(T.(fns{k})(i)))
                            % if any(ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k}))) % if not ismember, can't replace
                            % using 'str2sym' instead of arSym fixes an error
                            % when fp contains an equation (e.g. param_A*param_B)
                            % instead of simply one parameter param_A
                            % maybe this should be included in arSym directly?!
                            % if any(ismember(cellfun(@str2sym, ar.model(m).data(j).fp), arSym(fns{k})))
                            %     ar.model(m).data(j).fp{ismember(cellfun(@str2sym, ar.model(m).data(j).fp), arSym(fns{k}))} = ...
                            %        char(T.(fns{k})(i));
                            % end
                            ar.model(m).data(j).fp{ismember(cellfun(@str2sym, ar.model(m).data(j).fp), arSym(fns{k}))} = ...
                                char(T.(fns{k})(i));
                        
                        % condition parameter is replaced by an expression
                        else
                            % Check if parameter has already been replaced with the correct condition specific par
                            if sum(ismember(arSym(ar.model(m).data(j).fp),char(T.(fns{k}){i}))) == 0
                                % here there is an error sometimes
                                ar.model(m).data(j).fp{ismember(arSym(ar.model(m).data(j).fp), arSym(fns{k}))} = ...
                                    char(T.(fns{k}){i});
                            end
                        end
                    
                    % Check if initial state is set by condition (via stateID in PEtab)
                    elseif any(strcmp(ar.model(m).x, fns{k})) || any(strcmp(ar.model(m).x, [fns{k} '_state']))
                        if length(fns{k}) == 1
                            % if stateID in PEtab has only one character, d2d adds the suffix '_state'
                            InitialsSet = strcmp(ar.model(m).x, [fns{k}, '_state']);
                        else
                            InitialsSet = strcmp(ar.model(m).x, fns{k});
                        end
                        if sum(InitialsSet) > 1
                            error('Problem in finding right initial state set in condition table!')
                        end
                        initialVariable = ar.model(m).px0{InitialsSet};
                        InitialDataVariableIndex = strcmp(initialVariable, ar.model(m).data(j).p);
                        if isnumeric(T.(fns{k})(i))
                            ar.model(m).data(j).fp{InitialDataVariableIndex} = num2str(T.(fns{k})(i));
                        else
                            ar.model(m).data(j).fp{InitialDataVariableIndex} = char(T.(fns{k})(i));
                        end
                    
                    % initial state is set by condition (via init_StateID parameter)
                    % catch that inital parameter in sbml file can be named differently from the d2d convention
                    % note: this is not the recommended way to set initial states in PEtab (use stateID instead)
                    % elseif any(strcmp(ar.model(m).px0, fns{k})) || any(strcmp(ar.model(m).px0, [fns{k} '_state']))
                    elseif any(strcmp(ar.model(m).px0, fns{k})) || any(strcmp(ar.model(m).px0, [fns{k} '_state'])) || any(strcmp(ar.model(m).x, fns{k}))
                        if length(fns{k}) == 6
                            % all px0 have form ['init_' stateID]. If stateID in PEtab has only one character (px0 has 6 characters),
                            % d2d adds the suffix '_state'. E.g. the px0 'init_A' is replaced by 'init_A_state'
                            InitialsSet = strcmp(ar.model(m).px0, [fns{k}, '_state']);
                        else
                            InitialsSet = strcmp(ar.model(m).px0, fns{k});
                        end
                        if sum(InitialsSet) > 1
                            error('Problem in finding right initial state set in condition table!')
                        end
                        initialVariable = ar.model(m).px0{InitialsSet};
                        InitialDataVariableIndex = strcmp(initialVariable, ar.model(m).data(j).p);
                        if isnumeric(T.(fns{k})(i))
                            ar.model(m).data(j).fp{InitialDataVariableIndex} = num2str(T.(fns{k})(i));
                        else
                            ar.model(m).data(j).fp{InitialDataVariableIndex} = char(T.(fns{k})(i));
                        end
                        
                    % Check if compartment parameter is set condition-specific
                    elseif any(strcmp(ar.model(m).c,fns{k}))
                        CompartmentSet = strcmp(ar.model(m).c,fns{k});
                        
                        % Check if compartment volume parameter exists
                        % already, if not create it for all data sets
                        if isempty(str2num(ar.model(m).pc{CompartmentSet}))
                            CompartmentVariable = ar.model(m).pc{CompartmentSet};
                        else
                            CompartmentVariable = ['vol_' cell2mat(fns(k))];
                            StandardCompartmentSize = str2num(ar.model(m).pc{CompartmentSet});
                            ar.model(m).pc{CompartmentSet} = CompartmentVariable;
                            ar.model(m).px(end+1) = {CompartmentVariable};
                            for jData = 1:length(ar.model(m).data)
                                ar.model(m).data(jData).p(end+1) = {['vol_' cell2mat(fns(k))]};
                                ar.model(m).data(jData).fp(end+1) = {StandardCompartmentSize};
                            end
                        end
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