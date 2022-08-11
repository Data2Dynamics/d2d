function [algebraicEqns,ODEs,stateNames,stateValues,ICStateEqns,...
    parNames,parValues,ParIndex,inputNames,inputEqns,outputNames,outputEqns,...
    outputEqnsOriginal] = process_model_file(tree)

% remove semicolons and spaces from input string tree
tree=strrep(tree,';',''); tree=strrep(tree,' ','');
% remove tabs from input string tree
tree=strrep(tree,'\t','');
% remove comments from input string tree
tree=RemoveCommentsFromCellStringArray(tree);
% find indices for different sections in the input file
Index=strfind(tree,'!');
MasterIndex=find(~cellfun(@isempty,Index'));

if not(isequal(numel(MasterIndex),7))
    msg={'Input file is missing a section \n';
    'Order of the sections in the input file is\n';
    'section 1: Algebraic relations such as reaction rates\n';
    'section 2: ODEs\n';
    'section 3: Inputs\n';
    'section 4: Outputs\n';
    'section 5: Parameters\n';
    'section 6: States\n';
    'section 7: States/Parameters to be analysed for observability\n'};
    error('u:stuffed:it',strjoin(msg));
end

% Read algebraic equations
nv=MasterIndex(2)-MasterIndex(1)-1; algebraicEqns=cell(nv,1);
algebraicNames=cell(nv,1);
for j=1:nv
    d=strsplit(tree{MasterIndex(1)+j},'=');
    algebraicNames{j}=strtrim(d{1});
    algebraicEqns{j}=strtrim(d{2});
end

% Remove comments from equations
% algebraicEqns = RemoveCommentsFromCellStringArray(algebraicEqns);

% Find states and parameters, so we know nx and np
nx=MasterIndex(7)-MasterIndex(6)-1; stateNames=cell(nx,1); ICStateEqns=cell(nx,1);
icNames=cell(nx,1); stateValues=NaN(nx,1);
for j=1:nx
    d=strtrim(strsplit(tree{MasterIndex(6)+j},'='));
    stateNames{j}=strtrim(d{1});
    icNames{j}=['IC',num2str(j)];
    ICStateEqns{j}=['IC',num2str(j)];
    if isequal(numel(d),2)
        if ~isempty(d{2})
            value=str2double(d{2});
            if isnan(value)
                ICStateEqns{j}=d{2};
            else
                stateValues(j)=value;
            end
        end
    end
end

np=MasterIndex(6)-MasterIndex(5)-1; parNames=cell(np,1); parValues=NaN(np,1);
for j=1:np
    d=strtrim(strsplit(tree{MasterIndex(5)+j},'='));
    parNames{j}=d{1};
    if isequal(numel(d),2)
        if ~isempty(d{2}), parValues(j)=str2double(d{2});
        end
    end
end

% Find outputs available for measurements
ny=MasterIndex(5)-MasterIndex(4)-1; outputNames=cell(ny,1); outputEqns=cell(ny,1);
for j=1:ny
    d=strtrim(strsplit(tree{MasterIndex(4)+j},'='));
    outputNames{j}=d{1};
    outputEqns{j}=d{2};
end
% outputEqnsRaw=outputEqns;

% Find inputs
nu=MasterIndex(4)-MasterIndex(3)-1;
if nu~=0
    inputNames=cell(nu,1); inputEqns=cell(nu,1);
    for j=1:nu
        d=strtrim(strsplit(tree{MasterIndex(3)+j},'='));
        inputNames{j}=d{1};
        if ~isempty(d{2})
            inputEqns{j}=d{2};
        else
            inputEqns{j}=num2str(rand);
        end
%         if isequal(numel(d),2)
%             if ~isempty(d{2}), inputEqns(j)=str2double(d{2});
%             end
%         end
    end
    inputEqns=strcat(inputEqns,';');
else
    inputNames={};
    inputEqns={};
end

% Find parameters (including ICs) to be analysed from the final section of the input file
varNamesOld=cat(1,parNames,stateNames,inputNames,icNames,algebraicNames);
varNamesAug=cat(1,parNames,stateNames);
parAnalysed=strtrim(tree(MasterIndex(7)+1:end));
if isempty(parAnalysed)
    parAnalysed=varNamesAug; % analyse ALL variables in case nothing is in this section of the input file
end
[~,ia,~]=intersect(varNamesAug,parAnalysed);
ParIndex=sort(ia)';

% Next, find the ODEs in string-format from the second section
ODEs=cell(nx,1);
for j=1:nx
    d=strsplit(tree{MasterIndex(2)+j},'=');
    ODEs{j}=strtrim(d{2});
end

% Remove comments from ODEs
% ODEs = RemoveCommentsFromCellStringArray(ODEs);

% Substitute standard symbols x and p for states and parameters in the ODEs
% First make ONE sorted list of the stateNames and parNames and start the
% replacement with the longest name. Otherwise errors occur, since a short name
% can be part of a large name!
StringLengthVar=cellfun(@length,varNamesOld);
[~,IndexVar]=sort(StringLengthVar,'descend');

stateNamesNew=cell(nx,1); for i=1:nx, stateNamesNew{i}=['x(',num2str(i),')']; end
icNamesNew=cell(nx,1); for i=1:nx, icNamesNew{i}=['p(',num2str(np+i),')']; end
parNamesNew=cell(np,1); for i=1:np, parNamesNew{i}=['p(',num2str(i),')']; end
inputNamesNew=cell(nu,1); for i=1:nu, inputNamesNew{i}=['u(',num2str(i),')']; end
algebraicNamesNew=cell(nv,1); for i=1:nv, algebraicNamesNew{i}=['v(',num2str(i),')']; end

varNamesNew=cat(1,parNamesNew,stateNamesNew,inputNamesNew,icNamesNew,algebraicNamesNew);
nvar=numel(varNamesNew);
% now perform the substitution in all ODEs, from large to small strings
for j=1:nvar
    ind=IndexVar(j);
    ODEs=strrep(ODEs,varNamesOld{ind},varNamesNew{ind});
end
% add a semi colon if not present
%     ind=(~endsWith(ODEs,';'));
ODEs=strcat(ODEs,';');

% Do the same for the input equations
for j=1:nvar
    ind=IndexVar(j);
    inputEqns=strrep(inputEqns,varNamesOld{ind},varNamesNew{ind});
end
% add a semi colon if not present
%     ind=(~endsWith(ODEs,';'));
inputEqns=strcat(inputEqns,';');

% Do the same for the algebraic equations
% algebraicEqns=RemoveCommentsFromCellStringArray(algebraicEqns);
for j=1:nvar
    ind=IndexVar(j);
    algebraicEqns=strrep(algebraicEqns,varNamesOld{ind},varNamesNew{ind});
end
% ind=(~endsWith(algebraicEqns,';'));
algebraicEqns=strcat(algebraicEqns,';');

% Do the same for the outputEqns
% Remove comments
% outputEqns = RemoveCommentsFromCellStringArray(outputEqns);
outputEqnsOriginal=outputEqns;
for j=1:nvar
    ind=IndexVar(j);
    outputEqns=strrep(outputEqns,varNamesOld{ind},varNamesNew{ind});
end
%     ind=(~endsWith(outputEqns,';'));
outputEqns=strcat(outputEqns,';');

% And also for the Initial Conditions
% ICStateEqns=RemoveCommentsFromCellStringArray(ICStateEqns);
for j=1:nvar
    ind=IndexVar(j);
    ICStateEqns=strrep(ICStateEqns,varNamesOld{ind},varNamesNew{ind});
end
%     ind=(~endsWith(ICStateEqns,';'));
ICStateEqns=strcat(ICStateEqns,';');
