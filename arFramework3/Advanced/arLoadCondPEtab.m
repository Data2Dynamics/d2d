% arLoadCondPEtab(expcondfilename)
%
% This function can be used to process experimental condition files in the format of PEtab.
%
%   expcondfilename    name of file.
%
% In this data format, there is one .tsv-file that contains the info about
% the experimental conditions. This should be called before arCompile.
% This data format shall allow easier transitions between modeling
% tools.

% Example:
% arImportSBML(Model1.xml)
% arLoadModel(Model1.def)
% arLoadDataPEtab(measurements.tsv)
% arLoadCondPEtab(expCond.tsv)
% arCompileAll
% arLoadParsPEtab(parameters.tsv)
% 
% See also arLoadDataPEtab
%
% References
%   - https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md
%

function arLoadCondPEtab(expcondfilename)

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

for m = 1:length(ar.model)
    for j = 1:length(ar.model(m).data) 
        for i = 1:length(T.conditionId)
            if strcmp(ar.model(m).data(j).name,T.conditionId(i))
                for k = 3:(length(fns))
                    ar.model(m).data(j).fp{strcmp(ar.model(m).data(j).p,fns{k})} = num2str(T.(fns{k})(i));
                end
            end
        end
    end
end

end