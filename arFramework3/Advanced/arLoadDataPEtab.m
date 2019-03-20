% arLoadDataPEtab(dataname, [m])
% 
% This function can be used to process data files in the format of PEtab.
% 
%   dataname  name of file.
% 
%   m         model that shall be linked to. Name or ID [length(ar.model)]
%
% In this data format, there is one single .tsv-file that contains all data 
% points. This data format shall allow easier transitions between modeling 
% tools.
% This load-data-function utilizes a different approach than the usual data.def file: 
% Firstly, an empty data struct is created which is subsequently appended to the ar struct.
% See also arCreateDataStruct arAddDataStruct
% References
%   - https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md

function arLoadDataPEtab(dataname, m)

global ar;

if(~exist('m','var') || isempty(m))
    m = length(ar.model);
end
if(exist('m','var') && ischar(m))
    for jm=1:length(ar.model)
        if(strcmp(m, ar.model(jm).name))
            m = jm;
        end
    end
    if(ischar(m))
        error('Model %s was not found', m);
    end
end


%% Read in tsv file
T = tdfread(dataname); % all data of the model
fns = fieldnames(T);
for i = 1:length(fns)
    if ischar(T.(fns{i}))
        T.(fns{i}) = regexprep(string(T.(fns{i})),' ','');
    end
end
T = struct2table(T);
controlstr = {};
for i =1:height(T) % find unique combination to put in one data struct afterwards
    controlstr{i} = string(strcat(char(T.simulationConditionId(i))));%,char(T.observableParameters(i)),char(T.noiseParameters(i))));
end
controlstr = string(controlstr);
[uniSim,~,iCSim] = unique(controlstr);

%% Use condition specific experiments and distribute over data struct
for iSim = 1:length(uniSim)
    Sd2d = struct();
    args= {};
    % extract important info for data struct from .tsv file
    Tsub = T(iCSim == iSim,:);
    [uniObs,~,iCObs] = unique(cellstr(Tsub.observableId));
    [uniTimes,~,iTExp] = unique(Tsub.time);
    uniObs = regexprep(uniObs,' ','');
    Sd2d.tExp = uniTimes;
    Sd2d.y = uniObs';
    Sd2d.yNames = uniObs';
    Sd2d.yExpRaw = nan(length(uniTimes),length(uniObs));
    Sd2d.yExpStd = nan(length(uniTimes),length(uniObs));
    Sd2d.yExpStdRaw = nan(length(uniTimes),length(uniObs));
    for it = 1:length(uniTimes)
        for iobs = 1:length(uniObs)
            Sd2d.yExpRaw(it,iobs) = Tsub.measurement(it == iTExp & iobs == iCObs);
            Sd2d.yExp(it,iobs) = strcmp(Tsub.observableTransformation(it == iTExp & iobs == iCObs),'log10')*log10(Tsub.measurement(it == iTExp & iobs == iCObs)) + strcmp(Tsub.observableTransformation(it == iTExp & iobs == iCObs),'lin')*Tsub.measurement(it == iTExp & iobs == iCObs);
            if isnumeric(Tsub.noiseParameters(it == iTExp & iobs == iCObs))
                Sd2d.yExpStdRaw(it,iobs) = Tsub.noiseParameters(it == iTExp & iobs == iCObs);
                Sd2d.yExpStd(it,iobs) = strcmp(Tsub.observableTransformation(it == iTExp & iobs == iCObs),'log10')*log10(Tsub.noiseParameters(it == iTExp & iobs == iCObs)) + strcmp(Tsub.observableTransformation(it == iTExp & iobs == iCObs),'lin')*Tsub.noiseParameters(it == iTExp & iobs == iCObs);
            end
        end
    end
    % prepare info for creating data struct
    fns = fieldnames(Sd2d);
    for i = 1:length(fns)
        args(end+1) = fns(i);
        args{end+1} = Sd2d.(fns{i});
    end
    args{end+1} = 'doseresponse'; args{end+1} = 0;
    if rem(length(args),2)~=0
        error('arguments args has to be provided in pairs.')
    end
    
    % get correct substitutions of obs and error parameters
    pold = {}; fp = {};
    for i = 1:length(Sd2d.y)
        newpold = ar.model.fp(contains(ar.model.p,'noiseParameter') & contains(ar.model.p,Sd2d.y{i}));
        newPars = Tsub.noiseParameters(contains(Tsub.observableId,Sd2d.y{i}));
        if ~(isempty(newpold) || isnumeric(newPars) || isempty(newPars))
            newPars = strsplit(str2mat(unique(newPars)),';');
            pold = [pold(:);newpold];fp = [fp(:);newPars];    
        end
        newpold = ar.model.fp(contains(ar.model.p,'observableParameter') & contains(ar.model.p,Sd2d.y{i}));
        newPars = Tsub.observableParameters(contains(Tsub.observableId,Sd2d.y{i}));
        if ~(isempty(newpold) || isnumeric(newPars) || isempty(newPars))
            newPars = strsplit(str2mat(unique(newPars)),';');
            pold = [pold(:);newpold];fp = [fp(:);newPars];
        end
        if length(pold) ~= length(fp)
            error('Length of parameter names must be the same!')
        end
    end
    
    D = arCreateDataStruct(m,pold,fp,args{:});
    arAddDataStruct(D,m)
end

end