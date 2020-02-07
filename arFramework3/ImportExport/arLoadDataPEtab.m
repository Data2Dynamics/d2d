% arLoadDataPEtab(datafilename, [m])
%
% This function can be used to process data files in the format of PEtab.
%
%   datafilename    name of file.
%
%   m               model that shall be linked to. Name or ID [length(ar.model)]
%
% In this data format, there is one single .tsv-file that contains all data
% points. This data format shall allow easier transitions between modeling
% tools.
% This load-data-function utilizes a different approach than the usual data.def file:
% Firstly, an empty data struct is created which is subsequently appended to the ar struct.
%
% See also arCreateDataStruct arAddDataStruct
%
% References
%   - https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md

function arLoadDataPEtab(datafilename, obsfilename, m)

global ar;

if ~contains(datafilename,'.tsv')
    if ~contains(datafilename,'.')
        datafilename = [datafilename '.tsv'];
    else
        error('this file type is not supported!')
    end
end

if ~contains(obsfilename,'.tsv')
    if ~contains(obsfilename,'.')
        obsfilename = [obsfilename '.tsv'];
    else
        error('this file type is not supported!')
    end
end

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
Tdat = tdfread(datafilename); % all data of the model
fns = fieldnames(Tdat);
for i = 1:length(fns)
    if ischar(Tdat.(fns{i}))
        Tdat.(fns{i}) = regexprep(string(Tdat.(fns{i})),' ','');
    end
end
Tdat = struct2table(Tdat);

Tobs = tdfread(obsfilename); % all data of the model
fns = fieldnames(Tobs);
for i = 1:length(fns)
    if ischar(Tobs.(fns{i}))
        Tobs.(fns{i}) = regexprep(string(Tobs.(fns{i})),' ','');
    end
end
Tobs = struct2table(Tobs);

[uniCond,~,iCCond] = unique(Tdat.simulationConditionId);

%%


%% Use condition specific experiments and distribute over data struct
for iCond = 1:length(uniCond)
    Sd2d = struct();
    args= {};
    % extract important info for data struct from .tsv file
    Tsub = Tdat(iCCond == iCond,:);
    [uniObs,~,iCObs] = unique(cellstr(Tsub.observableId),'stable');
    
    [uniTimes,~,iTExp] = unique(Tsub.time);
    %    [~,ia,ic] = unique([iCobs,iTExp],'rows');
    if length(unique(iTExp(iCObs==1)))<length(iTExp(iCObs==1))
        uniTimes = Tsub.time(1:sum(iCObs==1));
        iTExp = [];
        for i=1:length(uniObs)
            iTExp = [iTExp 1:sum(iCObs==i)];
        end
        iTExp = iTExp';
    end
    uniObs = regexprep(uniObs,' ','');
    Sd2d.name = char(uniCond(iCond));
    Sd2d.tExp = uniTimes;
    Sd2d.tUnits = ar.model.tUnits;
    Sd2d.y = uniObs';
    Sd2d.yNames = uniObs';
    for iObs = 1:length(uniObs)
        idx = strcmp(Tobs.observableId,uniObs{iObs});
        Sd2d.fy{iObs} = char(string(Tobs.observableFormula(idx)));
        Sd2d.fystd{iObs} = char(string(Tobs.noiseFormula(idx)));
        Sd2d.logfitting(iObs) = double(strcmp(Tobs.observableTransformation(idx),'log10'));
    end
    Sd2d.yExpRaw = nan(length(uniTimes),length(uniObs));
    Sd2d.yExpStd = nan(length(uniTimes),length(uniObs));
    Sd2d.yExpStdRaw = nan(length(uniTimes),length(uniObs));
    for it = 1:length(uniTimes)
        for iobs = 1:length(uniObs)
%             [iSim it  iobs]
            if any(it==iTExp & iobs == iCObs)
                Sd2d.yExpRaw(it,iobs) = Tsub.measurement(it == iTExp & iobs == iCObs);
                Sd2d.yExp(it,iobs) = Sd2d.logfitting(iobs) * log10(Tsub.measurement(it == iTExp & iobs == iCObs)) + (1 - Sd2d.logfitting(iobs)) *Tsub.measurement(it == iTExp & iobs == iCObs);
                if isnumeric(Tsub.noiseParameters(it == iTExp & iobs == iCObs))
                    Sd2d.yExpStdRaw(it,iobs) = Tsub.noiseParameters(it == iTExp & iobs == iCObs);
                    Sd2d.yExpStd(it,iobs) =  Sd2d.logfitting(iobs) *log10(Tsub.noiseParameters(it == iTExp & iobs == iCObs)) + (1 - Sd2d.logfitting(iobs))*Tsub.noiseParameters(it == iTExp & iobs == iCObs);
                end
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
    % Deprecated: changes in PEtab
%     for i = 1:length(Sd2d.y)
%         for j = 1:length(ar.model)
%             newpold = ar.model(j).fp(~cellfun('isempty',regexp(ar.model.p,['noiseParameters\d_' Sd2d.y{i} '$'])));
%             newPars = Tsub.noiseParameters(contains(Tsub.observableId,Sd2d.y{i}));
%             if ~(isempty(newpold) || isnumeric(newPars) || isempty(newPars))
%                 newPars = strsplit(str2mat(unique(newPars)),';');
%                 pold = [pold;newpold];fp = [fp;newPars];
%             end
%             newpold = ar.model(j).fp(~cellfun('isempty',regexp(ar.model.p,['observableParameters\d_' Sd2d.y{i} '$'])));
%             newPars = Tsub.observableParameters(contains(Tsub.observableId,Sd2d.y{i}));
%             if ~(isempty(newpold) || isnumeric(newPars) || isempty(newPars))
%                 newPars = strsplit(str2mat(unique(newPars)),';');
%                 pold = [pold;newpold];fp = [fp;newPars];
%             end
%             if length(pold) ~= length(fp)
%                 error('Length of parameter names must be the same!')
%             end
%         end
%     end
    
    D = arCreateDataStruct(m,pold,fp,args{:});
    arAddDataStruct(D,m)
end

end