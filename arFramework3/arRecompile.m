%   arRecompile
%   arRecompile(sameParameterSettings)
%   arRecompile(sameParameterSettings,varargin)
% 
%   This function can be used if the global variable 'ar' is available
%   instead of calling the a setup-routine. The function loads the model(s)
%   and data sets and compiles everything.
% 
%   Model properties cannot be completely taken over!
%   Only model parameters, bounds, qFit, qLog10 etc are taken over by
%   default using arImportPars.m
%   Alternatively set argument sameParameterSettings=false
%   
%   sameParameterSettings   Should parameter Info like (ar.p, ar.lb, ... 
%                           be conserved?
%                           Default: true 
%                           If true, then arImportPars is called to copy
%                           parameter infos.
% 
%   varargin                These arguments are handed over to arCompile
% 
% 
function arRecompile(sameParameterSettings, varargin)
if ~exist('sameParameterSettings','var') || isempty(sameParameterSettings)
    sameParameterSettings = true;
end

global ar

ms = {ar.model.name};
ds = cell(size(ms));
emptyObs = zeros(size(ms));
for m=1:length(ar.model)
    if(isfield(ar.model(m).data,'name'))
        ds{m} = {ar.model(m).data.name};        
        [uni,ia,ib]= unique(regexprep(ds{m},'_nExpID(\d)+',''));
        ds{m} = uni(ib);  % replace zurueck
        ds{m} = ds{m}(sort(ia)); % nur die unique, aber in alter reihenfolge

        for d=1:length(ar.model(m).data)
            emptyObs(m) = emptyObs(m) || sum(sum(~isnan(ar.model(m).data(d).yExp)==0,1))>0;
        end
    else
        ds{m} = [];
    end
end

arIn = arDeepCopy(ar);

try
    
    arInit
    for m=1:length(ms)
        arLoadModel(ms{m});
    end
    
    for m=1:length(ds)
        for d=1:length(ds{m})
            %         arLoadData_withoutNormalization(ds{m}{d}, 1,[],[],[]);
            if emptyObs(m)
                disp('There are empty observations ...')
                arLoadData(ds{m}{d},m,'xls', 0);  % there are no empty observations
            else
                arLoadData(ds{m}{d},m,'xls', 1);  % there are empty observations
            end
        end
    end
    
    arCompileAll(varargin{:})
    
    if sameParameterSettings
        arImportPars(arIn); % use same parameter values, lb, ub, qLog10 etc as in the existing model arIn
    end
    
    % Check, wether data is the same:
    different = false;

    
    for m=1:length(ar.model)
        if(length(ar.model(m).data) ~= length(arIn.model(m).data))
            fprintf('length(ar.model(m).data): %i\nlength(arIn.model(m).data): %i\n', length(ar.model(m).data), length(arIn.model(m).data));

            onlyAr = setdiff({ar.model(m).data.name},{arIn.model(m).data.name});
            onlyArIn = setdiff({arIn.model(m).data.name},{ar.model(m).data.name});
            
            if ~isempty(onlyAr)
                fprintf('data names only in the recompiled ar:\n');
                fprintf('%s, ',onlyAr{:});
                fprintf('\n');
            end
            if ~isempty(onlyArIn)
                fprintf('data names only in the previous variable:\n');
                fprintf('%s, ',onlyArIn{:});
                fprintf('\n');
            end
            warning('Number of data structs is different for ar.model(%i)\nTherefore, data cannot be taken over!\n',m);
            break
        end
        for d=1:length(ar.model(m).data)
            if max(abs(ar.model(m).data(d).yExp-arIn.model(m).data(d).yExp))>1e-10
                different = true;
                break;
            end
            if max(abs(ar.model(m).data(d).yExpStd-arIn.model(m).data(d).yExpStd))>1e-10
                different = true;
                break;
            end
            if max(abs(ar.model(m).data(d).tExp-arIn.model(m).data(d).tExp))>1e-10
                different = true;
                break;
            end
        end
    end
    
    if different
        warining('The data is different! \nThis can occur because of manual normalization.');   
    end
    
    arSave('Recompile');
    
catch ERR
    ar = arIn;
    rethrow(ERR)
end
