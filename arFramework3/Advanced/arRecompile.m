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

if isfield(ar,'setup')  % only availabel in higher verions
    arIn = arDeepCopy(ar);
    
    try
        for i=1:length(arIn.setup.commands)
            if ~isempty(arIn.setup.arguments{i}) && iscell(arIn.setup.arguments{i})
                feval(arIn.setup.commands{i},arIn.setup.arguments{i}{:}); % function call
            else
                eval(arIn.setup.commands{i})  % no function call
            end
        end
        
        if sameParameterSettings
            arImportPars(arIn); % use same parameter values, lb, ub, qLog10 etc as in the existing model arIn
        end
        
    catch ERR
        ar = arIn;        
        fprintf('\narRecompile fails. global ar is restored as before the function has been called. Error is now rethrown:\n')
        rethrow(ERR)
    end
    
    
else  % setup commands are not stored/not available, e.g. because of older code version
    
    ms = {ar.model.name};
    ds = cell(size(ms));
    emptyObs = zeros(size(ms));
    for m=1:length(ar.model)
        if(isfield(ar.model(m).data,'name'))
            ds{m} = {ar.model(m).data.name};
            try % new code (unfortunately not tested, sorry)
                prands = [ar.model(m).data.prand];
                for ii=1:length(prands)
                    ds{m} = regexprep(ds{m},['_',prands{ii},'(\d)+'],'');
                end
                [uni,ia,ib]= unique(ds{m});
            catch            % old code
                disp('Please check replacement code of random parameters in arRecompile.m')
                [uni,ia,ib]= unique(regexprep(ds{m},'_nExpID(\d)+',''));  % nExpID should in general be replaced by
            end
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
            fprintf('arLoadModel(%s);\n',ms{m});
            arLoadModel(ms{m});
        end
        
        for m=1:length(ds)
            for d=1:length(ds{m})
                %         arLoadData_withoutNormalization(ds{m}{d}, 1,[],[],[]);
                if emptyObs(m)
                    disp('There are empty observations ...')
                    fprintf('arLoadData(%s,%i,''xls'', 0);\n',ds{m}{d},m);
                    arLoadData(ds{m}{d},m,'xls', 0);  % there are no empty observations
                else
                    fprintf('arLoadData(%s,%i,''xls'', 1);\n',ds{m}{d},m);
                    arLoadData(ds{m}{d},m,'xls', 1);  % there are empty observations
                end
            end
        end
        
        fprintf('arCompile(varargin{:});\n')
        arCompileAll(varargin{:});
        
        if sameParameterSettings
            arImportPars(arIn); % use same parameter values, lb, ub, qLog10 etc as in the existing model arIn
        end
        
        
        
    catch ERR
        ar = arIn;
        fprintf('\narRecompile fails. global ar is restored as before the function has been called. Error is now rethrown:\n')
        rethrow(ERR)
    end
    
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
