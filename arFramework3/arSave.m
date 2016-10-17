% save model struct and last ple results
% or return base path of last arSave
%
% basepath = arSave(name, withSyms)
% 
% Special case:
% arSave('current')
% saves to the current repository (i.e. current value of ar.config.savepath)


function basepath = arSave(name, withSyms)

global ar pleGlobals

if(~exist('withSyms','var'))
    withSyms = false;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

if(isempty(ar.config.savepath)) % never saved before, ask for name
    if(~exist('name','var'))
        name = input('enter new repository name addition: ', 's');
    end
    if(~isempty(name))
        ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
    else
        ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
    end
    arSaveFull(withSyms);
    arSaveParOnly(ar, ar.config.savepath, pleGlobals);
else
    if(exist('name','var')) % saved before, but new name give
        if(~strcmp(name, 'current'))
            if(~isempty(name))
                ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
            else
                ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
            end
        end
        arSaveFull(withSyms);
        arSaveParOnly(ar, ar.config.savepath, pleGlobals);
        
    else
        if(nargout == 0) % saved before, ask for new name give
            name = input(sprintf('enter new repository name addition [%s]: ', ...
                ar.config.savepath), 's');
            
            if(~isempty(name))
                ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
            end
            
            arSaveFull(withSyms);     
            arSaveParOnly(ar, ar.config.savepath, pleGlobals);
        else
            if(~exist(ar.config.savepath, 'dir')) 
                % saved before, path output requested, 
                % however save path does not exist anymore
                arSaveFull(withSyms);
                arSaveParOnly(ar, ar.config.savepath, pleGlobals);
            end
        end
    end
end

if(nargout>0)
    basepath = ar.config.savepath;
end

% full save
function arSaveFull(withSyms)
global ar
global pleGlobals  %#ok<NUSED>

% remove storage-consuming fields in global ar
% that are non-essential
arCompress;

% check is dir exists
if(~exist(ar.config.savepath, 'dir'))
    mkdir(ar.config.savepath)
end

% remove symy for saving
if(~withSyms)
    for jm = 1:length(ar.model)
        ar.model(jm).sym = [];
        for jc = 1:length(ar.model(jm).condition)
            ar.model(jm).condition(jc).sym = [];
        end
        if(isfield(ar.model(jm), 'data'))
            for jd = 1:length(ar.model(jm).data)
                ar.model(jm).data(jd).sym = [];
            end
        end
    end
end

warning off MATLAB:Figure:FigureSavedToMATFile
save([ar.config.savepath '/workspace.mat'],'ar','pleGlobals','-v7.3');
warning on MATLAB:Figure:FigureSavedToMATFile

fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);



% save only parameters
function arSaveParOnly(ar2, savepath, pleGlobals2)
if nargin>2
    if(isfield(pleGlobals2,'chi2s'))
        pleGlobals.chi2s = pleGlobals2.chi2s;
    else
        pleGlobals = [];
    end
else
    pleGlobals = [];
end
ar = struct([]);
ar(1).pLabel = ar2.pLabel;
ar.p = ar2.p;
ar.qLog10 = ar2.qLog10;
ar.qFit = ar2.qFit;
ar.lb = ar2.lb;
ar.ub = ar2.ub;
ar.type = ar2.type;
ar.mean = ar2.mean;
ar.std = ar2.std;
try %#ok<TRYNC>
    ar.chi2fit = ar2.chi2fit;
    ar.ndata = ar2.ndata;
    ar.nprior = ar2.nprior;
end
try %#ok<TRYNC>
    ar.ps = ar2.ps;
    ar.ps_errors = ar2.ps_errors;
    ar.chi2s = ar2.chi2s;
    ar.chi2sconstr = ar2.chi2sconstr;
    ar.timing = ar2.timing;
    ar.exitflag = ar2.exitflag;
    ar.fun_evals = ar2.fun_evals;
    ar.ps_start = ar2.ps_start;
    ar.chi2s_start = ar2.chi2s_start;
    ar.chi2sconstr_start = ar2.chi2sconstr_start;
    ar.optim_crit = ar2.optim_crit;
    
    ar.chi2s_sorted = ar2.chi2s_sorted;
    ar.chi2sconstr_sorted = ar2.chi2sconstr_sorted;
    ar.ps_sorted = ar2.ps_sorted;
    ar.chi2s_start_sorted = ar2.chi2s_start_sorted;
    ar.chi2sconstr_start_sorted = ar2.chi2sconstr_start_sorted;
    ar.ps_start_sorted = ar2.ps_start_sorted;
end
if(ar2.config.useFitErrorMatrix == 0)
    ar.config.fiterrors = ar2.config.fiterrors; %#ok<STRNU>
else
    ar.config.useFitErrorMatrix = true;
    ar.config.fiterrors_matrix = ar2.config.fiterrors_matrix;
    ar.config.ploterrors_matrix = ar2.config.ploterrors_matrix; 
end
save([savepath '/workspace_pars_only.mat'],'ar','pleGlobals','-v7.3');
