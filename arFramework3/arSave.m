% save model struct and last ple results
% or return base path of last arSave
%
% basepath = arSave(name, withSyms)

function basepath = arSave(name, withSyms)

global ar
global pleGlobals 

if(~exist('withSyms','var'))
    withSyms = false;
end

% remove storage-consuming fields in global ar if ar is huge:
tmp = whos('ar');
if(tmp.bytes > 5e8 || true)
    arUncompressed = arCompress; % large matrices are removed, e.g. sx, su, ...
else 
    arUncompressed = [];
end


if(isempty(ar.config.savepath))
    if(~exist('name','var'))
        name = input('enter new repository name addition: ', 's');
    end
    if(~isempty(name))
        ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
    else
        ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
    end
    
    if(~exist(ar.config.savepath, 'dir'))
        mkdir(ar.config.savepath)
    end
    
    if(~withSyms)
        for jm = 1:length(ar.model)
            ar.model(jm).sym = [];
            if(isfield(ar.model(jm),'condition'))
                for jc = 1:length(ar.model(jm).condition)
                    ar.model(jm).condition(jc).sym = [];
                end
            end
            if(isfield(ar.model(jm), 'data'))
                for jd = 1:length(ar.model(jm).data)
                    ar.model(jm).data(jd).sym = [];
                end
            end
        end
    end
    
    save([ar.config.savepath '/workspace.mat'],'ar','pleGlobals','-v7.3');
    fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
    
    % save only parameters
    ar2 = struct([]);
    ar2(1).pLabel = ar.pLabel;
    ar2.p = ar.p;
    ar2.qLog10 = ar.qLog10;
    ar2.qFit = ar.qFit;
    ar2.lb = ar.lb;
    ar2.ub = ar.ub;
    ar2.type = ar.type;
    ar2.mean = ar.mean;
    ar2.std = ar.std;
    ar2.chi2fit = ar.chi2fit;
    ar2.ndata = ar.ndata;
    ar2.nprior = ar.nprior;
    ar2.config.fiterrors = ar.config.fiterrors;
    
    arSaveParOnly(ar2, ar.config.savepath);
else
    if(exist('name','var'))
        if(~strcmp(name, 'current'))
            if(~isempty(name))
                ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
            else
                ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
            end
        end
        
        if(~exist(ar.config.savepath, 'dir'))
            mkdir(ar.config.savepath)
        end
        
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
        
        save([ar.config.savepath '/workspace.mat'],'ar','pleGlobals','-v7.3');
        fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
        
        % save only parameters
        ar2 = struct([]);
        ar2(1).pLabel = ar.pLabel;
        ar2.p = ar.p;
        ar2.qLog10 = ar.qLog10;
        ar2.qFit = ar.qFit;
        ar2.lb = ar.lb;
        ar2.ub = ar.ub;
        ar2.type = ar.type;
        ar2.mean = ar.mean;
        ar2.std = ar.std;
        ar2.chi2fit = ar.chi2fit;
        ar2.ndata = ar.ndata;
        ar2.nprior = ar.nprior;
        ar2.config.fiterrors = ar.config.fiterrors;
       
        arSaveParOnly(ar2, ar.config.savepath);
        
    else
        if(nargout == 0)
            name = input(sprintf('enter new repository name addition [%s]: ', ...
                ar.config.savepath), 's');
            
            if(~isempty(name))
                ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
            end
            
            if(~exist(ar.config.savepath, 'dir'))
                mkdir(ar.config.savepath)
            end
            
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
            
            save([ar.config.savepath '/workspace.mat'],'ar','pleGlobals','-v7.3');
            fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
            
            % save only parameters
            ar2 = struct([]);
            ar2(1).pLabel = ar.pLabel;
            ar2.p = ar.p;
            ar2.qLog10 = ar.qLog10;
            ar2.qFit = ar.qFit;
            ar2.lb = ar.lb;
            ar2.ub = ar.ub;
            ar2.type = ar.type;
            ar2.mean = ar.mean;
            ar2.std = ar.std;
            ar2.chi2fit = ar.chi2fit;
            ar2.ndata = ar.ndata;
            ar2.nprior = ar.nprior;
            ar2.config.fiterrors = ar.config.fiterrors;

            arSaveParOnly(ar2, ar.config.savepath);
        else
            if(~exist(ar.config.savepath, 'dir'))
                mkdir(ar.config.savepath)
                
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
                
                save([ar.config.savepath '/workspace.mat'],'ar','pleGlobals','-v7.3');
                fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
                
                % save only parameters
                ar2 = struct([]);
                ar2(1).pLabel = ar.pLabel;
                ar2.p = ar.p;
                ar2.qLog10 = ar.qLog10;
                ar2.qFit = ar.qFit;
                ar2.lb = ar.lb;
                ar2.ub = ar.ub;
                ar2.type = ar.type;
                ar2.mean = ar.mean;
                ar2.std = ar.std;
                ar2.chi2fit = ar.chi2fit;
                ar2.ndata = ar.ndata;
                ar2.nprior = ar.nprior;
                ar2.config.fiterrors = ar.config.fiterrors;
                
                arSaveParOnly(ar2, ar.config.savepath);
            end
        end
    end
end

if(~isempty(arUncompressed))
    arUncompressed.config.savepath = ar.config.savepath;
    ar = arUncompressed;
end

if(nargout>0)
    basepath = ar.config.savepath;
end

function arSaveParOnly(ar, savepath) %#ok<INUSL>
save([savepath '/workspace_pars_only.mat'],'ar','-v7.3');



function arIn = arCompress(file)
% Removes some fields from ar to save memory
% 
% if a file is provided, then the compressed struct is only saved 
% i.e. it is overwritten by the uncompressed ar after saving

if(~exist('file','var'))
    file = '';
end

global ar
if(~isempty(file) || nargout>0)
    arIn = ar;
end

if(~isfield(ar.model(1).condition(1),'suFineSimu'))
    arSimu(true,true,true);
end

for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        ar.model(m).condition(c).suFineSimu_dim = size(ar.model(m).condition(c).suFineSimu);
        ar.model(m).condition(c).suFineSimu = [];
        ar.model(m).condition(c).svFineSimu_dim = size(ar.model(m).condition(c).svFineSimu);
        ar.model(m).condition(c).svFineSimu = [];
        ar.model(m).condition(c).sxFineSimu_dim = size(ar.model(m).condition(c).sxFineSimu);
        ar.model(m).condition(c).sxFineSimu = [];
        ar.model(m).condition(c).suExpSimu_dim = size(ar.model(m).condition(c).suExpSimu);
        ar.model(m).condition(c).suExpSimu = [];
        ar.model(m).condition(c).svExpSimu_dim = size(ar.model(m).condition(c).svExpSimu);
        ar.model(m).condition(c).svExpSimu = [];
        ar.model(m).condition(c).sxExpSimu_dim = size(ar.model(m).condition(c).sxExpSimu);
        ar.model(m).condition(c).sxExpSimu = [];
    end
end

ar.isCompressed = 1;

if(~isempty(file))
    if(~isempty(ar.config.savepath))
        pfad = ar.config.savepath;
    else
        pfad = '.';
    end
    [pfad,'/',file]
    save([pfad,'/',file],'ar');
    ar = arIn;
end


