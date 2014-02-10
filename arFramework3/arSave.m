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
                arSaveParOnly(ar2, ar.config.savepath);
            end
        end
    end
end


if(nargout>0)
    basepath = ar.config.savepath;
end

function arSaveParOnly(ar, savepath) %#ok<INUSL>
save([savepath '/workspace_pars_only.mat'],'ar','-v7.3');
