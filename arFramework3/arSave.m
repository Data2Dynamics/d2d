% save model struct and last ple results
% or return base path of last arSave
%
% basepath = arSave(name, withSyms)

function basepath = arSave(name, withSyms)

global ar
global pleGlobals 
global mledist

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
    
    save([ar.config.savepath '/workspace.mat'],'-v7.3');
    fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
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
        
        save([ar.config.savepath '/workspace.mat']);
        fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
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
            
            save([ar.config.savepath '/workspace.mat']);
            fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
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
                
                save([ar.config.savepath '/workspace.mat']);
                fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);
            end
        end
    end
end


if(nargout>0)
    basepath = ar.config.savepath;
end

