% load PLE workspace if existing

function ple = pleLoad(ar)

savepath = [ar.config.savepath '/PLE'];

if(exist(savepath,'dir')==7 && exist([savepath '/results.mat'],'file')==2)
    S = load([savepath '/results.mat']);
    if isfield(S,'pleGlobals') % old format
        ple = S.pleGlobals; %
        ple.merit = ple.chi2; % pleGlobals.chi2 is now termed ar.ple.merit
    else
        ple = S.ple;
    end
    fprintf(1,'PLE workspace loaded from file %s\n', [savepath '/results.mat']);
else
    ple = [];
end


