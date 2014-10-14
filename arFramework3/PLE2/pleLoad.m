% load PLE workspace if existing

function ple = pleLoad(ar)

savepath = [ar.config.savepath '/PLE'];

if(exist(savepath,'dir')==7)
    S = load([savepath '/results.mat']);
    ple = S.pleGlobals;
    fprintf(1,'PLE workspace loaded from file %s\n', [savepath '/results.mat']);
else
    ple = [];
end
    

