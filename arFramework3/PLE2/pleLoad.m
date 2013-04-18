% load PLE workspace if existing

function pleLoad

global ar

savepath = [ar.config.savepath '/PLE'];

if(exist(savepath,'dir')==7)
    load([savepath '/results.mat'])
    fprintf(1,'PLE workspace loaded from file %s\n', [savepath '/results.mat']);
else
    fprintf(1,'No saved PLE workspace found!\n');
end
    

