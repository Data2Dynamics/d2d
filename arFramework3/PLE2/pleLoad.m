% load PLE workspace if existing

function pleLoad

global ar

savepath = [ar.config.savepath '/PLE'];

if(exist(savepath,'dir')==7)
    load([savepath '/results.mat'])
else
    fprintf(1,'\nNo saved workspace for ple found!\n');
end
    

