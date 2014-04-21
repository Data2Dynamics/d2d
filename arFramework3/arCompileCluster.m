function arCompileCluster

global ar

n = parpool('size');

filename = ar.fkt;
fkt_does_exist = nan(1,n);
parfor j=1:n
    fkt_does_exist(j) = exist(filename,'file')>0;
end

if(sum(~fkt_does_exist)>0)
    fprintf('Require compiling on %i workers...\n', sum(~fkt_does_exist));
    fprintf('Transfering files to workers...');
    
    try
        parpool('addAttachedFiles',{'.\Compiled'})
    catch 
        parpool('updateAttachedFiles');
    end
    
    if(~ispc)
        ar_path = strrep(which('arInit.m'),'/arInit.m','');
    else
        ar_path = strrep(which('arInit.m'),'\arInit.m','');
    end
    try
        parpool('addAttachedFiles',{ar_path})
    catch 
        parpool('updateAttachedFiles');
    end
    
    fprintf('done\n');
    
    % TODO pctRunOnAll
    ar1 = ar;
    parfor j=1:n
        if(exist(filename,'file')==0)
            ar2 = ar1;
            arCompile(ar2);
        end
    end
else
    fprintf('%i workers are ready.\n',n);
end

