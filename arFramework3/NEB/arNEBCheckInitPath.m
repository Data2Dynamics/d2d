function q_pathfound = arNEBCheckInitPath
% evaluates the initial direct path at the nodes and checks if the initial 
% path is already a connecting path

global ar

fprintf(['checking initial NEB path with ' num2str(ar.merger.neb.steps) '\n'] )

%save config
ar.merger.neb.ini.springconst = 0;
    
% set initial path 
ar.p(ar.merger.neb.ps_index) = ar.merger.neb.ps_init;

tic 
% this is where the magic happens
arSimu

fprintf('\n')
    
ar.merger.neb.ini.time = toc;

% collect reults ps
ar.merger.neb.init.ps_result = ar.p(ar.merger.neb.ps_index);

ar.merger.neb.ini.chi2_step = [];

% collect chi2s for each step on the path after fitting the whole path
for im = 1:length(ar.model)
    tmpchi2 = 0;
    for id = 1:length(ar.model(im).data)
        % ToDo: ar.constr
        tmpchi2 = tmpchi2 + sum( ar.model(im).data(id).reserr(:).^2) + sum( ar.model(im).data(id).res(:).^2);
    end
    ar.merger.neb.ini.chi2_step(im) = tmpchi2;
end
    
    % for Checkpath
    ar.merger.neb.chi2_step = ar.merger.neb.ini.chi2_step;
    ar.merger.neb.ps_result = ar.merger.neb.init.ps_result;


%analyze path
q_pathfound = arNEBCheckPath;
ar.merger.neb.ini.q_pathfound = q_pathfound;
    
    if q_pathfound
        ar.merger.neb.spring_out = 0;
    end

%fill up field in struct for saving in ar.merger.lhsrsults(i,j)

ar.merger.neb.springconst = 0;
ar.merger.neb.state = 'off';
ar.merger.neb.res_start = [];
ar.merger.neb.sres_start = [];
ar.merger.neb.res_spring = [];
ar.merger.neb.sres_spring = [];
ar.merger.neb.res_end = [];
ar.merger.neb.sres_end = [];


end
