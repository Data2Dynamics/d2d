function arNEBFit(springconst)
% NEB fitting: Add residuals for spring forces and relax (fit) band from
% initial path
% projection of gradients is done in arCollectRes

global ar

%save config
ar.merger.neb.springconst = springconst^2 * 2; % corrected spring constant (because of quadratic terms in S, etc.)

% residuals for springs
arAddCustomResidual( 'NEB_res_node_start1', @user_res_start, 1);
arAddCustomResidual( 'NEB_res_node_start2', @user_res_start, 1); % double spring strength for fixed end note to first free node
arAddCustomResidual( 'NEB_res_node_start3', @user_res_start, 1); % double spring strength for fixed end note to first free node

 for inode = 1:ar.merger.neb.steps    
     arAddCustomResidual( ['NEB_res1_node_' sprintf('%03d', inode)], @()user_res1(inode) );
     arAddCustomResidual( ['NEB_res2_node_' sprintf('%03d', inode)], @()user_res2(inode) );
 end
 
arAddCustomResidual( 'NEB_res_node_end1', @user_res_end, 1);
arAddCustomResidual( 'NEB_res_node_end2', @user_res_end, 1); % double spring strength for fixed end note to first free node
arAddCustomResidual( 'NEB_res_node_end3', @user_res_end, 1); % double spring strength for fixed end note to first free node

try 
    
    ar.merger.neb.state = 'on';
    ar.merger.neb.chi2_step = 0;
    
    % set initial path 
    ar.p(ar.merger.neb.ps_index) = ar.merger.neb.ps_init;
    
    % this is where the magic happens
    arFit(1) % silent
    
    % collect results ps
    ar.merger.neb.ps_result = ar.p(ar.merger.neb.ps_index);
    
    ar.merger.neb.chi2_step = [];
    % collect chi2s for each step on the path after fitting the whole path
    for im = 1:length(ar.model)
        tmpchi2 = 0;
        for id = 1:length(ar.model(im).data)
            % ToDo: ar.constr
            tmpchi2 = tmpchi2 + sum( ar.model(im).data(id).reserr(:).^2) + sum( ar.model(im).data(id).res(:).^2);
        end
        ar.merger.neb.chi2_step(im) = tmpchi2;
    end
    
catch ERR
    resetResidualsAndState;
    rethrow(ERR)
end
resetResidualsAndState;
   
end



function resetResidualsAndState

global ar

    ar.merger.neb.state = 'off';

    arRemoveCustomResidual( 'NEB_res_node_start1');
    arRemoveCustomResidual( 'NEB_res_node_start2');
    arRemoveCustomResidual( 'NEB_res_node_start3');

    for inode = 1:ar.merger.neb.steps    
        arRemoveCustomResidual( ['NEB_res1_node_' sprintf('%03d', inode)]);
        arRemoveCustomResidual( ['NEB_res2_node_' sprintf('%03d', inode)]);
    end

    arRemoveCustomResidual( 'NEB_res_node_end1');
    arRemoveCustomResidual( 'NEB_res_node_end2');
    arRemoveCustomResidual( 'NEB_res_node_end3');

end


function [res_user, res_type, sres_user] = user_res_start
global ar

sres_user(1,1:length(ar.p)) = 0;
if strcmp(ar.merger.neb.state,'on')

    res_type = 1; % like a data residual
    
    % start p_0 (or p)
    p_start = ar.p(ar.merger.neb.ps_index(1,:));
    % next p_1
    p_poststart = ar.p(ar.merger.neb.ps_index(2,:));

    % res 
    tmpres = sqrt(sum((p_poststart - p_start).^2));
    res_user =  ar.merger.neb.springconst * tmpres;
    
    %all others 0
    sres_user(1,ar.qFit==1) = 0;
    
    % start
    sres_user(1, ar.merger.neb.ps_index(1,:)) = ...
        -2*(p_poststart-p_start) ./ (2*sqrt(sum((p_poststart-p_start).^2)));
    
    % next
    sres_user(1, ar.merger.neb.ps_index(2,:)) = ...
        2*(p_poststart-p_start) ./ (2*sqrt(sum((p_poststart-p_start).^2)));

    sres_user(isnan(sres_user)) = 0;
    
else % Number of data points should be the same
    res_user = 0;
    sres_user = zeros(1,length(ar.p));
    res_type = 1;
end

ar.merger.neb.res_start = res_user;
ar.merger.neb.sres_start = sres_user;
end


function [res_user, res_type, sres_user] = user_res1(inode)
global ar

% indexing in ar.merger.neb.ps_index is different from node indexing
inode = inode+1;

sres_user(1,1:length(ar.p)) = 0;
if strcmp(ar.merger.neb.state,'on')

    res_type = 1; % like a data residual
    
    % current
    p_c = ar.p(ar.merger.neb.ps_index(inode,:));
    % next
    p_n = ar.p(ar.merger.neb.ps_index(inode+1,:));
    % previous
    p_p = ar.p(ar.merger.neb.ps_index(inode-1,:));
    
    %naive tangente
    tmpres = sum((p_c-p_p).*(p_n-p_p)) ./ sqrt( sum( (p_n-p_p).^2 ) );
    res_user =  ar.merger.neb.springconst * tmpres;
    
    %all others 0
    sres_user(1,ar.qFit==1) = 0;
    
    % current
    sres_user(1, ar.merger.neb.ps_index(inode,:)) = (p_n - p_p)./ sqrt( sum( (p_n-p_p).^2 )) ;
    
    % next
    sres_user(1, ar.merger.neb.ps_index(inode+1,:)) = ...
        (p_c-p_p) ./ sqrt(sum((p_n-p_p).^2 )) ... 
        - (2.*(p_n-p_p).*sum((p_c-p_p).*(p_n-p_p)) ...
        / (2*sqrt(sum((p_n-p_p).^2 ))^3)) ;
    
    % previous
    sres_user(1, ar.merger.neb.ps_index(inode-1,:)) = ...
        2.*(p_n-p_p).* sum((p_c-p_p).*(p_n-p_p))  ...
        / (2 * sqrt(sum((p_n-p_p).^2 ))^3) ...
        - (p_c + p_n - 2.*p_p) ./ sqrt(sum((p_n-p_p).^2 )) ;
    
    sres_user(isnan(sres_user)) = 0;
    
else % Number of data points should be the same
    res_user = 0;
    sres_user = zeros(1,length(ar.p));
    res_type = 1;
end
 ar.merger.neb.res_spring(inode*2-1) = res_user;
 ar.merger.neb.sres_spring(inode*2-1,:) = sres_user;
end




function [res_user, res_type, sres_user] = user_res2(inode)
global ar

% indexing in ar.merger.neb.ps is different from node indexing
inode = inode+1;

sres_user(1,1:length(ar.p)) = 0;
if strcmp(ar.merger.neb.state,'on')

    res_type = 1; % like a data residual
    
    % current
    p_c = ar.p(ar.merger.neb.ps_index(inode,:));
    % next
    p_n = ar.p(ar.merger.neb.ps_index(inode+1,:));
    % previous
    p_p = ar.p(ar.merger.neb.ps_index(inode-1,:));
       
    %naive tangente
    tmpres = sum((p_c-p_n).*(p_n-p_p)) ./ sqrt( sum( (p_n-p_p).^2 ) );
    res_user =  ar.merger.neb.springconst * tmpres;
    
    %all others 0
    sres_user(1,ar.qFit==1) = 0;
    
    % current
    sres_user(1, ar.merger.neb.ps_index(inode,:)) = (p_n - p_p)./ sqrt(sum((p_n-p_p).^2));
    
    % next
    sres_user(1, ar.merger.neb.ps_index(inode+1,:)) = ...
        (p_c - 2.*p_n + p_p) ./ sqrt(sum((p_n-p_p).^2 )) ... 
        - 2.*(p_n-p_p).*sum((p_c-p_n).*(p_n-p_p)) ...
        / (2*sqrt(sum((p_n-p_p).^2 ))^3);
    
    % previous
    sres_user(1, ar.merger.neb.ps_index(inode-1,:)) = ...
        2.*(p_n-p_p).* sum((p_c-p_n).*(p_n-p_p)) ...
        / (2 * sqrt(sum((p_n-p_p).^2 ))^3) ...
        - (p_c - p_n) ./ sqrt(sum((p_n-p_p).^2 )) ;
    
    sres_user(isnan(sres_user)) = 0;
    
else % Number of data points should be the same
    res_user = 0;
    sres_user = zeros(1,length(ar.p));
    res_type = 1;
end
 ar.merger.neb.res_spring(inode*2) = res_user;
 ar.merger.neb.sres_spring(inode*2,:) = sres_user;
end


function [res_user, res_type, sres_user] = user_res_end
global ar

sres_user(1,1:length(ar.p)) = 0;
if strcmp(ar.merger.neb.state,'on')

    res_type = 1; % like a data residual
    
    % end p_end
    p_end = ar.p(ar.merger.neb.ps_index(end,:));
    % last before end: p_end-1
    p_preend = ar.p(ar.merger.neb.ps_index(end-1,:));

    % res end
    tmpres = sqrt(sum((p_preend - p_end).^2));
    res_user =  ar.merger.neb.springconst * tmpres;
    
    %all others 0
    sres_user(1,ar.qFit==1) = 0;

    
    % end
    sres_user(1, ar.merger.neb.ps_index(end,:)) = ...
        2*(p_end - p_preend) ./ (2*sqrt(sum((p_end-p_preend).^2)));
    
    % preend
    sres_user(1, ar.merger.neb.ps_index(end-1,:)) = ...
        -2*(p_end - p_preend) ./ (2*sqrt(sum((p_end-p_preend).^2)));

    sres_user(isnan(sres_user)) = 0;
    
else % Number of data points should be the same
    res_user = 0;
    sres_user = zeros(1,length(ar.p));
    res_type = 1;
end
ar.merger.neb.res_end = res_user;
ar.merger.neb.sres_end = sres_user;
end

