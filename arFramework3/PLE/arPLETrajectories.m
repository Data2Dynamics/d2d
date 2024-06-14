% arPLETrajectories([jks], [n], [show_hit_bound], [saveToFile])
%
% Plot trajectories of ple
%
%   jks               parameter index or vector of par indices [all if not specified]
%                     or alternatively the name of a parameter or a cell of
%                     names can be provided
%   n                 trajectories per parameter (0=all) [10] 
%   show_hit_bound    show hitting boundary of parameters  [false]
%   saveToFile        [false]
%   plotAsRegion      plotting as region [false] 
%                     otherwise (by default) lines are plotted
%
% Plot trajectories for all parameter vectors that were compatible with the
% the confidence intervals determined by ple. pleTrajectories are plotted
% for all variables that were selected for regular plotting, e.g. in
% arPlotter.
%
% See also: ple, arPlotter

function arPLETrajectories(jks, n, show_hit_bound, saveToFile, plotAsRegion)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end
if(~isfield(ar,'ple') || isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 

if(~exist('jks','var') || isempty(jks))
    jks = 1:length(ar.pLabel);
end
if(~exist('n','var') || isempty(n))
    n = 10;
end
if(~exist('show_hit_bound','var') || isempty(show_hit_bound))
    show_hit_bound = false;
end
if(~exist('saveToFile','var') || isempty(saveToFile))
    saveToFile = false;
end
if(~exist('plotAsRegion','var') || isempty(plotAsRegion))
    plotAsRegion = false;
end

ps = [];
chi2s = [];

if ~isfield( ar.ple, 'p_labels' )
    error( 'No ple struct available in ar.ple. Did you perform profile likelihood?' );
end

if(ischar(jks))
    jks = find(ismember(ar.ple.p_labels, jks)); %R2013a compatible
else
    jks = find(ismember(ar.ple.p_labels, ar.pLabel(jks))); %R2013a compatible
end

if ( isempty( jks ) )
    error( 'No profiled parameter selected' );
end

if(length(jks)>2)
    response_str = '_multi';
else
    response_str = '';
end

% collect
for j = jks
    if(j<=length(ar.ple.ps) && ~isempty(ar.ple.ps{j}))
        if(length(jks)<=2)
            if(j==1)
                response_str = ar.ple.p_labels{j};
            else
                response_str = [response_str '_' ar.ple.p_labels{j}];
            end
        end
        
        % all non nans
        chi2stmp = ar.ple.chi2s{j}(~isnan(ar.ple.chi2s{j}));
        pstmp = ar.ple.ps{j}(~isnan(ar.ple.chi2s{j}),:);
        
%         % cut out bounds
%         jk = j;
%         qCloseToUB = pstmp > ones(length(chi2stmp),1) * (ar.ub - ar.ple.dist_thres) & ...
%             ones(length(chi2stmp),1) * ar.qFit==1;
%         qCloseToLB = pstmp < ones(length(chi2stmp),1) * (ar.lb + ar.ple.dist_thres) & ...
%             ones(length(chi2stmp),1) * ar.qFit==1;
%         
%         qhitbound = false(size(pstmp));
%         qhitbound(:,ar.qFit==1) = ar.ple.gradient{jk}(~isnan(ar.ple.chi2s{j}),ar.qFit==1) > ar.ple.grad_thres & qCloseToLB(:,ar.qFit==1) | ...
%             ar.ple.gradient{jk}(~isnan(ar.ple.chi2s{j}),ar.qFit==1) < -ar.ple.grad_thres & qCloseToUB(:,ar.qFit==1);
%         
%         if(~show_hit_bound)
%             chi2stmp = chi2stmp(sum(qhitbound,2)==0);
%             pstmp = pstmp(sum(qhitbound,2)==0,:);
%         end
        
        % cut below threshold
        qbelow = chi2stmp < ar.ple.merit + ar.ple.dchi2_point;
        chi2stmp = chi2stmp(qbelow);
        pstmp = pstmp(qbelow,:);
        
        % make n subset
        if(n>0 && length(chi2stmp) > n)
            iselectn = floor(linspace(1, length(chi2stmp), n));
            chi2stmp = chi2stmp(iselectn);
            pstmp = pstmp(iselectn,:);
        end
        chi2s = [chi2s chi2stmp]; %#ok<*AGROW>
        ps = [ps; pstmp];
    end
end

% % weights
% % minimum chi^2 = 50% alpha, dchi2 = 10%;
% weigths = (chi2s-ar.ple.merit)/ar.ple.dchi2_point;
% weigths = -weigths; 
% weigths = 10.^(weigths);
% weigths = weigths / max(weigths);
% % subplot(2,1,1)
% % plot(chi2s)
% % subplot(2,1,2)
% % plot(weigths)
% % return 

weigths = ones(size(chi2s));

% plot
psnew = ones(size(ps,1),1) * ar.p;
psnew(:,ismember(ar.pLabel, ar.ple.p_labels)) = ps(:,ismember(ar.ple.p_labels, ar.pLabel)); %R2013a compatible


ar.ps_trajectories = psnew;
ar.chi2s_trajectories = chi2s;
arPlotMulti(psnew, weigths, saveToFile, ['_PLE' response_str], plotAsRegion)

