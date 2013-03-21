% plot trajectories of ple
%
% arPLETrajectories(jks, n, show_hit_bound, saveToFile)
%
% jks               par response                                    [all]
% n                 trajectories per parameter                      [10] (0=all)
% show_hit_bound    show hitting boundary of parameters             [false]
% saveToFile                                                        [false]

function arPLETrajectories(jks, n, show_hit_bound, saveToFile)

global ar
global pleGlobals

if(isempty(ar))
    error('please initialize by arInit')
end
if(isempty(pleGlobals))
    error('PLE ERROR: please initialize')
end 

if(~exist('jks','var') || isempty(jks))
    jks = 1:length(ar.pLabel);
end
if(~exist('n','var'))
    n = 10;
end
if(~exist('show_hit_bound','var'))
    show_hit_bound = false;
end
if(~exist('saveToFile','var'))
    saveToFile = false;
end

ps = [];
chi2s = [];

jks = find(ismember(pleGlobals.p_labels, ar.pLabel(jks)));

if(length(jks)>2)
    response_str = '_multi';
else
    response_str = '';
end

% collect
for j = jks
    if(j<=length(pleGlobals.ps) && ~isempty(pleGlobals.ps{j}))
        if(length(jks)<=2)
            if(j==1)
                response_str = pleGlobals.p_labels{j};
            else
                response_str = [response_str '_' pleGlobals.p_labels{j}];
            end
        end
        
        % all non nans
        chi2stmp = pleGlobals.chi2s{j}(~isnan(pleGlobals.chi2s{j}));
        pstmp = pleGlobals.ps{j}(~isnan(pleGlobals.chi2s{j}),:);
        
%         % cut out bounds
%         jk = j;
%         qCloseToUB = pstmp > ones(length(chi2stmp),1) * (pleGlobals.ub - pleGlobals.dist_thres) & ...
%             ones(length(chi2stmp),1) * pleGlobals.q_fit==1;
%         qCloseToLB = pstmp < ones(length(chi2stmp),1) * (pleGlobals.lb + pleGlobals.dist_thres) & ...
%             ones(length(chi2stmp),1) * pleGlobals.q_fit==1;
%         
%         qhitbound = false(size(pstmp));
%         qhitbound(:,pleGlobals.q_fit==1) = pleGlobals.gradient{jk}(~isnan(pleGlobals.chi2s{j}),pleGlobals.q_fit==1) > pleGlobals.grad_thres & qCloseToLB(:,pleGlobals.q_fit==1) | ...
%             pleGlobals.gradient{jk}(~isnan(pleGlobals.chi2s{j}),pleGlobals.q_fit==1) < -pleGlobals.grad_thres & qCloseToUB(:,pleGlobals.q_fit==1);
%         
%         if(~show_hit_bound)
%             chi2stmp = chi2stmp(sum(qhitbound,2)==0);
%             pstmp = pstmp(sum(qhitbound,2)==0,:);
%         end
        
        % cut below threshold
        qbelow = chi2stmp < pleGlobals.chi2 + pleGlobals.dchi2_point;
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
% weigths = (chi2s-pleGlobals.chi2)/pleGlobals.dchi2_point;
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
psnew(:,ismember(ar.pLabel, pleGlobals.p_labels)) = ps(:,ismember(pleGlobals.p_labels, ar.pLabel));

ar.ps_trajectories = psnew;
ar.chi2s_trajectories = chi2s;
arPlotMulti(psnew, weigths, saveToFile, ['_PLE' response_str])

