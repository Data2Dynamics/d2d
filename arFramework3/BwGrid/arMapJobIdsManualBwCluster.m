function arMapJobIdsManualBwCluster(cjId, jobIds, rescolfun)

% arMapJobIdManualBwCluster(cjId, jobIds)
%
% Manually add mapping from workspace ID to cluster job IDs to
% ar.config.cluster. This function is intended to be used when the cluster
% config struct is lost. Please keep in mind that the config is backuped in
% the clConf_[cjId].mat file!
%
%   cjId        computing job ID to assign jobIDs to
%   jobIds      vector of cluster job IDs
%   rescolfun   name of results collection function (including extension)
%
% See also: arHelpBwCluster

global ar

rescolcheck = strsplit(rescolfun, '.');
if rescolcheck{2} ~= 'm'
    error('Please include extension in rescolfun argument')
end

ar.config.cluster.mapping.(cjId).jobIds = jobIds;
ar.config.cluster.mapping.(cjId).rescolscriptfilename = rescolfun;

fprintf('Mapping updated successfully\n')
end