% Instructions for using the BwCluster functions located in BwGrid/
%
% This folder contains a collection of functions written for parallelizing
% multistart optimization (arFitLHS) over the individual multistart runs,
% as well as profile likelihood calculation (ple) over the two branches of
% all individual profiles on bwForCluster MLS&WISO, the high performance
% computing system for molecular life sciences as wenn as economic and
% social studies of the state of Baden-Württemberg [1], which uses the
% slurm queueing system.
%
% The functions in BwGrid/ can be divided into two categories:
%
%   (1) Functions that can be called from a login node of the cluster
%       - arFitLhsBwCluster
%       - pleBwCluster
%
%   (2) Functions that can called from a local machine with access to the
%   cluster via ssh
%       - arUploadBwCluster
%       - arSendJobBwCluster
%       - arJobStatusBwCluster
%       - arKillJobBwCluster
%       - arDownloadResultsBwCluster
%
% References
%[1] https://wiki.bwhpc.de/e/Category:BwForCluster_MLS%26WISO_Production
help arHelpBwCluster
 
