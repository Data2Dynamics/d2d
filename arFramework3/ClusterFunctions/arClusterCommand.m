% varargout = arClusterCommand(cluster, command, [folder])
% Prints out information on your cluster calculation (job, id, parallel,..)
%
% cluster:      MATLAB cluster object       (see help parcluster)

function varargout = arClusterCommand(cluster, command, folder)

global ar_command_cluster
global ar_command_cluster_command

if(isempty(ar_command_cluster)) % new job
    if(nargin==0)
        error('specify cluster!');
    end
    if(nargout>0)
        error('no job available');
    end
    
    if(~exist('folder', 'var'))
        folder = '.';
    end

    fprintf('arClusterCommand sending job...');
    ar_command_cluster = batch(cluster, @arClusterCommandFun, 2, ...
        {command}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', folder);
    ar_command_cluster_command = [folder '> ' command];
    
    fprintf('done\n');
    
elseif(isa(ar_command_cluster,'parallel.job.MJSIndependentJob') || ...
        isa(ar_command_cluster,'parallel.job.MJSCommunicatingJob') || ...
        isa(ar_command_cluster,'parallel.job.CJSIndependentJob') || ...
        isa(ar_command_cluster,'parallel.job.CJSCommunicatingJob')) % old job
    if(nargout>0)
        varargout{1} = ar_command_cluster;
        return
    end
    try
        fprintf('arClusterCommand (ID %i) ', ar_command_cluster.ID);
    catch
        fprintf('arClusterCommand invalid job ID, deleting...\n');
        clear global ar_command_cluster
        return
    end
    
    if(~strcmp(ar_command_cluster.State, 'finished')) % still running
        fprintf('%s, status %s...\n', ar_command_cluster_command, ar_command_cluster.State);
    else % finished
        try
            fprintf('retrieving results...\n\n');
            S = fetchOutputs(ar_command_cluster);
            fprintf('%s\n',ar_command_cluster_command);
            disp(S{2})
            fprintf('(exit code %i)\n',S{1});
            delete(ar_command_cluster);
            clear global ar_command_cluster
        catch err_id
            delete(ar_command_cluster);
            clear global ar_command_cluster
            rethrow(err_id);
        end
    end
else
    error('arClusterCommand global variable ar_command_cluster is invalid!\n');
end


function [exitcode, out] = arClusterCommandFun(command)
[exitcode, out] = system(command);
