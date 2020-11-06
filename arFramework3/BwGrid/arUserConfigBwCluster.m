function arUserConfigBwCluster(varargin)

% arUserConfigBwCluster([username, loginNodeUrl, matlabVersion, d2dpath, wd])
%
% Configure cluster access variables:
%   username            SSH username for the login node
%   loginNodeUrl        Login node url
%   matlabVersion       MATLAB version to use on the cluster (e.g. 'R2019b')
%   d2dpath             Absolute path to D2D on the cluster
%   wd                  Working directory on the cluster
%
% Must be called with either no argument (interactive mode) or all of the
% above.
%
% See also: arHelpBwCluster

global ar
if(isempty(ar))
    error('please initialize by arInit')
end

if nargin == 0
    
    fprintf(1,'\n%s','*** D2D cluster configuration ***');
    fprintf(1,'\n(Press return to skip and use value in ar.config.cluster or enter 1 to use default value)\n\n');
    
    if ~isfield(ar.config, 'cluster')
        ar.config.cluster = [];
    end
    
    clusterCreated = isfield(ar.config, 'cluster');
    currentValsStruct = getCurrentVals(ar.config.cluster, {'username', 'loginNodeUrl', 'matlabVersion', 'd2dpath', 'wd'});
    
    
    username = input(sprintf('Please enter your ssh username [return: %s]\n-> ', currentValsStruct.username),'s');
    loginNodeUrl = input(sprintf('Please enter the ssh server url [return: %s, 1: bwforcluster.bwservices.uni-heidelberg.de]\n-> ', currentValsStruct.loginNodeUrl),'s');
    matlabVersion = input(sprintf('Please enter MATLAB version to use on the cluster [return: %s, 1: R2019b]\n-> ', currentValsStruct.matlabVersion),'s');
    d2dpath = input(sprintf('Please enter the absolute path to D2D on the cluster [return: %s, 1: ~/d2d/arFramework3]\n-> ', currentValsStruct.d2dpath),'s');
    wd = input(sprintf('Please enter the working directory on the cluster [return: %s, 1: ~]\n-> ', currentValsStruct.wd),'s');
    
    if exist('username') && ~isempty(username)
        ar.config.cluster.username = username;
    end
    if exist('loginNodeUrl') && ~isempty(loginNodeUrl)
        if loginNodeUrl == '1'
            ar.config.cluster.loginNodeUrl = 'bwforcluster.bwservices.uni-heidelberg.de';
        else
            ar.config.cluster.loginNodeUrl = loginNodeUrl;
        end
    end
    if exist('matlabVersion') && ~isempty(matlabVersion)
        if matlabVersion == '1'
            ar.config.cluster.matlabVersion = 'R2019b';
        else
            if ~isempty(matlabVersion) && isempty(regexp(matlabVersion,'R20[0-9][0-9](a|b)'))
                error('Invalid MATLAB version specified. Use the format R2019b.')
            end
            ar.config.cluster.matlabVersion = matlabVersion;
        end
    end
    if exist('d2dpath') && ~isempty(d2dpath)
        if d2dpath == '1'
            ar.config.cluster.d2dpath = '~/d2d/arFramework3';
        else
            ar.config.cluster.d2dpath = d2dpath;
        end
    end
    if exist('wd') && ~isempty(wd)
        if wd == '1'
            ar.config.cluster.wd = '~';
        else
            ar.config.cluster.wd = wd;
        end
    end
    
elseif nargin == 5
    
    ar.config.cluster.username = varargin{1};
    ar.config.cluster.loginNodeUrl = varargin{2};
    ar.config.cluster.matlabVersion = varargin{3};
    ar.config.cluster.d2dpath = varargin{4};
    ar.config.cluster.wd = varargin{5};
    
else
    error('Either none or all input arguments must be specified.')
end

arCheckFields(ar.config.cluster, {'username', 'loginNodeUrl', 'matlabVersion', 'd2dpath', 'wd'}, 2);
fprintf(1,'\nConfiguration completed\n\n');
end

function out = checkFields(struct, fields)
n = length(fields);
status = zeros(n,1);

if ~isempty(struct)
    for i = 1:n
        if isfield(struct, fields{i}) && ~isempty(struct.(fields{i}))
            status(i) = 1;
        end
    end
end

if sum(status) == n
    out = true;
else
    out = fields(~status);
end
end

function out = getCurrentVals(struct, fields)
    n = length(fields);
if ~isempty(struct)
    for i = 1:n
        if isfield(struct, fields{i}) && ~isempty(struct.(fields{i}))
            out.(fields{i}) = struct.(fields{i});
        else
            out.(fields{i}) = '';
        end
    end
else
    for i = 1:n
        out.(fields{i}) = '';
    end
end
end