% revision_path = arGetOldRevision(sha)
%
% Download the D2D revision specified by ar.info.revision
%
%   sha     ID of the different commits/revisions on github
%           (eg '7182eff7eb6f723da7b45da5470a37fd13c9e1de')
%           see commit pages at github
%           the commit sha is stored in ar.info.revision
%           [ar.info.revision]
%
%   revision_path: path where to find OldRevisions
%
% This function downloads an old D2D revision. The target folder is
% .../arFramework3/Advanced/OldRevisions/
%
% For using an older revisions, set your Matlab paths properly, e.g. via
% addpath()

function revision_path = arGetOldRevision(sha)
global ar

if ~exist('sha','var') || isempty(sha)
    if isfield(ar.info,'revision')
        sha = ar.info.revision;
    else
        disp('Revision unknown, may because you are not using git.');
        return;
    end
end

arRemoveOldRevisionPaths

revision_path = [fileparts(which('arInit')),filesep,'Advanced',filesep,'OldRevisons'];
if ~exist(revision_path,'dir')
    mkdir(revision_path);
end


suc3 = false;
suc4 = false;
suc5 = false;
suc6 = false;
suc7 = false;

if ~exist([revision_path,filesep,'d2d-',sha],'dir') % revision not yet available
    fprintf('Downloading %s.zip (may take several minutes)...\n',sha(1:7));
    try
        websave(sprintf('%s.zip',sha(1:7)),sprintf('https://github.com/Data2Dynamics/d2d/archive/%s.zip',sha(1:7)));
        suc1 = true;
    catch
        suc1 = false;
    end
    
    if suc1
        fprintf('unzipping ...\n')
        suc3 = unzip([sha(1:7),'.zip'],revision_path);
        
        if  ~isempty(suc3) % delete large and unnecessary folders
            try
                rmdir([revision_path,filesep,'d2d-',sha,filesep,'arFramework3',filesep,'d2d-presenter'],'s');
            catch
                disp('Could not remove folder d2d-presenter');
            end
            
            try
                rmdir([revision_path,filesep,'d2d-',sha,filesep,'arFramework3',filesep,'Examples'],'s');
            catch
                disp('Could not remove folder Examples');
            end
            
%             try
%                 rmdir([revision_path,filesep,'d2d-',sha,filesep,'arFramework3',filesep,'ThirdParty'],'s');
%             catch
%                 disp('Could not remove folder ThirdParty ');
%             end
%             
%             try
%                 rmdir([revision_path,filesep,'d2d-',sha,filesep,'arFramework3',filesep,'pthreads*'],'s');
%             catch
%                 disp('Could not remove folder pthreads* ');
%             end
%             
%             try
%                 delete([revision_path,filesep,'d2d-',sha,filesep,'arFramework3',filesep,'*.tar']);
%             catch
%                 disp('Could not remove files *.tar');
%             end
            
            try
                delete([sha(1:7),'.zip']);
            catch
                disp('Could not remove sha file.');
            end
        else
            disp('could not unzip the downloaded revision.')
        end
    else
        disp('Could not downloaded the specified old revision. Is wget available on your system?');
    end
    if ~suc1
        fprintf('Could not download revision %s D2D version.\n',sha);
    end
else
    fprintf('Revision %s is already available.\n',sha);
end
