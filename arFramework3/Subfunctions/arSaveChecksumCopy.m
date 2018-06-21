% This function saves a struct into Folder 'Checksums'. The struct should
% contain all fields which are evaluated for creating the checksum.
% 
%   This function is called by arChecksumData, arChecksumPara and
%   arChecksumFitting.
% 
% 
%   type    'data', 'para', 'fitting'
% 

function arSaveChecksumCopy(arCopy,type,checkstr)

    if ~exist('CheckSums','dir')
        mkdir('CheckSums');
        fid = fopen(['CheckSums',filesep,'Readme.txt'],'w');
        fprintf(fid,'This folder contains the assignments from checksums to ar.\n\n');
        fprintf(fid,'The workspaces contain the fields of global ar which are evaluated for creating the checksum.\n');
        fclose(fid);        
    end
    
    file = ['CheckSums',filesep,type,'_',checkstr,'.mat'];
    
    if ~exist(file,'file')
        ar = arCopy;
        save(file,'ar');
    end
end
