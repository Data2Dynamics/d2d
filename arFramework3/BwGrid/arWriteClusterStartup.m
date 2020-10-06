% arWriteClusterStartup(conf)
% arWriteClusterStartup(conf, useSlurm)
% 
% arWriteClusterStartup writes the startup file required on the bwGrid for
% calling moab several times.
% 
%   useSlurm       [true]
%                  In 2020, the bwClusters switched from moab queuing to
%                  slurm. useSlurm=true will execute slurm commands. 
%                  useSlurm=false still writes moab-based scripts.
%
% 
% Example:
% conf = arClusterConfig;
% arWriteClusterStartup(conf)
% 
% See also arFitLhsBwCluster, arClusterConfig

function arWriteClusterStartup(conf,useSlurm)
if ~exist('useSlurm','var') || isempty(useSlurm)
    useSlurm = true;
end

fid = fopen(conf.file_startup,'w');

fprintf(fid,'%s\n','#!/bin/bash');
fprintf(fid,'%s\n',['cd ',conf.pwd]);
fprintf(fid,'%s\n',['for icall in {1..',num2str(conf.n_calls),'}']);
fprintf(fid,'%s\n','do');
if useSlurm
    fprintf(fid,'%s\n',['    ( sbatch ',conf.file_slurm,' $icall ) &']); % calling moab with option icall indicating the call index
else
    fprintf(fid,'%s\n',['    ( msub ',conf.file_moab,' $icall ) &']); % calling moab with option icall indicating the call index
end
fprintf(fid,'%s\n','    sleep 2');
fprintf(fid,'%s\n','done');


fclose(fid);