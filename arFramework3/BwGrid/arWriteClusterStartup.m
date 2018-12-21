% arWriteClusterStartup(conf)
% 
% arWriteClusterStartup writes the startup file required on the bwGrid for
% calling moab several times.
% 
% Example:
% conf = arClusterConfig;
% arWriteClusterStartup(conf)
% 
% See also arFitLhsBwCluster, arClusterConfig

function arWriteClusterStartup(conf)

fid = fopen(conf.file_startup,'w');

fprintf(fid,'%s\n','#!/bin/bash');
fprintf(fid,'%s\n',['cd ',conf.pwd]);
fprintf(fid,'%s\n',['for icall in {1..',num2str(conf.n_calls),'}']);
fprintf(fid,'%s\n','do');
fprintf(fid,'%s\n',['    ( msub ',conf.file_moab,' $icall ) &']); % calling moab with option icall indicating the call index
fprintf(fid,'%s\n','    sleep 2');
fprintf(fid,'%s\n','done');


fclose(fid);