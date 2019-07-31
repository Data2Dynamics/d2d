function  arNEBMakeModelFiles(name,steps)
% generates model files for NEB method
%
% name      filename of model file template in the Models folder without suffix '_NEB_XXX.def'
% setps     number of nodes of the band
%
% See also arNEBMakeDataFilesCSV arNEBMakeDataFilesXLS arNEBMakeDataFilesXLSX



path  = [ 'Models/' name ];
path_read  = [ 'Models/' name '_NEB_XXX.def'];    


fid = fopen(path_read,'rt');
X = fread(fid) ;
fclose(fid) ;
X = char(X.') ;

%OUT = {};

% replace strings
for i = 1:steps
    nebstring = sprintf('%03d', i);
    Y = strrep(X, '_NEB_XXX', ['_NEB_' nebstring ]) ;
    fid2 = fopen([path '_NEB_' nebstring '.def'],'wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    
   % OUT(end+1) = ['arLoadModel(''' name '_NEB_' nebstring ''');\n'];
    
end

% end
Y = strrep(X, '_NEB_XXX', ['_NEB_end' ]) ;
fid2 = fopen([path ['_NEB_end' ] '.def'],'wt') ;
fwrite(fid2,Y) ;
fclose (fid2) ;

 %   OUT(end+1) = ['arLoadModel(''' name '_NEB_end'');\n'] ;

end
