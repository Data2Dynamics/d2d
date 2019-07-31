function  arNEBMakeDataFilesCSV(name,steps)
% generates data def files and data tables for NEB method
%
% name      filename of data file template in the Data folder without suffix '_NEB_XXX.csv'
% setps     number of nodes of the band
%
% See also arNEBMakeDataFilesXLSX arNEBMakeDataFilesXLS arNEBMakeModelFiles


path  = [ 'Data/' name ];
path_read  = [ 'Data/' name '_NEB_XXX.def'];    

data_file = [ 'Data/' name '_NEB_XXX.csv'];    

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
    
    copyfile(data_file, [path '_NEB_' nebstring '.csv'])
    
end

% end
Y = strrep(X, '_NEB_XXX', ['_NEB_end' ]) ;
fid2 = fopen([path ['_NEB_end' ] '.def'],'wt') ;

 copyfile(data_file, [path '_NEB_end.csv'])

fwrite(fid2,Y) ;
fclose (fid2) ;


end
