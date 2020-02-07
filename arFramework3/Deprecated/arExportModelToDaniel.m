% arExportModelToDaniel(whichone,m,d,prefix)
function filenames = arExportModelToDaniel(whichone,m,d,prefix)
disp('arExportModelToDaniel is deprecated and will be removed in future releases. Use arExportModelToDmod instead.');

filenames = arExportModelToDmod(whichone,m,d,prefix);
end
