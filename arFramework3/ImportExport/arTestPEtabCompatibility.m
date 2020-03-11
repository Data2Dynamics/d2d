function arTestPEtabCompatibility()
% arTestPEtabCompatibility()
%
% Test PEtab compatibility of current PEtab folder by 
%     1 Importing PEtab model
%     2 Exporting D2D model to PEtab
%     3 Importing Exported model again
%     4 Comparing two d2d models after initial import and export-import cycle
arInit;
arImportPEtab;
arCalcMerit(false)
arSave('FirstImport')
FirstImportStruct = arDeepCopy(ar);
arExportPEtab('d2dExport')
clearvars -except FirstImportStruct;
arInit;
arImportPEtab('d2dExport');
arCalcMerit(false)
arSave('SecondImport');
SecondImportStruct = arDeepCopy(ar);
diary PEtabcheck.txt
same = arCompare(FirstImportStruct,SecondImportStruct,'Main');
if same == 1
   fprintf('The test passed! \nThe two structs after the initial import and the export-import cycle are identical!\n'); 
else
   fprintf('The test failed! \nThe two structs after the initial import and the export-import cycle are NOT identical!\n\n');
end
diary off

end

