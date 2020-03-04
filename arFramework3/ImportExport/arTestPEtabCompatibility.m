function arTestPEtabCompatibility()
% arTestPEtabCompatibility()
%
% Test PEtab compatibility of current PEtab folder by 
%     1 Importing PEtab model
%     2 Exporting D2D model to PEtab
%     3 Importing Exported model again
%     4 Comparing two d2d models after initial import and export-import cycle
arInit
arImportPEtab
arSimu(0,1,0)
arSave('FirstImport')
FirstImportStruct = arDeepCopy(ar);
arExportPEtab('d2dExport')
clearvars -except FirstImportStruct
arInit
arImportPEtab('d2dExport')
arSimu(0,1,0)
arSave('SecondImport')
SecondImportStruct = arDeepCopy(ar);
arCompare(FirstImportStruct,SecondImportStruct,'main')

end

