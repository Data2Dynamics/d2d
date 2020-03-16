function TestFeature()
fprintf( 2,  'INTEGRATION TEST FOR PETAB EXPORT & IMPORT (Boehm 2014)\n' );

fprintf( 'Compiling model...\n' );
Setup_FullModel_Boehm2014
arSimu(true,true,true)
arCalcMerit
ar.chi2
arPlot

d2d_ar = arDeepCopy(ar);

fprintf( '\nExporting to PEtab format...\n' );
arExportPEtab

fprintf( '\nImporting from PEtab format...\n' );
arInit
arImportPEtab('Boehm_JProteomeRes2014')
arSimu(true,true,true)
arCalcMerit
arGetMerit
ar.chi2
arPlot

ar.model.data.y = cellfun(@(x) x(1:end-15), ar.model.data.y, 'UniformOutput', false);
fprintf( '\nComparing ar1, ar2...\n' );
pass = arComparePEtab(ar, d2d_ar, true, true, true, true, false, false, true, false);

if pass
fprintf( '\nPASSED\n' );
else
fprintf( '\nFAILED\n' );
end

system('rm -r PEtab');
end