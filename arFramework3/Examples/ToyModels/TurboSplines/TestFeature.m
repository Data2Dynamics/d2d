function TestFeature()

fprintf( 2, 'Testing splines (normal)... ' );
arInit;
arLoadModel('monotone_longer');
arLoadData('test2', 1, 'csv');
ar.model(1).data(1).tExtra = ar.model(1).data(1).tExp;
arCompileAll(true);
save( 'monotone', 'ar' );
fprintf('PASSED\n');
fprintf( 2, 'Checking sensitivities... ' );
arCalcMerit;
arFiniteDifferences(1e-6);
dif = ar.sresFD - ar.sres;
diffmerge = zeros( size( dif ) );
for i=1:size(dif,1), for j=1:size(dif,2), diffmerge(i,j)=min([dif(i,j)/ar.sres(i,j), dif(i,j)]); end, end;
if ( max(max(diffmerge)) < 1e-5 )
    fprintf('PASSED\n');
else
    error( 'FAILED, SENSITIVITY ERROR TOO LARGE\n' );
end
fprintf( 2, 'Fitting... ' );
arFit;
if ( ar.chi2fit < 1e6 )
    fprintf('PASSED\n');
else
    error( 'FAILED, ERROR AFTER FITTING TOO LARGE\n' );
end

fprintf( 2, 'Testing splines (cached splines)... ' );
arInit;
ar.config.turboSplines = 1;
arLoadModel('monotone_turbo');
arLoadData('test2', 1, 'csv');
ar.model(1).data(1).tExtra = ar.model(1).data(1).tExp;
arCompileAll(true);
save( 'monotone_turbo', 'ar' );
fprintf('PASSED\n');
fprintf( 2, 'Checking sensitivities... ' );
arCalcMerit;
arFiniteDifferences(1e-6);
dif = ar.sresFD - ar.sres;
for i=1:size(dif,1), for j=1:size(dif,2), diffmerge(i,j)=min([dif(i,j)/ar.sres(i,j), dif(i,j)]); end, end;
if ( max(max(diffmerge)) < 1e-5 )
    fprintf('PASSED\n');
else
    error( 'FAILED, SENSITIVITY ERROR TOO LARGE\n' );
end
fprintf('PASSED\n');
fprintf( 2, 'Fitting... ' );
arFit;
if ( ar.chi2fit < 1e6 )
    fprintf('PASSED\n');
else
    error( 'FAILED, ERROR AFTER FITTING TOO LARGE\n' );
end
