function TestFeature()

fprintf( 'INTEGRATION TEST FOR SBML IMPORT (Dalla Man 2007)\n' );

fprintf( 2, 'Convert model from SBML... ' );
arInit;
arImportSBML('BIOMD0000000379', 'tend', 100);
fprintf( 'PASSED\n' );

fprintf( 2, 'Compiling model from SBML... ' );
ar.config.checkForNegFluxes = false;

arLoadModel('BIOMD0000000379');
arLoadData('BIOMD0000000379_data');
arCompileAll(true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Testing concordance with biomodels simulation result... ' );

arSimu(true,true,true); arChi2;

if ( arCompareWithBiobaseSimulation('SIMU1448030646732.dat', true ) )
    fprintf('PASSED\n');
else
    error( 'FINAL ERROR TOO LARGE' );
end

