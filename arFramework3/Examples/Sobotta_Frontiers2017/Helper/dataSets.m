% Function to generate some variables which help handle the model

% Data groups
coreDatasets = arFindData({'bohl', 'nucSTAT3_Ratio', 'braun_hep', 'braun_app_hep'}, 'names');
app_calibrationdataSets = arFindData('calibration', 'names');
app_validationdataSets = arFindData('validation', 'names');
app_predictiondataSets = arFindData( {'ntSTAT3_triple', 'ntSTAT3_single', 'prediction', 'steady'}, 'names');

% Parameter groups
app = { 'apcs', 'cxcl10', 'fgg', 'hamp', 'hp', 'hpx', 'il33' };
appPars = [];
for a = 1 : length( app )
    names = strcat( {app{a}, app{a}, app{a}, app{a}, app{a}, 'init', 'n'}, {'rna_hill', 'rna_ka', 'rna_pro', 'rna_des', 'rna_delay', app{a}, app{a} } );
    appPars = [appPars arFindPar( ar, names )];
end

calibrationPars = arFindPar( 'calibration');
calibrationPars = union( calibrationPars, appPars );
calibrationPars = union( calibrationPars, arFindPar('knot') );
calibrationPars = union( calibrationPars, arFindPar('qpcr_tc') );
for a = 1 : length( app )
    calibrationPars = union( calibrationPars, arFindPar(sprintf('sd_%s_qpcr', app{a}) ) );
end

validationPars  = arFindPar( 'validation');
relPars         = arFindPar( 'rel' );
predictionPars  = arFindPar( {'ntSTAT3', 'triple', 'predict', 'loc', 'var_il6', 'input_rux', 'steady'} );
