% Load model
arInit;
arLoadModel('reduction');
arLoadData( 'out', 1, 'csv' );
arLoadData( 'out2', 1, 'csv' );
arCompileAll(true);

