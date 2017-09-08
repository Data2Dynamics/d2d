if ( 0 )
    arInit;
    arLoadModel('reduction');
    arLoadData( 'out', 1, 'csv' );
    arCompileAll(true);
    arSimu(false, false, true);
    x  = ar.model.condition.xExpSimu;
    sx = ar.model.condition.sxExpSimu;
end

arInit;
arLoadModel('reduction');
arReduce;
arLoadData( 'out', 1, 'csv' );
arCompileAll(true);
arSimu(false, true, true);
x_r  = ar.model.condition.xExpSimu;
sx_r = ar.model.condition.sxExpSimu;
