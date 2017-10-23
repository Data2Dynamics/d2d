function testFeature;
    global ar;

    % This file tests whether the description of some relevant fields in
    % arStruct.m is still correct.
    compile = 1;

    if ( compile )
        arInit
        arLoadModel('fieldtest');
        arLoadData('fieldtest');
        arCompileAll();
    end

    initial = 2;
    par = 1;
    gamma = 2.7183;
    beta = 3.1416;

    stateID = find( strcmp( ar.model.x, 'A_state' ) );
    initID = arFindPar( 'init_B_state' );
    parID = arFindPar( 'p1' );
    ar.p( initID ) = initial;
    ar.p( parID ) = par;
   
    ar.qLog10 = ones(size(ar.qLog10));
    ar.model.data.logfitting = 1;
    test( 'Data in log10 space,  parameters in log10 space ', stateID, parID, initID, 10^par * log(10) * gamma, gamma, 1, ( 1 / ( log(10) * 10^initial ) ) );
 
    ar.model.data.logfitting = 0;
    test( 'Data in linear space, parameters in log10 space ', stateID, parID, initID, 10^par * log(10) * gamma, gamma, beta * 10^initial * log(10), beta );

    ar.qLog10 = ones(size(ar.qLog10))*0;
    ar.model.data.logfitting = 1;
    test( 'Data in log10 space,  parameters in linear space', stateID, parID, initID, gamma, gamma, 1 / ( initial * log(10) ), 1 / ( initial * log(10 ) ) );

    ar.model.data.logfitting = 0;
    test( 'Data in linear space, parameters in linear space', stateID, parID, initID, gamma, gamma, beta, beta );
    
end

function result = test( name, stateID, parID, initID, sxfineref, sxexpref, syfineref, syexpref )
    global ar;
    
    arSimu(true,true,true); arCalcMerit(true, ar.p);
    computeRC = @(x, t)(x(4, stateID, parID)-x(1, stateID, parID))/(t(4)-t(1));
    sxfine = computeRC( ar.model.condition.sxFineSimu, ar.model.condition.tFine );
    sxexp  = computeRC( ar.model.condition.sxExpSimu, ar.model.condition.tExp );
    syfine = ar.model.data.syFineSimu(1,:,initID);
    syexp  = ar.model.data.syExpSimu(1,:,initID);
    
    errormsg = @(fname, val, valref) error( 'Failed while testing the case [%s]. %s different than expected (obtained: %g vs expected: %g).\nDid you change D2D in such a way that sxFineSimu is different? Please adjust the test and arStruct.m accordingly before pushing your changes!', name, fname, val, valref );
    
    fprintf( 2, 'Testing case [%s]\t[ ', name );
    
    tol = 1e-8;
    singletest( errormsg, 'sxFineSimu ', sxfine, sxfineref, tol );
    singletest( errormsg, 'syFineSimu ', syfine, syfineref, tol );
    singletest( errormsg, 'sxExpSimu ', sxexp, sxexpref, tol );
    singletest( errormsg, 'syExpSimu ', syexp, syexpref, tol );
    
    fprintf( 2, '] ... [ OK ]\n' );
end

function singletest( errormsg, field, val, ref, tol )
    if ( val - ref ) > tol
        errormsg( field, val, ref );
    else
        fprintf( 2, field );
    end
end