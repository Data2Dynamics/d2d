Setup

% Fit "normally" (everything below detection limit set to DL/2)
% I don't know why this procedure is ever used, but it often occurs in
% literature.
arFit; arFit; arFit;

states = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I' };
target = [ 5 1 0 10^5 10^1 10^0 5 1 0 ];
for a = 1  : numel( states )
    method1(a) = arGetPars( ['state' states{a}], 0 );
end

% Use the detection limit system (Maximum likelihood approach)
arDetectionLimit(1,1,[1,1,1,1,1,1,1,1,1]);
arFit; arFit; arFit;
for a = 1  : numel( states )
    method2(a) = arGetPars( ['state' states{a}], 0 );
end

% Show the results
disp( '     Abs err (DL/2)     Rel err (DL/2)  Abs err (ML)    Rel err (ML)');  
for a = 1  : numel( states )
    if ( target(a) == 0 )
        fprintf( '%s    %.8g\t\t%.8g\t\t%.8g\t\t%.8g\t\n', states{a}, method1(a)-target(a), (method1(a)-target(a)/target(a)), method2(a)-target(a), (method2(a)-target(a)/target(a)) );
    else
        fprintf( '%s    %.8g\t\t%.8g\t%.8g\t%.8g\t\n', states{a}, method1(a)-target(a), (method1(a)-target(a)/target(a)), method2(a)-target(a), (method2(a)-target(a)/target(a)) );
    end
end

% Plot the sensitivities
checksensis = 0;
if ( checksensis )
    figure;
    arFiniteDifferences;
    plot(ar.sresFD), hold on;
    plot(ar.sres, 'o')
    title( 'Sensitivities' );
    ylabel( 'Sensitivity' );
end