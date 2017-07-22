% Small benchmark to demonstrate speed difference when caching spline
% coefficients. To enable this mode, set ar.config.turboSplines to 1. Note
% that this enforces recompilation of conditions every time, since this is
% used to determine number of unique splines appearing in the sensitivity
% equations.

load('monotone');
tic;
for a = 1 : 5
    ar.p=-1*ones(size(ar.p)) + a;
    ar.config.maxstepsize=1e-1; % Set a lower stepsize to pretend it's like a model that requires many steps
    arCalcMerit;
end
c = toc;
fprintf( 'Time taken without caching %d\n', c );

load( 'monotone_turbo' );
tic;
for a = 1 : 5
    ar.p=-1*ones(size(ar.p)) + a;
    ar.config.maxstepsize=1e-1;
    arCalcMerit;
end
c = toc;
fprintf( 'Time taken with caching %d\n', c );