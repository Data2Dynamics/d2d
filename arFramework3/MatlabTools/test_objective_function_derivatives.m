function test_objective_function_derivatives(fkt, pars, dp)

if(~exist('dp','var'))
    dp = 1e-6;
end

pars_reset = pars + 0;

[llh_reset, dllh_returned] = feval(fkt, pars);

dllh_fdiff = nan(size(dllh_returned));

arWaitbar(0);
for j=1:length(pars)
    arWaitbar(j, length(pars));
    pars = pars_reset + 0;
    pars(j) = pars(j) + dp;
    
    llh = feval(fkt, pars) + 0;
    
    dllh_fdiff(:,j) = (llh - llh_reset)/dp;
end
arWaitbar(-1);

figure(1); clf;

subplot(1,2,1)
semilogy(eps + abs(dllh_returned - dllh_fdiff)')
title('Derivatives (absolute differences)')

subplot(1,2,2)
plot(asinh(dllh_returned), '.-')
hold on
plot(asinh(dllh_fdiff), 'o')
hold off
title('Derivatives asinh(values)')
legend('User', 'FD')