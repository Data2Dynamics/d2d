% arLinReg(b, labels, y, x, y_est, covar)
% calculate regression statistics for linear models
%
% after Backhaus, Erichson, Plinke, Weiber "Multivariate Analysemethoden",
% Springer 2006 (11. Auflage)
% chapter 1.2.3, pages: 63-

function arLinReg(b, labels, y, x, y_est, covar)

fprintf('regression statistics for linear models\n\n');

R2 = sum((y_est - mean(y)).^2) / sum((y - mean(y)).^2);
fprintf('coefficient of determination R^2 = %f\n', R2);

K = length(y);
J = length(b);
R2corr = R2 - (J*(1-R2)/(K-J-1));
fprintf('coefficient of determination R^2(corrected) = %f\n\n', R2corr);

F = (R2/J) / ((1-R2) / (K-J-1)); 
F_pval = 1-fcdf(F, J, K-J-1);
fprintf('f-test for non-zero parameters: %f (p-value = %g)\n\n', F, F_pval);

bstd = sqrt(diag(covar))';
t = b./bstd;
t_pval = 1-tcdf(t, K);

fprintf('regression parameters\n');
fprintf('%20s %20s %20s %20s (p-value)\n', 'parameter','value','error','t-test');
for j=1:length(b)
    fprintf('%20s %20g %20g %20g (%f)\n', labels{j}, b(j), bstd(j), t(j), t_pval(j));
end

y_std = std(y);
x_std = std(x);
beta = b .* x_std/y_std;

fprintf('\nregression parameters (standardized)\n');
fprintf('%20s %20s %20s %20s\n', 'parameter','value','value^2', '%');
for j=1:length(b)
    fprintf('%20s %20g %20g %20g\n', labels{j}, beta(j), beta(j)^2, beta(j)^2/sum(beta.^2));
end
fprintf('sum(value^2) = %f\n', sum(beta.^2));







