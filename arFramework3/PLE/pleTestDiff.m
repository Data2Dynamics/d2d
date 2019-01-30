% pleTestDiff(jk, n, dp)
% Test differences on parameters if changed by dp
%
% jk    which parameter
% n     how many times to change the parameter jk
% dp    parameter changes
%
% Example:
%    pleTestDiff(3,2,1)
%    plots chi2s and diffmerits of parameter 3, and 2 times adds +1 to
%    value of parameter 3

function pleTestDiff(jk, n, dp)

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end

ps = zeros(1,n);

ps(1) = ar.ple.p(jk);

for j=2:n
	ps(j) = ps(j-1)+dp;
end

chi2s = zeros(1, n);
betas = zeros(n, sum(ar.qFit));
alphas = zeros(n, sum(ar.qFit)^2);

h = waitbar(0, 'Please wait...');
for j=1:n
    h = waitbar(j/n, h);
    
    ptmp = ar.ple.p;
    ptmp(jk) = ps(j);
    
    feval(ar.ple.integrate_fkt, ptmp);
    chi2s(j) = feval(ar.ple.merit_fkt);
    
    [beta, alpha] = feval(ar.ple.diffmerit_fkt);
    betas(j,:) = beta;
    alphas(j,:) = alpha(:);
end
close(h);

%% Plot
figure(1)

subplot(3,1,1)
plot(ps, chi2s, 'x-')
ylabel('chi2s')

subplot(3,1,2)
% plot(ps, betas(:,jk), 'x-')
plot(ps, betas, 'x-')
ylabel('beta')

subplot(3,1,3)
% plot(ps, alphas(:,jk), 'x-')
plot(ps, alphas, 'x-')
ylabel('alpha')
