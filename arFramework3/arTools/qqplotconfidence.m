function qqplotconfidence(p, runConfInts)

if(~exist('runConfInts','var'))
    runConfInts = false;
end

if(~exist('p','var'))
    p = randn(1,100);
    figure(1);
    runConfInts = true;
end

[p,q] = gety(p);
plot(q, p, 'k.');
hold on
l = max(abs([q(:); p(:)]));
plot([-l l]*1.1, [-l l]*1.1, 'r');

if(runConfInts)
    N = 1000;
    pps = [];
    for j=1:N
        pp = randn(size(p));
        [pp,~] = gety(pp);
        pps(j,:) = pp; %#ok<AGROW>
    end
    ppquants = quantile(pps, [0.025 0.975]);
    patch([q' fliplr(q')], [ppquants(1,:) fliplr(ppquants(2,:))], ...
        -ones(size([q' fliplr(q')])), ...
        'FaceColor','r','FaceAlpha',0.3,'EdgeColor','none');
end

hold off

xlim([-l l]*1.1);
ylim([-l l]*1.1);

xlabel('Standard normal quantiles');
ylabel('Quantiles of parameters');

function [p,q] = gety(p)

p = zscore(p);
p = sort(p);

[n, m] = size(p);
if n == 1
   p = p';
   n = m;
   m = 1;
end

nvec = sum(~isnan(p));
q = repmat((1:n)', 1, m);
q = (q-.5) ./ repmat(nvec, n, 1);
q(isnan(p)) = NaN;

q = norminv(q);

