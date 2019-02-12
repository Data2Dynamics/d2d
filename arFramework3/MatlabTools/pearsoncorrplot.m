% plot correlation of two variables x and y

function [pval, pval2] = pearsoncorrplot(x,y,l,s,~,sy, N)

pval = nan;
pval2 = nan;

if(~exist('s','var') || isempty(s))
    s = 6;
end

if(~exist('N','var') || isempty(N))
    N = 1;
end

% generate test data
if(nargin==0)
    x = 0.5; % correlation
    y = 10; % sample size
    
    figure(1);
    R = mvnrnd([0 0], [1 x; x 1], y);
    x = R(:,1);
    y = R(:,2);
end

% make column vectors
x = x(:);
if(isvector(y))
    y = y(:);
    if(exist('sy','var') && isempty(sy))
        sy = sy(:);
    end
end

% remove nan's
qnonnan = ~isnan(x) & ~isnan(y);
x = x(qnonnan);
y = y(qnonnan);

% plot scatter plots
colors = lines(size(y,2));

for j=1:size(y,2)
    if(~exist('sy','var') || isempty(sy))
        plot(x, y(:,j), '.', 'MarkerSize', s);
        hold on
    else
        qstdzero = sy(:,j)==0;
        errorbar(x(~qstdzero), y(~qstdzero,j), sy(~qstdzero,j), '.', 'MarkerSize', s, 'Color', colors(j,:));
        hold on
        plot(x(qstdzero), y(qstdzero,j), '.', 'MarkerSize', s, 'Color', colors(j,:));
    end
end

% add regression line(s)
if(size(y,1)>2)
    for j=1:size(y,2)
        %     [~,S0,~] = polyfit(x, y(:,j),0);
        [P,S,mu] = polyfit(x, y(:,j),N);
        
        %     pval = 1-chi2cdf(S0.normr - S.normr, 1)
        
        X = linspace(min(x)-(max(x)-min(x))*0.1, max(x)+(max(x)-min(x))*0.1, 50);
        [Y,DELTA] = polyconf(P,X,S, 'alpha', 1-0.68, 'mu', mu);
        
        plot3(X,Y,zeros(size(X))-1,'Color', colors(j,:));
        patch([X, fliplr(X)], [Y+DELTA, fliplr(Y-DELTA)], zeros(size([X, fliplr(X)]))-1, ...
            ones(size([X, fliplr(X)])), 'EdgeColor', 'none', ...
            'FaceColor', colors(j,:), 'FaceAlpha', 0.2);
    end
end

hold off

% add additional labels
if(isvector(y))
    % plot correlation values
    if(length(y)>2)
        [corrval,pval] = corr(x, y);
        [corrval2,pval2] = corr(x, y, 'type', 'Spearman');
        text(0.01,1, sprintf('Pearson %4.2f (p-val %.2g)\nSpearman %4.2f (p-val %.2g)', ...
            corrval, pval, corrval2, pval2), 'Units', 'normalized', ...
            'VerticalAlignment', 'top');
    end
    
    % add labels
    if(nargin>2 && ~isempty(l))
        l = l(qnonnan);
        for j=1:length(l)
            text(x(j),y(j), l{j}, 'FontSize', 8, 'VerticalAlignment', 'top');
        end
    end
end

% adjust axis
arSpacedAxisLimits

xlabel('Prediction');
ylabel('Outcome');

