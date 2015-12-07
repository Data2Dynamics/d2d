% plot correlation of two variables x and y

function pearsoncorrplot(x,y,l,s)

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
end

% plot scatter plots
for j=1:size(y,2)
    h = scatter(x, y(:,j),'filled');
    h.SizeData = s;
    hold on
end
hold off

% add regression line(s)
lsline

% add additional labels
if(isvector(y))
    % remove nan's
    qnonnan = ~isnan(x) & ~isnan(y);
    x = x(qnonnan);
    y = y(qnonnan);

    % plot correlation values
    [corrval,pval] = corr(x, y);
    text(0.01,1, sprintf('Pearson %4.2f (p-val %4.2f)', corrval,pval), 'Units', 'normalized', ...
        'VerticalAlignment', 'top');
    
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