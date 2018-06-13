function Y = logcolormap(X)

x = linspace(0, 1, size(X,1))';

Y = [interp1(x, X(:,1), sqrt(x)) interp1(x, X(:,2), sqrt(x)) interp1(x, X(:,3), sqrt(x))];

