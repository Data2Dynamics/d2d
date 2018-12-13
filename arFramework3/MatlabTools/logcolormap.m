function Y = logcolormap(X, h)

x = linspace(0, 1, size(X,1))';

if(nargin<2)
    h = 1/2;
else
    h = 1/h;
end

Y = [interp1(x, X(:,1), x.^h) interp1(x, X(:,2), x.^h) interp1(x, X(:,3), x.^h)];

