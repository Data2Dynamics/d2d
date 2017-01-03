function inv = chi2inv (x, n)
if (nargin ~= 2)
    error ('chi2inv: you must give two arguments');
end

if (~isscalar (n))
    [retval, x, n] = common_size(x, n);
    if (retval > 0)
        error ('chi2inv: x and n must be of common size or scalar');
    end
end

inv = gaminv(x, n / 2, 2);

function inv = gaminv (x, a, b)
if (nargin ~= 3)
    error ('gaminv: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gaminv: x, a and b must be of common size or scalars');
    end
end

sz = size (x);
inv = zeros (sz);

k = find ((x < 0) | (x > 1) | isnan (x) | ~(a > 0) | ~(b > 0));
if (any (k))
    inv (k) = NaN;
end

k = find ((x == 1) & (a > 0) & (b > 0));
if (any (k))
    inv (k) = Inf;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
    if (~isscalar(a) || ~isscalar(b))
        a = a (k);
        b = b (k);
        y = a .* b;
    else
        y = a * b * ones (size (k));
    end
    x = x (k);
    l = find (x < eps);
    if (any (l))
        y(l) = sqrt (eps) * ones (length (l), 1);
    end
    
    y_old = y;
    for i = 1 : 100

        h     = (gamcdf (y_old, a, b) - x) ./ gampdf (y_old, a, b);
        y_new = y_old - h;
        ind   = find (y_new <= eps);
        if (any (ind))
            y_new (ind) = y_old (ind) / 10;
            h = y_old - y_new;
        end
        if (max (abs (h)) < sqrt (eps))
            break;
        end
        y_old = y_new;
    end
    
    inv (k) = y_new;
end

function pdf = gampdf (x, a, b)
if (nargin ~= 3)
    error ('gampdf: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gampdf: x, a and b must be of common size or scalars');
    end
end

sz = size(x);
pdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    pdf (k) = NaN;
end

k = find ((x > 0) & (a > 0) & (a <= 1) & (b > 0));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = (x(k) .^ (a - 1)) ...
            .* exp(- x(k) ./ b) ./ gamma (a) ./ (b .^ a);
    else
        pdf(k) = (x(k) .^ (a(k) - 1)) ...
            .* exp(- x(k) ./ b(k)) ./ gamma (a(k)) ./ (b(k) .^ a(k));
    end
end

k = find ((x > 0) & (a > 1) & (b > 0));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = exp (- a .* log (b) + (a-1) .* log (x(k)) ...
            - x(k) ./ b - gammaln (a));
    else
        pdf(k) = exp (- a(k) .* log (b(k)) + (a(k)-1) .* log (x(k)) ...
            - x(k) ./ b(k) - gammaln (a(k)));
    end
end




function cdf = gamcdf (x, a, b)
if (nargin ~= 3)
    error ('gamcdf: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gamcdf: x, a and b must be of common size or scalars');
    end
end

sz = size (x);
cdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    cdf (k) = NaN;
end

k = find ((x > 0) & (a > 0) & (b > 0));
if (any (k))
    if (isscalar (a) && isscalar(b))
        cdf (k) = gammainc (x(k) ./ b, a);
    else
        cdf (k) = gammainc (x(k) ./ b(k), a(k));
    end
end
