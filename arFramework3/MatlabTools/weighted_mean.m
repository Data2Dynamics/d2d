% weighted mean by their standard deviation

function wmx = weighted_mean(x, sx)

wmx = sum(x./sx.^2) / sum(1./sx.^2);
