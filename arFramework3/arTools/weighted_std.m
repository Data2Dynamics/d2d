% weighted standard deviation

function wsx = weighted_std(sx)

wsx = 1 / sum(1./sx.^2);
