% [tmpres, tmpsres] = arRandomEffect(p)
% random effect function, implemention a normal distribution assumption
% with variable mean and std

function [tmpres, tmpsres] = arRandomEffect(p)

tmpres = calclogL(p);

tmpsres = zeros(1,length(p));
dp = 1e-6;
for j=1:length(p)
    pmod = p;
    pmod(j) = pmod(j)+dp;
    
    tmpsres(j) = (calclogL(pmod) - tmpres)/dp;
end

function L = calclogL(p)

pmean = mean(p);
pstd = std(p);

L = sqrt(sum(2*log(pstd) + 50 + ((pmean - p)/pstd).^2));
