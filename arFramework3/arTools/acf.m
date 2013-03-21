function y = acf(t, data, maxlag)
n = length(data);
maxtau = n/2-1;
if(exist('maxlag')) maxtau = sum(t<maxlag); end

for tau = 0:maxtau
	waitbar(tau/maxtau);
	data1 = data(1:n-tau);
	data2 = data((tau+1):n);
	y(tau+1) = sum(data1.*data2) / n;
end
