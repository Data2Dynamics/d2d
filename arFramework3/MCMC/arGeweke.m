% [z,p] = arGeweke([a],[b],[chain]) 
%
% Performs Geweke test on chain and returns test values
%
%   a       Bound for lower fraction of chain       [0.1]
%   b       Bound for upper fraction of chain       [0.5]
%   chain   sample of MCMC chain analysed           [ar.ps]
%
%   z       Z-values of Geweke test
%   p       p-Values of Geweke test
%
% See also arGewekePlot
%
% Code based on published code by Marko Laine (2006): MCMC toolbox for Matlab
% http://helios.fmi.fi/~lainema/mcmc/#sec-2
%
% See:
% Stephen P. Brooks and Gareth O. Roberts.
% Assessing convergence of Markov chain Monte Carlo algorithms.
% Statistics and Computing, 8:319--335, 1998.



function [z,p]=arGeweke(varargin)
%GEWEKE Geweke's MCMC convergence diagnostic
% [z,p] = arGeweke(a, b, chain)
% Test for equality of the means of the first a% (default 10%) and
% last b% (50%) of a Markov chain.
% See:
% Stephen P. Brooks and Gareth O. Roberts.
% Assessing convergence of Markov chain Monte Carlo algorithms.
% Statistics and Computing, 8:319--335, 1998.

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

global ar

a = 0.1;
b = 0.5;

% Read out given variables
if nargin > 2
	chain = varargin{3};
else
    chain = ar.ps;
end
if nargin > 1
	b = varargin{2};
end
if nargin > 0
    a = varargin{1};    
end


[nsimu,~]=size(chain);

na = floor(a*nsimu);
nb = nsimu-floor(b*nsimu)+1;

if (na+nb-1)/nsimu > 1
  error('Error with na and nb');
end

m1 = mean(chain(1:na,:));
m2 = mean(chain(nb:end,:));

%%% Spectral estimates for variance
sa = spectrum0(chain(1:na,:));
sb = spectrum0(chain(nb:end,:));

z = (m1-m2)./(sqrt(sa/na+sb/(nsimu-nb+1)));
p = 2*(1-nordf(abs(z)));




function s=spectrum0(x)
%SPECTRUM0 Spectral density at frequency zero
% spectrum0(x) spectral density at zero for columns of x

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

[m,n]= size(x);
s = zeros(1,n);
for i=1:n
  spec = spectrum(x(:,i),m);
  s(i) = spec(1);
end

function [y,f]=spectrum(x,nfft,nw)
%SPECTRUM Power spectral density using Hanning window
%  [y,f]=spectrum(x,nfft,nw) 

% See also: psd.m in Signal Processing Toolbox 

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:40 $

if nargin < 2 || isempty(nfft)
  nfft = min(length(x),256);
end
if nargin < 3 || isempty(nw)
  nw = fix(nfft/4);
end
noverlap = fix(nw/2);

% Hanning window
w = .5*(1 - cos(2*pi*(1:nw)'/(nw+1)));
% Daniel
%w = [0.5;ones(nw-2,1);0.5];
n = length(x);
if n < nw
    x(nw)=0;  n=nw;
end
x = x(:);

k = fix((n-noverlap)/(nw-noverlap)); % no of windows
index = 1:nw;
kmu = k*norm(w)^2; % Normalizing scale factor
y = zeros(nfft,1);
for i=1:k
% xw = w.*detrend(x(index),'linear');
  xw = w.*x(index);
  index = index + (nw - noverlap);
  Xx = abs(fft(xw,nfft)).^2;
  y = y + Xx;
end

y = y*(1/kmu); % normalize

n2 = floor(nfft/2);
y  = y(1:n2);
f  = 1./n*(0:(n2-1));

function y=nordf(x,mu,sigma2)
% NORDF the standard normal (Gaussian) cumulative distribution.
% NORPF(x,mu,sigma2) x quantile, mu mean, sigma2 variance

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.6 $  $Date: 2012/09/27 11:47:38 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end

%y = 0.5*erf(-Inf,sqrt(2)*0.5*x);
y = 0.5+0.5*erf((x-mu)/sqrt(sigma2)/sqrt(2));

