% Code based on published code by Marko Laine (2006): MCMC toolbox for Matlab
% http://helios.fmi.fi/~lainema/mcmc/#sec-2
% Plots Geweke's diagnostic for increasing number o
%
% function arGewekePlot(a,b,chain)

function zz=arGewekePlot(varagin)
%GEWEKEPLOT Plot Geweke's diagnostic
% gewekeplot(chain) plots Geweke's diagnostic for increasing number of 
% iterations. See arGeweke.m

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


[nsimu,npar]=size(chain);

n  = 40;
e  = fix(nsimu/2);
l  = fix(e/n);
ii = 1:l:e;

z = zeros(length(ii),npar);
for i=1:length(ii)
  z(i,:) = arGeweke(a,b,chain(ii(i):end,:));
end
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
for i=1:npar
  h = subplot(npar,1,i);
  plot(ii,z(:,i));
  set(h,'XLim',[1 e]);
  if i==1
    title('Geweke diagnostics');
  end
  if i==npar
    xlabel('first iteration used');
  end
end  

if nargout > 0
  zz=z;
end