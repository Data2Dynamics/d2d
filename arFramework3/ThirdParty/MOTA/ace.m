% [psi,phi]=ace(x)
% title: Alternating Conditional Expectation Algorithm
% arguments:
%               x : (p,n)-matrix for p parameters and n data points
% values:
%               phi: optimal transformations
% author: Henning Voss (with updates from Stefan Hengl)
function [psi,phi]=ace(x)

    ll=length(x(1,:));
    dim=length(x(:,1))-1;
    wl=3;           % width of smoothing kernel
    oi=100;         % maximum number of outer loop iterations
    ii=10;          % maximum number of inner loop iterations
    ocrit=10*eps;   % numerical zeroes for convergence test
    icrit=1e-4; 
    shol=0;         % 1-> show outer loop convergence, 0-> do not
    shil=0;         % same for inner loop
    
for d=1:dim+1 [xsort(d,:),ind(d,:)]=sort(x(d,:)); end;
for d=1:dim+1 ranks(d,ind(d,:))=[1:ll]; end;
phi=(ranks-(ll-1)/2.)/ sqrt(ll*(ll-1)/12.);
ieps=1.; oeps=1.; oi1=1; ocrit1=1;
while oi1<=oi & ocrit1>ocrit
ii1=1; icrit1=1;
while ii1<=ii & icrit1>icrit
for d=2:dim+1; sum0=0;
  for dd=2:dim+1 if dd ~=d sum0=sum0+phi(dd,:); end; end;
  phi(d,:)=cef(phi(1,:)-sum0,ind(d,:),ranks(d,:),wl,ll);
end;
icrit1=ieps;
if dim==1 sum0=phi(2,:); else sum0=sum(phi(2:dim+1,:)); end;
ieps=sum((sum0-phi(1,:)).^2)/ll;
icrit1=abs(icrit1-ieps);
if shil disp(num2str([ii1 ieps icrit1])); end;
ii1=ii1+1;
end;
phi(1,:)=cef(sum0,ind(1,:),ranks(1,:),wl,ll);
phi(1,:)=(phi(1,:)-mean(phi(1,:)))/std(phi(1,:));
ocrit1=oeps; oeps=sum((sum0-phi(1,:)).^2)/ll; ocrit1=abs(ocrit1-oeps);
if shol disp(num2str([oi1 oeps ocrit1])); end;
oi1=oi1+1;
end;
psi=corrcoef(phi(1,:),sum0); psi=psi(1,2);

function r=cef(y,ind,ranks,wl,ll)
cey=win(y(ind),wl,ll);
r=cey(ranks);

function r=win(y,wl,ll)
% wl=0,1,2...
r=conv(y,ones(2*wl+1,1));
r=r(wl+1:ll+wl)/(2*wl+1);
r(1:wl)=r(wl+1); r(ll-wl+1:ll)=r(ll-wl);
