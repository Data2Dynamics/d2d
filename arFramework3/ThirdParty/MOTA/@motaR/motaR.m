%
% motaR.m
% des:    create mota object 'motaR' with mota-input 'K' and mota-output 'S', 'rSquared'
%         the object can be investigated with the two functions 'print' and 'plot'
% usage:  m=motaR(K,S,rSquared)
% author: Stefan Hengl 
% year:   2007
%
function  m=motaR(K,S,rSquared)
    
%% checks
if~exist('K','var')||isempty(K)
   error('please provide a (n X p) Matrix as first argument')
end
if~exist('S','var')||isempty(S)
   error('please provide the MOTA-output-matrix S as second argument')
end
if~exist('rSquared','var')||isempty(rSquared)
   error('please provide the MOTA-rSquared-matrix as third argument')
end

% Is S a valid matrix?
for i=1:length(S(:,1))
    for j=1:length(S(:,1))
       if (S(i,j)~=1)&&(S(i,j)~=0) 
          error('no valid MOTA-output-matrix') 
       end
    end    
end
%%

    m.K=K;
    m.S=S;
    m.r2=rSquared;
    
    m=class(m,'motaR');
end