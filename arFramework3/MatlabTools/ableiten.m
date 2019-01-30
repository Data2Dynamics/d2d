% out = ableiten(in, dx)
% Berechnet die erste Ableitung

function out = ableiten(in, dx)
out = in*0;
if(exist('dx')~=1) dx = 1; end
n = length(in);

abl = (in(2:n)-in(1:n-1))/dx;

out(1) = abl(1);
out(n) = abl(n-1);

out(2:n-1) = (abl(1:n-2)+abl(2:n-1))/2;
