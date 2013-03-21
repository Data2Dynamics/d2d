%% Berechnet die zweite Ableitung

function out = ableiten2(in, dx)
xsize = length(in);
out = in*0;
if(exist('dx')~=1) dx = 1; end

for j=2:xsize-1
	out(j) = (in(j-1) - 2*in(j) + in(j+1)) / dx^2;
end
out(1) = out(2);
out(xsize) = out(xsize-1);
