function [f g]= rosen10(x)
f=0;
n=10;
for i=1:n-1
    f = f + 100*(x(i)^2 - x(i+1))^2 + (x(i)-1)^2;
end
g=[];
return