function [f,g]=ex3(x,k1,k2,k3,k4)

% k1=0.09755988;
% k3=0.0391908;
% k2=0.99*k1;
% k4=0.9*k3;

f=-x(4);

%Equality constraints
g(1)=x(4)-x(3)+x(2)-x(1)+k4*x(4).*x(6);
g(2)=x(1)-1+k1*x(1).*x(5);
g(3)=x(2)-x(1)+k2*x(2).*x(6);
g(4)=x(3)+x(1)-1+k3*x(3).*x(5);

%Inequality constraint
g(5)=x(5).^0.5+x(6).^0.5;

return
