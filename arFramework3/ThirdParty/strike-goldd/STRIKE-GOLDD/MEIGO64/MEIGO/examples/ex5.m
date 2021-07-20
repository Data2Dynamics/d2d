function [J,g,R]=ex5(x,texp,yexp)

[tout,yout] = ode15s(@ex5_dynamics,texp,[100 0 0 0 0],[],x);

R=(yout-yexp);
R=reshape(R,numel(R),1);


J = sum(sum((yout-yexp).^2));
g=0;
return

%***************************************************
%Function of the dynamic system
function dy=ex5_dynamics(t,y,p)

dy=zeros(5,1);  %Initialize the state variables


dy(1)=-(p(1)+p(2))*y(1);
dy(2)=p(1)*y(1);
dy(3)=p(2)*y(1)-(p(3)+p(4))*y(3)+p(5)*y(5);
dy(4)=p(3)*y(3);
dy(5)=p(4)*y(3)-p(5)*y(5);

return
%***************************************************



   
     
   
 
 
 