function yOut=yReduceTemplate(t,z,uInput,th,h,vAlgebra,V,k)
% This file combines the algebraic relations (e.g. reaction rates) with the
% ODEs that include these relations with v(1), v(2), etc
x=V(:,1:k)*z;
u=uInput(t,th);
v=vAlgebra(x,th,u);
yOut=h(t,x,u,th,v);