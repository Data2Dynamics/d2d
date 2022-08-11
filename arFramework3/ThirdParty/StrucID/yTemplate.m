function yOut=yTemplate(t,x,uInput,th,h,vAlgebra)
% This file combines the algebraic relations (e.g. reaction rates) with the
% ODEs that include these relations with v(1), v(2), etc
u=uInput(t,th);
v=vAlgebra(x,th,u);
yOut=h(t,x,u,th,v);