function xdot=xdotTemplate(t,x,uInput,th,f,vAlg)
% This file combines the algebraic relations (e.g. reaction rates) with the
% ODEs that include these relations with v(1), v(2), etc
if isa(uInput,'function_handle')
    u=uInput(t,th);
else
    u=uInput;
end
v=vAlg(x,th,u);
xdot=f(t,x,u,th,v);