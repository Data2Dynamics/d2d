function zdot=zdotReduceTemplate(t,z,uInput,th,f,vAlg,V,k)
% This file combines the algebraic relations (e.g. reaction rates) with the
% ODEs that include these relations with v(1), v(2), etc
% This is the REDUCED model, T is the transformation matrix, k is the order
% of the reduced model (dimension of z)

T=V(:,1:k); % dominant modes of observed dynamics
x=T*z; % back transform from z to x
u=uInput(t,th);
v=vAlg(x,th,u);
zdot=T'*f(t,x,u,th,v);