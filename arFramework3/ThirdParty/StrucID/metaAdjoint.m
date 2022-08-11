function dxdt=metaAdjoint(t,X,uInput,theta,f,vAlg,dim)

nx=dim(1);
% np=dim(2)-nx;

% states
x=X(1:nx);

% state dynamics
xdot=xdotTemplate(t,x,uInput,theta,f,vAlg);

% sensitivities
dxdth=reshape(X((nx+1):end),nx,nx);

dfdx = admDiffComplex(@(xVar) xdotTemplate(t,xVar,uInput,theta,f,vAlg),1,x);
dxdthdot = -dfdx'*dxdth;

dxdt=[xdot; dxdthdot(:)];