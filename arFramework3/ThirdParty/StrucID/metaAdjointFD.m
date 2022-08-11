function dxdt=metaAdjointFD(t,X,uInput,theta,f,vAlg,dim)

nx=dim(1);

% states
x=X(1:nx);

% state dynamics
xdot=xdotTemplate(t,x,uInput,theta,f,vAlg);

% sensitivities
dxdth=reshape(X((nx+1):end),nx,nx);

dfdx = admDiffFD(@(xVar) xdotTemplate(t,xVar,uInput,theta,f,vAlg),1,x);
dxdthdot = -dfdx'*dxdth;

dxdt=[xdot; dxdthdot(:)];