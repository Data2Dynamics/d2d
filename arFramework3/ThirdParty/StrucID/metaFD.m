function dxdt=metaFD(t,X,uInput,theta,f,vAlg,dim)

nx=dim(1); nth=dim(2);

% states
x=X(1:nx);
% v=vAlg(x,theta,u);

% state dynamics
xdot=xdotTemplate(t,x,uInput,theta,f,vAlg);

% sensitivities
dxdth=reshape(X((nx+1):end),nx,nth);

dfdx = admDiffFD(@(xVar) xdotTemplate(t,xVar,uInput,theta,f,vAlg),1,x);
dfdth = admDiffFD(@(th) xdotTemplate(t,x,uInput,th,f,vAlg),1,theta);

dxdthdot=dfdx*dxdth+dfdth;

dxdt=[xdot; dxdthdot(:)];