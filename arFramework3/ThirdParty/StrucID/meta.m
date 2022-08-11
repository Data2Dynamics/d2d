function dxdt=meta(t,X,uInput,theta,f,vAlg,dim)

nx=dim(1); nth=dim(2);

% states
x=X(1:nx);

% state dynamics
xdot=xdotTemplate(t,x,uInput,theta,f,vAlg);

% sensitivities
dxdth=reshape(X((nx+1):end),nx,nth);

dfdx = admDiffComplex(@(xVar) xdotTemplate(t,xVar,uInput,theta,f,vAlg),1,x);
% dfdx = AutoDiffJacobian(@(xVar) xdotTemplate(t,xVar,uInput,theta,f,vAlg),x,'AutoDiff');
dfdth = admDiffComplex(@(th) xdotTemplate(t,x,uInput,th,f,vAlg),1,theta);
% dfdth = AutoDiffJacobian(@(th) xdotTemplate(t,x,uInput,th,f,vAlg),theta,'AutoDiff');

dxdthdot=dfdx*dxdth+dfdth;

dxdt=[xdot; dxdthdot(:)];