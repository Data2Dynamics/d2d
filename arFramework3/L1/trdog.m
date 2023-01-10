function[s,snod,qpval,posdef,pcgit,Z] = trdog(x,g,H,D,delta,dv,...
    mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l,u,Z,dnewt,preconflag,varargin)

[dofit,out] = determine_dofit(x,g,H,D,delta,dv,...
    mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l,u,Z,dnewt,preconflag,varargin);

s = zeros(size(x));
snod = zeros(size(x));

if(isempty(Z))
    z = Z;
    Z = zeros(length(x),2);
else
    z = Z(dofit,:);
    z = z./(ones(size(z,1),1)*sum(z.^2));
    z(isnan(z)) = 0;
end

ind_L1 = [];
ind_L1DiffPen = [];
ind_grplas = [];
grplas_grps = {};
try
    global ar
    fitted = ar.qFit == 1;
    l1s = find(ar.type(fitted) == 3);
    grplass = find(ar.type(fitted) == 5);
    if ~isempty(l1s)
        ind_L1 = find(ismember(dofit,l1s));
    end
    if ~isempty(grplass)
        grplas_grps = ar.grplas.groupinds;
        ind_grplas = find(ismember(dofit,grplass));
    end
    
    %% CUSTOM
    if isfield(ar, 'L1DiffPen_activate') && ar.L1DiffPen_activate == 1 && ...
            isfield(ar, 'L1DiffPen_useCustomRes') && ar.L1DiffPen_useCustomRes
        
        diffLookup = ar.L1DiffPen_diffs;
        ind_L1DiffPen = diffLookup(sum(ismember(diffLookup,find(~fitted)),2) == 0,:); % remove columns that contain non-fitted parameters;

        % Adjust indices in ind_L1DiffPen by subtracting number of
        % non-fitted parameters with ind smaller than diffs
        
        numNotFitted = [0,cumsum(~fitted)];
        for ii = 1:size(ind_L1DiffPen,1)
            for ij = 1:size(ind_L1DiffPen,2)
                myind = ind_L1DiffPen(ii,ij);
                ind_L1DiffPen(ii,ij) = myind - numNotFitted(myind);
            end
        end
    end
end

if(size(H,1)==size(H,2))
    [s2,snod2,qpval,posdef,pcgit,Z2] = trdog_matlab(x(dofit),g(dofit),H(dofit,dofit),D(dofit,dofit),delta,dv(dofit),...
        mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l(dofit),u(dofit),z,dnewt,preconflag,ind_L1,ind_L1DiffPen,grplas_grps,ind_grplas,varargin{:});
else
    [s2,snod2,qpval,posdef,pcgit,Z2] = trdog_matlab(x(dofit),g(dofit),H(:,dofit),D(dofit,dofit),delta,dv(dofit),...
        mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l(dofit),u(dofit),z,dnewt,preconflag,ind_L1,ind_L1DiffPen,grplas_grps,ind_grplas,varargin{:});
end

s(dofit) = s2;
snod(dofit) = snod2;
Z(dofit,1:size(Z2,2)) = Z2;


function [dofit,out] = determine_dofit(x,g,H,D,delta,dv,...
    mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l,u,Z,dnewt,preconflag,varargin)

thresh = 1e-2;

% Initialization
n = length(g);
grad = D*g;
DM = D;
DG = sparse(1:n,1:n,full(abs(g).*dv));
posdef = 1;
pcgit = 0;
tol2 = sqrt(eps);
v1 = dnewt;
qpval1 = Inf; %#ok  - Kept for clarity
qpval2 = Inf;
qpval3 = Inf;

% DETERMINE A 2-DIMENSIONAL SUBSPACE
if isempty(Z)
    if isempty(v1)
        switch preconflag
            case 'hessprecon'
                % preconditioner based on H, no matter what it is
                [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
            case 'jacobprecon'
                [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
            otherwise
                error(message('optimlib:trdog:InvalidPreconflag'));
        end
        % We now pass kmax in from calling function
        %kmax = max(1,floor(n/2));
        if tol <= 0
            tol = 0.1;
        end
        
        [v1,posdef,pcgit] = pcgr(DM,DG,grad,kmax,tol,...
            mtxmpy,H,R,permR,preconflag,pcoptions,varargin{:});
    end
    if norm(v1) > 0
        v1 = v1/norm(v1);
    end
    Z(:,1) = v1;
    if n > 1
        if (posdef < 1)
            v2 = D*sign(grad);
            if norm(v2) > 0
                v2 = v2/norm(v2);
            end
            v2 = v2 - v1*(v1'*v2);
            nrmv2 = norm(v2);
            if nrmv2 > tol2
                v2 = v2/nrmv2;
                Z(:,2) = v2;
            end
        else
            if norm(grad) > 0
                v2 = grad/norm(grad);
            else
                v2 = grad;
            end
            v2 = v2 - v1*(v1'*v2);
            nrmv2 = norm(v2);
            if nrmv2 > tol2
                v2 = v2/nrmv2;
                Z(:,2) = v2;
            end
        end
    end
end

%  REDUCE TO THE CHOSEN SUBSPACE
W = DM*Z;
switch preconflag
    case 'hessprecon'
        WW = feval(mtxmpy,H,W,varargin{:});
    case 'jacobprecon'
        WW = feval(mtxmpy,H,W,0,varargin{:});
    otherwise
        error(message('optimlib:trdog:InvalidPreconflag'));
end

W = DM*WW;
MM = full(Z'*W + Z'*DG*Z);
rhs = full(Z'*grad);

%  Determine 2-D TR soln
st = trust(rhs,MM,delta);
ss = Z*st;
s = abs(diag(D)).*ss;
s = full(s);
ssave = s;
sssave = ss;
stsave = st;

% Check direction for NaNs
if any(isnan(s))
    error(message('optimlib:trdog:NaNInStep'))
end

% Truncate the TR solution?
arg = (abs(s) > 0);
% No truncation if s is an all zero vector
if ~any(arg)
    alpha = 1;
    mmdis = 1;
    ipt = [];  % Set to empty if all-zero step so that it won't be undefined
else
    dis = max((u(arg)-x(arg))./s(arg), (l(arg)-x(arg))./s(arg));
    [mmdis,ipt] = min(dis);
    mdis = theta*mmdis;
    alpha = min(1,mdis);
end
s = alpha*s;
st = alpha*st;
ss = full(alpha*ss);
qpval1 = rhs'*st + (.5*st)'*MM*st;
if n > 1
    %   Evaluate along the reflected direction?
    qpval3 = Inf;
    ssssave = mmdis*sssave;
    if norm(ssssave) < .9*delta
        r = mmdis*ssave;
        nx = x+r;
        stsave = mmdis*stsave;
        qpval0 = rhs'*stsave + (.5*stsave)'*MM*stsave;
        switch preconflag
            case 'hessprecon'
                ng = feval(mtxmpy,H,r,varargin{:});
            case 'jacobprecon'
                ng = feval(mtxmpy,H,r,0,varargin{:});
            otherwise
                error(message('optimlib:trdog:InvalidPreconflag'));
        end
        
        ng = ng + g;
        ngrad = D*ng;
        ngrad = ngrad + DG*ssssave;
        
        %      nss is the reflected direction
        nss = sssave;
        nss(ipt) = -nss(ipt);
        ZZ(:,1) = nss/norm(nss);
        W = DM*ZZ;
        
        switch preconflag
            case 'hessprecon'
                WW = feval(mtxmpy,H,W,varargin{:});
            case 'jacobprecon'
                WW = feval(mtxmpy,H,W,0,varargin{:});
            otherwise
                error(message('optimlib:trdog:InvalidPreconflag'));
        end
        
        
        W = DM*WW;
        MM = full(ZZ'*W + ZZ'*DG*ZZ);
        nrhs=full(ZZ'*ngrad);
        [nss,tau] = quad1d(nss,ssssave,delta);
        nst = tau/norm(nss);
        ns = abs(diag(D)).*nss;
        ns = full(ns);
        
        % Check direction for NaNs
        if any(isnan(ns))
            error(message('optimlib:trdog:NaNInRflctdStep'))
        end
        
        
    end
    
    %   Evaluate along gradient direction
    gnorm = norm(grad);
    ZZ(:,1) = grad/(gnorm + (gnorm == 0)); % Protect against norm of 0
    W = DM*ZZ;
    
    switch preconflag
        case 'hessprecon'
            WW = feval(mtxmpy,H,W,varargin{:});
        case 'jacobprecon'
            WW = feval(mtxmpy,H,W,0,varargin{:});
        otherwise
            error(message('optimlib:trdog:InvalidPreconflag'));
    end
    
    
    W = DM*WW;
    MM = full(ZZ'*W + ZZ'*DG*ZZ);
    rhs = full(ZZ'*grad);
    st = trust(rhs,MM,delta);
    ssg = ZZ*st;
    sg = abs(diag(D)).*ssg;
    sg = full(sg);
    
    % Check direction for NaNs
    if any(isnan(sg))
        error(message('optimlib:trdog:NaNInGrad'))
    end
    
    
end


out = union(find(x+ssave>u & abs(u-x)<=thresh),find(x+ssave<l & abs(x-l)<=thresh));
if n > 1
    out = union(out,find(x+sg>u & abs(u-x)<=thresh));
    out = union(out,find(x+sg<l & abs(x-l)<=thresh));
end
if exist('ns')
    out = union(out,find(x+(ns + r)>u & abs(u-x)<=thresh));
    out = union(out,find(x+(ns + r)<l & abs(x-l)<=thresh));
end
dofit= setdiff(1:length(x),out);


function[s,snod,qpval,posdef,pcgit,Z] = trdog_matlab(x,g,H,D,delta,dv,...
    mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l,u,Z,dnewt,preconflag,...
    ind_L1,ind_L1DiffPen,grplas_grps,ind_grplas,varargin)
%

%TRDOG Reflected (2-D) trust region trial step (box constraints)
%
% [s,snod,qpval,posdef,pcgit,Z] = TRDOG(x,g,H,D,delta,dv,...
%                 mtxmpy,pcmtx,pcoptions,tol,theta,l,u,Z,dnewt,preconflag);
%
%   Determine the trial step `s', an approx. trust region solution.
%   `s' is chosen as the best of 3 steps: the scaled gradient
%   (truncated to  maintain strict feasibility),
%   a 2-D trust region solution (truncated to remain strictly feas.),
%   and the reflection of the 2-D trust region solution,
%   (truncated to remain strictly feasible).
%
%   The 2-D subspace (defining the trust region problem) is defined
%   by the scaled gradient direction and a CG process (returning
%   either an approximate Newton step of a direction of negative curvature.
%   Driver functions are: SNLS, SFMINBX
%   SNLS actually calls TRDOG with the Jacobian matrix (and a special
%   Jacobian-matrix multiply function in MTXMPY).

%   Copyright 1990-2011 The MathWorks, Inc.

% Initialization
n = length(g);
grad = D*g;
DM = D;
DG = sparse(1:n,1:n,full(abs(g).*dv));
posdef = 1;
pcgit = 0;
tol2 = sqrt(eps);
v1 = dnewt;
qpval1 = Inf; %#ok  - Kept for clarity
qpval2 = Inf;
qpval3 = Inf;

% DETERMINE A 2-DIMENSIONAL SUBSPACE
if isempty(Z)
    if isempty(v1)
        switch preconflag
            case 'hessprecon'
                % preconditioner based on H, no matter what it is
                [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
            case 'jacobprecon'
                [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
            otherwise
                error(message('optimlib:trdog:InvalidPreconflag'));
        end
        % We now pass kmax in from calling function
        %kmax = max(1,floor(n/2));
        if tol <= 0
            tol = 0.1;
        end
        
        [v1,posdef,pcgit] = pcgr(DM,DG,grad,kmax,tol,...
            mtxmpy,H,R,permR,preconflag,pcoptions,varargin{:});
    end
    if norm(v1) > 0
        v1 = v1/norm(v1);
    end
    Z(:,1) = v1;
    if n > 1
        if (posdef < 1)
            v2 = D*sign(grad);
            if norm(v2) > 0
                v2 = v2/norm(v2);
            end
            v2 = v2 - v1*(v1'*v2);
            nrmv2 = norm(v2);
            if nrmv2 > tol2
                v2 = v2/nrmv2;
                Z(:,2) = v2;
            end
        else
            if norm(grad) > 0
                v2 = grad/norm(grad);
            else
                v2 = grad;
            end
            v2 = v2 - v1*(v1'*v2);
            nrmv2 = norm(v2);
            if nrmv2 > tol2
                v2 = v2/nrmv2;
                Z(:,2) = v2;
            end
        end
    end
end

%  REDUCE TO THE CHOSEN SUBSPACE
W = DM*Z;
switch preconflag
    case 'hessprecon'
        WW = feval(mtxmpy,H,W,varargin{:});
    case 'jacobprecon'
        WW = feval(mtxmpy,H,W,0,varargin{:});
    otherwise
        error(message('optimlib:trdog:InvalidPreconflag'));
end

W = DM*WW;
MM = full(Z'*W + Z'*DG*Z);
rhs = full(Z'*grad);

%  Determine 2-D TR soln
st = trust(rhs,MM,delta);
ss = Z*st;
s = abs(diag(D)).*ss;
s = full(s);
ssave = s;
sssave = ss;
stsave = st;

% Check direction for NaNs
if any(isnan(s))
    error(message('optimlib:trdog:NaNInStep'))
end

% Truncate the TR solution?
arg = (abs(s) > 0);
% No truncation if s is an all zero vector
if ~any(arg)
    alpha = 1;
    mmdis = 1;
    ipt = [];  % Set to empty if all-zero step so that it won't be undefined
else
    dis = max((u(arg)-x(arg))./s(arg), (l(arg)-x(arg))./s(arg));
    [mmdis,ipt] = min(dis);
    mdis = theta*mmdis;
    alpha = min(1,mdis);
end
s = alpha*s;
st = alpha*st;
ss = full(alpha*ss);

% if ~isempty(ind_L1)
%     myind = zeros(size(x));
%     myind(ind_L1) = 1;
%     myind = logical(myind);
%     newx = x+s;
%     signchange = (abs(sign(newx) - sign(x)) == 2) & (abs(newx) > 1e-10) & (abs(x) > 1e-10) & (myind);
%     if sum(signchange) > 0
%         distzero = abs(x) ./ abs(s);
%         alpha = min(distzero(signchange));
%         s = alpha*s;
%         st = alpha*st;
%         ss = alpha*ss;
%     end
% end

%%% CUSTOM

if ~isempty(ind_L1)
    myind = zeros(size(x));
    myind(ind_L1) = 1;
    myind = logical(myind);
    newx = x+s;
    signchange = (abs(sign(newx) - sign(x)) == 2) & (abs(newx) > 1e-10) & (abs(x) > 1e-10) & (myind);
    signchangeDiff = 0;
    
    if ~isempty(ind_L1DiffPen)
        xDiffs = bsxfun(@minus, x(ind_L1DiffPen(:,1)), x(ind_L1DiffPen(:,2)));
        sDiffs = bsxfun(@minus, s(ind_L1DiffPen(:,1)), s(ind_L1DiffPen(:,2)));
        xDiffsNew = bsxfun(@minus, newx(ind_L1DiffPen(:,1)), newx(ind_L1DiffPen(:,2)));
        signchangeDiff = (abs(sign(xDiffs) - sign(xDiffsNew)) == 2) & (abs(xDiffsNew) > 1e-10) & (abs(xDiffs) > 1e-10);
    end
    
    if sum(signchange) + sum(signchangeDiff) > 0
        distzero = abs(x) ./ abs(s);
        alpha = min(distzero(signchange));
        if sum(signchangeDiff) > 0
            distzeroDiffs = abs(xDiffs) ./ abs(sDiffs);
            alphaDiff = min(distzeroDiffs(signchangeDiff));
            alpha = min([alpha, alphaDiff]);
        end
        
        s = alpha*s;
        st = alpha*st;
        ss = alpha*ss;
        %fprintf('diff: %d\n',(x(3) - x(5)))
        %fprintf('alpha: %d\n', alpha)
        %fprintf('newdiff: %d\n',(x(3) + s(3) - x(5) - s(5)))
        %fprintf('newdiff: %d\n',(xDiffsNew(5)))
    end
end

if ~isempty(ind_grplas)
    distto0 = nan(1,length(grplas_grps));
    lam = nan(1,length(grplas_grps));
    for g = 1:length(grplas_grps)
        gind = ind_grplas(ismember(ind_grplas,grplas_grps{g}));
        if ~isempty(gind)
            xgind = x(gind);
            sgind = s(gind);
            lam(g) = -dot(xgind,sgind)/dot(sgind,sgind);
            distto0(g) = norm(xgind+ lam(g).* sgind);
        end
    end
    [closestto0,Iclose] = min(distto0);
    if (~isempty(closestto0) && closestto0 < 1e-10 && lam(Iclose) >=0 && lam(Iclose) < 1)
        alpha = lam(Iclose);
        s = alpha*s;
        st = alpha*st;
        ss = alpha*ss;
    end
end

qpval1 = rhs'*st + (.5*st)'*MM*st;
if n > 1
    %   Evaluate along the reflected direction?
    qpval3 = Inf;
    ssssave = mmdis*sssave;
    if norm(ssssave) < .9*delta
        r = mmdis*ssave;
        nx = x+r;
        stsave = mmdis*stsave;
        qpval0 = rhs'*stsave + (.5*stsave)'*MM*stsave;
        switch preconflag
            case 'hessprecon'
                ng = feval(mtxmpy,H,r,varargin{:});
            case 'jacobprecon'
                ng = feval(mtxmpy,H,r,0,varargin{:});
            otherwise
                error(message('optimlib:trdog:InvalidPreconflag'));
        end
        
        ng = ng + g;
        ngrad = D*ng;
        ngrad = ngrad + DG*ssssave;
        
        %      nss is the reflected direction
        nss = sssave;
        nss(ipt) = -nss(ipt);
        ZZ(:,1) = nss/norm(nss);
        W = DM*ZZ;
        
        switch preconflag
            case 'hessprecon'
                WW = feval(mtxmpy,H,W,varargin{:});
            case 'jacobprecon'
                WW = feval(mtxmpy,H,W,0,varargin{:});
            otherwise
                error(message('optimlib:trdog:InvalidPreconflag'));
        end
        
        
        W = DM*WW;
        MM = full(ZZ'*W + ZZ'*DG*ZZ);
        nrhs=full(ZZ'*ngrad);
        [nss,tau] = quad1d(nss,ssssave,delta);
        nst = tau/norm(nss);
        ns = abs(diag(D)).*nss;
        ns = full(ns);
        
        % Check direction for NaNs
        if any(isnan(ns))
            error(message('optimlib:trdog:NaNInRflctdStep'))
        end
        
        % Truncate the reflected direction?
        arg = (abs(ns) > 0);
        % No truncation if s is zero length
        if ~any(arg)
            alpha = 1;
        else
            dis = max((u(arg)-nx(arg))./ns(arg), (l(arg)-nx(arg))./ns(arg));
            mdis = min(dis);
            mdis = theta*mdis;
            alpha = min(1,mdis);
        end
        ns = alpha*ns;
        nst = alpha*nst;
        nss = full(alpha*nss);
        
        %         if ~isempty(ind_L1)
        %             myind = zeros(size(x));
        %             myind(ind_L1) = 1;
        %             myind = logical(myind);
        %             newx = x+(ns + r);
        %             signchange = (abs(sign(newx) - sign(x)) == 2) & (abs(newx) > 1e-10) & (abs(x) > 1e-10) & (myind);
        %             if sum(signchange) > 0
        %                 distzero = abs(x) ./ abs(ns + r);
        %                 alpha = min(distzero(signchange));
        %                 ns = alpha*ns;
        %                 nst = alpha*nst;
        %                 nss = full(alpha*nss);
        %             end
        %         end
        
        %%% CUSTOM
        
        if ~isempty(ind_L1)
            myind = zeros(size(x));
            myind(ind_L1) = 1;
            myind = logical(myind);
            newx = x+(ns+r);
            signchange = (abs(sign(newx) - sign(x)) == 2) & (abs(newx) > 1e-10) & (abs(x) > 1e-10) & (myind);
            signchangeDiff = 0;
            
            if ~isempty(ind_L1DiffPen)
                xDiffs = bsxfun(@minus, x(ind_L1DiffPen(:,1)), x(ind_L1DiffPen(:,2)));
                sDiffs = bsxfun(@minus, ns(ind_L1DiffPen(:,1))+r(ind_L1DiffPen(:,1)), ns(ind_L1DiffPen(:,2))+r(ind_L1DiffPen(:,2)));
                xDiffsNew = bsxfun(@minus, newx(ind_L1DiffPen(:,1)), newx(ind_L1DiffPen(:,2)));
                signchangeDiff = (abs(sign(xDiffs) - sign(xDiffsNew)) == 2) & (abs(xDiffsNew) > 1e-10) & (abs(xDiffs) > 1e-10);
            end
            
            if sum(signchange) + sum(signchangeDiff) > 0
                distzero = abs(x) ./ abs(ns + r);
                alpha = min(distzero(signchange));
                if sum(signchangeDiff) > 0
                    distzeroDiffs = abs(xDiffs) ./ abs(sDiffs);
                    alphaDiff = min(distzeroDiffs(signchangeDiff));
                    alpha = min([alpha, alphaDiff]);
                end
                
%                 distzero = abs(x - r) ./ abs(ns);
%                 alpha = -min(distzero(signchange));
%                 x - r + alpha*ns
%                 r = -r;
                
                ns = alpha*ns;
                r = alpha*r; %????
                nst = alpha*nst;
                nss = alpha*nss;
            end
        end
        
        if ~isempty(ind_grplas)
            distto0 = nan(1,length(grplas_grps));
            lam = nan(1,length(grplas_grps));
            for g = 1:length(grplas_grps)
                gind = ind_grplas(ismember(ind_grplas,grplas_grps{g}));
                if ~isempty(gind)
                    xgind = x(gind);
                    nsrgind = ns(gind)+r(gind);
                    lam(g) = -dot(xgind,nsrgind)/dot(nsrgind,nsrgind);
                    distto0(g) = norm(xgind+ lam(g).* nsrgind);
                end
            end
            [closestto0,Iclose] = min(distto0);
            if (~isempty(closestto0) && closestto0 < 1e-10 && lam(Iclose) >=0 && lam(Iclose) < 1)
                alpha = lam(Iclose);
                ns = alpha*s;
                nst = alpha*st;
                nss = alpha*ss;
            end
        end
        
        qpval3 = qpval0 +  nrhs'*nst + (.5*nst)'*MM*nst;
    end
    
    %   Evaluate along gradient direction
    gnorm = norm(grad);
    ZZ(:,1) = grad/(gnorm + (gnorm == 0)); % Protect against norm of 0
    W = DM*ZZ;
    
    switch preconflag
        case 'hessprecon'
            WW = feval(mtxmpy,H,W,varargin{:});
        case 'jacobprecon'
            WW = feval(mtxmpy,H,W,0,varargin{:});
        otherwise
            error(message('optimlib:trdog:InvalidPreconflag'));
    end
    
    
    W = DM*WW;
    MM = full(ZZ'*W + ZZ'*DG*ZZ);
    rhs = full(ZZ'*grad);
    st = trust(rhs,MM,delta);
    ssg = ZZ*st;
    sg = abs(diag(D)).*ssg;
    sg = full(sg);
    
    % Check direction for NaNs
    if any(isnan(sg))
        error(message('optimlib:trdog:NaNInGrad'))
    end
    
    %   Truncate the gradient direction?
    arg = (abs(sg) > 0);
    if ~any(arg)
        % No truncation if s is zero length
        alpha = 1;
    else
        dis = max((u(arg)-x(arg))./sg(arg), (l(arg)-x(arg))./sg(arg));
        mdis = min(dis);
        mdis = theta*mdis;
        alpha = min(1,mdis);
    end
    sg = alpha*sg;
    st = alpha*st;
    ssg = full(alpha*ssg);
    
    %     if ~isempty(ind_L1)
    %         myind = zeros(size(x));
    %         myind(ind_L1) = 1;
    %         myind = logical(myind);
    %         newx = x+sg;
    %         signchange = (abs(sign(newx) - sign(x)) == 2) & (abs(newx) > 1e-10) & (abs(x) > 1e-10) & (myind);
    %         if sum(signchange) > 0
    %             distzero = abs(x) ./ abs(sg);
    %             alpha = min(distzero(signchange));
    %             sg = alpha*sg;
    %             st = alpha*st;
    %             ssg = alpha*ssg;
    %         end
    %     end
    
    %%% CUSTOM
    
    if ~isempty(ind_L1)
        myind = zeros(size(x));
        myind(ind_L1) = 1;
        myind = logical(myind);
        newx = x+sg;
        signchange = (abs(sign(newx) - sign(x)) == 2) & (abs(newx) > 1e-10) & (abs(x) > 1e-10) & (myind);
        signchangeDiff = 0;
        
        if ~isempty(ind_L1DiffPen)
            xDiffs = bsxfun(@minus, x(ind_L1DiffPen(:,1)), x(ind_L1DiffPen(:,2)));
            sgDiffs = bsxfun(@minus, sg(ind_L1DiffPen(:,1)), sg(ind_L1DiffPen(:,2)));
            xDiffsNew = bsxfun(@minus, newx(ind_L1DiffPen(:,1)), newx(ind_L1DiffPen(:,2)));
            signchangeDiff = (abs(sign(xDiffs) - sign(xDiffsNew)) == 2) & (abs(xDiffsNew) > 1e-10) & (abs(xDiffs) > 1e-10);
        end
        
        if sum(signchange) + sum(signchangeDiff) > 0
            distzero = abs(x) ./ abs(sg);
            alpha = min(distzero(signchange));
            if sum(signchangeDiff) > 0
                distzeroDiffs = abs(xDiffs) ./ abs(sgDiffs);
                alphaDiff = min(distzeroDiffs(signchangeDiff));
                alpha = min([alpha, alphaDiff]);
            end
            
            sg = alpha*sg;
            st = alpha*st;
            ssg = alpha*ssg;
        end
    end
    
    
    if ~isempty(ind_grplas)
        distto0 = nan(1,length(grplas_grps));
        lam = nan(1,length(grplas_grps));
        for g = 1:length(grplas_grps)
            gind = ind_grplas(ismember(ind_grplas,grplas_grps{g}));
            if ~isempty(gind)
                xgind = x(gind);
                sggind = sg(gind);
                lam(g) = -dot(xgind,sggind)/dot(sggind,sggind);
                distto0(g) = norm(xgind+ lam(g).* sggind);
            end
        end
        [closestto0,Iclose] = min(distto0);
        if (~isempty(closestto0) && closestto0 < 1e-10 && lam(Iclose) >=0 && lam(Iclose) < 1)
            alpha = lam(Iclose);
            sg = alpha*s;
            st = alpha*st;
            ssg = alpha*ss;
        end
    end
    
    qpval2 = rhs'*st + (.5*st)'*MM*st;
end

% Choose the best of s, sg, ns.
if qpval2 <= min(qpval1,qpval3)
    qpval = qpval2;
    s = sg;
    snod = ssg;
elseif qpval1 <= min(qpval2,qpval3)
    qpval = qpval1;
    snod = ss;
else
    qpval = qpval3;
    s = ns + r;
    snod = nss + ssssave;
end

%%% CUSTOM DEBUG
% newx = x+s;
% xDiffsNew = bsxfun(@minus, newx(ind_L1DiffPen(:,1)), newx(ind_L1DiffPen(:,2)));
% 
% fprintf('---\n')
% xDiffsNew
% jjj = 1;
% 
%fprintf('Start: %d\n', x(3) - x(5))
%fprintf('End: %d\n', x(3)+s(3) - (x(5)+s(5)))
%xxx = 1; 

%-----------------------------------------------------------
function[nx,tau] = quad1d(x,ss,delta)
%QUAD1D	1D quadratic zero finder for trust region step
%
% [nx,tau] = quad1d(x,ss,delta) tau is min(1,step-to-zero)
% of a 1-D quadratic ay^2 + b*y + c.
% a = x'*x; b = 2*(ss'*x); c = ss'*ss-delta^2). nx is the
% new x value, nx = tau*x;

% Algorithm:
% numer = -(b + sign(b)*sqrt(b^2-4*a*c));
% root1 = numer/(2*a);
% root2 = c/(a*root1);   % because root2*root1 = (c/a);

a = x'*x;
b = 2*(ss'*x);
c = ss'*ss-delta^2;

numer = -(b + sign(b)*sqrt(b^2-4*a*c));
r1 = numer/(2*a);
r2 = c/(a*r1);

tau = max(r1,r2);
tau = min(1,tau);
if tau <= 0
    error(message('optimlib:trdog:SqrRt'));
end
nx = tau*x;



