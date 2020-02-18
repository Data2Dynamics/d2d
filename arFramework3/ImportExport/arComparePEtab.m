function pass = arComparePEtab(ar1,ar2,silent,chi2,p,X,V,Z,Y,conf)

if(nargin==0)
    filenames = fileChooserMulti('./Results', true); 
    if length(filenames)>2
       error('Error: Comparison of more than two models is not supported.') 
    end
    for j=1:length(filenames)
        fname = ['./Results/' filenames{j} '/workspace.mat'];
        if(exist(fname,'file'))
            S=load(fname);
            if j==1
                ar1 = S.ar;
            elseif j==2
                ar2 = S.ar;
            end
        else
            error('Error: No workspace found in %s',fname) 
        end    
    end
end

if ~exist('chi2','var') || isempty(chi2)
    chi2 = true;
end
if ~exist('p','var') || isempty(p)
    p = true;
end
if ~exist('X','var') || isempty(X)
    X = true;
end
if ~exist('V','var') || isempty(V)
    V = true;
end
if ~exist('Z','var') || isempty(Z)
    Z = true;
end
if ~exist('Y','var') || isempty(Y)
    Y = true;
end
if ~exist('conf','var') || isempty(conf)
    conf = true;
end

cnt = 0;
if conf
    b = arCompareConf(ar1,ar2);
    cnt = cnt+b;
end
if V
    b = arCompareV(ar1,ar2,silent);
    cnt = cnt+b;
    if ~b
        warning('arComparePEtab.m: Simulations are not the same.')
    end    
end
if Z
    b = arCompareZ(ar1,ar2,silent);
    cnt = cnt+b;
    if ~b
        warning('arComparePEtab.m: Derived simulations are not the same.')
    end 
end
if X
    b = arCompareX(ar1,ar2,silent);
    cnt = cnt+b;
    if ~b
        warning('arComparePEtab.m: State simulations are not the same.')
    end
end
if Y
    b = arCompareY(ar1,ar2,silent);
    cnt = cnt+b;
    if ~b
        warning('arComparePEtab.m: Observables are not the same.')
    end    
end
if p
    b = arComparePars(ar1,ar2);
    cnt = cnt+b;
    if ~b
        warning('arComparePEtab.m: Parameters are not the same.')
    end
end
if chi2
    b = arCompareChi2(ar1,ar2);
    cnt = cnt+b;
end
cnt = cnt ./ (chi2+p+X+V+Z+Y+conf);

if cnt<1
    pass = 0;
else
    pass = 1;
end