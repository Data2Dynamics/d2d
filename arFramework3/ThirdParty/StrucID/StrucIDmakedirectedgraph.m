function StrucIDmakedirectedgraph(model,TypeOfAnalysis)
% Creates a directed graph of the model in the App
np=numel(model.parValues);
nu=numel(model.inputNames);
ny=numel(model.outputNames);

if nargin<2
    Observability=true;
    Controllability=false;
else
    Observability=(upper(TypeOfAnalysis)=='O');
    Controllability=~Observability;
end

varValues=cat(1,model.parValues,model.stateValues);
% parValues=varValues(1:np); % not needed because included in varValues
uInput=model.U; %inputValues=InputTable.Data.inputValues;
% compute values for algebraic equations
vAlgebra=model.vAlg;
% vAlgebraValues=vAlgebra(stateValues,parValues,inputValues);

% compute adjacency matrix for directed graph construction
fDyn=model.fDyn;
metafDyn=str2func('xdotTemplate');
hObs=model.Y;
metahObs=str2func('yTemplate');

% NaN's are substituted with random values
index=find(isnan(varValues));
varValues(index)=0.5+rand(size(index)); stateValues=varValues(np+1:end);
A = admDiffComplex(@(x) feval(metafDyn,1,x,uInput,varValues,fDyn,vAlgebra),1,stateValues);
Aadj=double(A~=0);

% Find observable/controllable nodes in the model
if Observability
    C = admDiffComplex(@(x) feval(metahObs,1,x,uInput,varValues,hObs,vAlgebra),1,stateValues);
    Cadj=(C~=0); if ny>1, Cadj=any(Cadj); end
    ObsNodes=find(Cadj);
elseif Controllability
    u0=rand(nu,1);
    C = admDiffComplex(@(u) feval(metafDyn,1,stateValues,u,varValues,fDyn,vAlgebra),1,u0)';
    Cadj=(C~=0); if nu>1, Cadj=any(Cadj); end
    ObsNodes=find(Cadj);
end

w=digraph(Aadj'); figure; H=plot(w);
highlight(H,ObsNodes,'NodeColor','r');