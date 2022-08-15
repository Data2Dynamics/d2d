function [CorParIndex,NumberOfZeroSingularValues,GapSize,dydthAnalysisAllTraj,dudthAnalysisAllTraj, unidentifiable_params] = ...
    StrucIDanalysis(model,TypeOfAnalysis)
% main analysis
% Iden = 1 system is observable/identifiable/reachable
% Iden = 2 system is NOT observable/identifiable/reachable
% Iden = 3 No sensors/inputs available. Change input file to include at least one

% first set default values for arguments that are not included
if nargin<2
    Observability=true; % default is observability analysis
    Controllability=false;
else
    Observability=(upper(TypeOfAnalysis)=='O');
    Controllability=~Observability;
end

Tf = model.Tf;
integrator = model.integrator;
JacobiComplex=model.JacobiComplex; 
nTraject=model.nTraject; % number of trajectories to be included in the analysis

% Set more defaults (can be used at a later stage as in the GUI app)
% for this text version of StrucID
StateReduction=false;
Concatenate=(gt(nTraject,1));

if Concatenate
    fprintf("**Concatenation with %d trajectories! \n",nTraject);
    varNames=cat(1,model.parNames,model.stateNames);
    ind=model.FixedVarIndex;
    fprintf("Parameters that are fixed for ALL trajectories: ");
    for i=1:numel(ind)
        fprintf("%s ",varNames{ind(i)});
    end
    fprintf("\n");
end

% Yout = app.Yout(OutputIndex);
YoutText = model.outputNames; %YoutText(OutputIndex);
ny=numel(YoutText);
nx=numel(model.stateNames); np=numel(model.parNames);
nth=nx+np;
nu=numel(model.inputNames); % number of inputs available

ChoiceAnalysis='1'; %app.TypeOfAnalysisDropDown.Value;
switch ChoiceAnalysis
    case '2' % state reduction
        StateReduction=true;
        kRed=app.InputNumberEditField.Value; % order reduced model kRed<nx
    case '3' % Concatenate
        Concatenate=true;
        %         nTraject=app.InputNumberEditField.Value; % number of trajectories
end

% return if no inputs are present
if (nu==0 && Controllability)
    NumberOfZeroSingularValues=0; GapSize=0; CorParIndex=[];
    return;
end
% return if no outputs are present
if (ny==0 && Observability)
    NumberOfZeroSingularValues=0; GapSize=0; CorParIndex=[];
    return;
end

% parIndex is a ROW vector of indices that indicate the unknown parameters
% in the total set of parameters
UnknownVarIndex=model.UnknownVarIndex;
nthA=numel(UnknownVarIndex); % number of variables (states/parameters) to be analysed
% varValues contains nominal values of both states and parameters
varNames=cat(1,model.parNames,model.stateNames);
thetaNom=cat(1,model.parValues,model.stateValues);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 BEGIN ANALYSIS
% Randomly draw theta from a uniform distribution
% FIX the parameters that have been specified numerically in model file
thetaLow=0.5*ones(nth,1); thetaHigh=1.5*ones(nth,1);
ind=find(isnan(thetaNom)); nNaN=numel(ind);
thetaNom(ind)=thetaLow(ind)+rand(nNaN,1).*(thetaHigh(ind)-thetaLow(ind));

if Concatenate
    thetaReal=zeros(nth,nTraject);
    FixedVarIndex=model.FixedVarIndex;
    for j=1:nTraject
        randomMultipliers=0.9+0.2*rand(nth,1);
        randomMultipliers(FixedVarIndex)=1;
        thetaReal(:,j)=randomMultipliers.*thetaNom;
    end
else
    thetaReal=thetaNom;
end

% integration times
T0=0; % Tf=app.T_finalEditField.Value;
Nt=max(7,nth+1);
timeForward=linspace(T0,Tf,Nt);
timeBackward=linspace(Tf,T0,Nt);

% 6.2 INTEGRATION OF ODEs + SENSITIVITIES
% tol=1e-20; % absolute tolerance on integration
opts = odeset; %('AbsTol',tol);

% function handles for dynamics, IC, input, and algebraic relations
fDyn=model.fDyn; % function handle to xdot=f(x,u,th,v)
modelx0=model.x0; % function handle to initial conditions that may depend on the parameters in the model
vAlgebra=model.vAlg; % function handle to algebraic eqns
uInput=model.U; % function handle to input signal/function

% storage variables for sensitivities in forward and adjoint eqns
dxdthRows=cell(nTraject,1);
dxdthAdjointRows=cell(nTraject,1);
StateTraj=cell(nTraject,1);

% switch true
%     case JacobiComplex
%         disp('Using AdiMATComplex for sensitivity computations...');
%     case ~JacobiComplex
%         disp('Using AdiMATfd for sensitivity computations...');
% end

% disp('Clock started');
tic;

for j=1:nTraject
    if (ny>0)
        theta=thetaReal(:,j);
        x0=modelx0(theta);
        x0Sens=admDiffComplex(@(th) modelx0(th),1,theta);
        x0Aug=cat(1,x0,x0Sens(:));

        if JacobiComplex
                [~,X] = integrator(@(t,x) meta(t,x,uInput,theta,fDyn,...
                    vAlgebra,[nx nth]),timeForward,x0Aug,opts);
        else
                [~,X] = integrator(@(t,x) metaFD(t,x,uInput,theta,fDyn,...
                    vAlgebra,[nx nth]),timeForward,x0Aug,opts);
        end

        dxdthRows{j}=X(:,(nx+1):end);
        StateTraj{j}=X(:,1:nx);
    end

    if (nu>0)
        xf=X(end,1:nx)'; % start backwards integration with final state
        theta(np+1:end)=xf;
        %         x0Sens=x0Sens(:,(np+1):end); % only states are relevant
        xfSens=admDiffComplex(@(th) modelx0(th),1,theta);
        xfSens=xfSens(:,np+1:end); % only states are considered for the controllability problem
        %         x0Aug=cat(1,x0,x0Sens(:));
        xfAug=cat(1,xf,xfSens(:));

        if JacobiComplex
                [~,X] = integrator(@(t,x) metaAdjoint(t,x,uInput,theta,fDyn,...
                    vAlgebra,[nx nth]),timeBackward,xfAug,opts);
        else
                [~,X] = integrator(@(t,x) metaAdjointFD(t,x,uInput,theta,fDyn,...
                    vAlgebra,[nx nth]),timeBackward,xfAug,opts);
        end

        dxdthAdjointRows{j}=X(:,(nx+1):end);
    end
end

% Analysis of sensitivities (from dxdth to dydth)
% define observation function
hObs=model.Y; 

% Build the output sensitivity matrix
dydthAnalysisAllTraj=[];
dudthAnalysisAllTraj=[];
for k=1:nTraject
    dydthAnalysis=[];
    dudthAnalysis=[];

    theta=thetaReal(:,k);

    for j=1:Nt
        if (ny>0) % output signals are available
                X=StateTraj{k}(j,:);
                dhdx=admDiffComplex(@(x) yTemplate(timeForward(j),x,uInput,theta,hObs,vAlgebra),1,X);
                dxdth=reshape(dxdthRows{k}(j,:),nx,nth);
                dhdth = admDiffComplex(@(th) yTemplate(timeForward(j),X,uInput,th,hObs,vAlgebra),1,theta);
                dydth = dhdx*dxdth + dhdth;
                dydthAnalysis=cat(1,dydthAnalysis,dydth(:,UnknownVarIndex));
        end

        if (nu>0) % input signals are available
            X=StateTraj{k}(j,:); u0=uInput(timeForward(j),theta);
            dfdu=admDiffComplex(@(u) xdotTemplate(timeForward(j),X,u,theta,fDyn,vAlgebra),1,u0);
            dxdth=reshape(dxdthRows{k}(j,:),nx,nth);
            dudth=dfdu'*dxdth;
            dudthAnalysis=cat(1,dudthAnalysis,dudth);
        end
    end
    % append sensitivities for each trajectory
    if (ny>0)
        dydthAnalysisAllTraj=cat(1,dydthAnalysisAllTraj,dydthAnalysis);
    end
    if (nu>0)
        dudthAnalysisAllTraj=cat(1,dudthAnalysisAllTraj,dudthAnalysis);
    end
end

switch ChoiceAnalysis
    case {'1','3'}
        if Observability
            [~,S,V]=svd(dydthAnalysisAllTraj,'econ');
        elseif Controllability
            [~,S,V]=svd(dudthAnalysisAllTraj,'econ');
        end
    case '2' % 'balanced' reduction
        %             if (ny>0)
        %                 [~,S,V2]=svd(dydthAnalysisAllTraj);
        %             else
        %                 V2=eye(nx);
        %             end
        %
        %             if (nu>0)
        %                 [~,~,V1]=svd(dudthAnalysisAllTraj);
        %             else
        %                 V1=eye(nx);
        %             end
        %             V=V1*V2';
        [~,S,V]=svd(cat(1,dydthAnalysisAllTraj,dudthAnalysisAllTraj),'econ');
        % [~,S,V]=svd(dudthAnalysisAllTraj,'econ');
end

% Process SVD results
SingularValues=diag(S)';
rankdydth=sum(SingularValues>(max(size(dydthAnalysis))*eps(S(1)))); % from Matlab rank documentation
NumberOfZeroSingularValues=nthA-rankdydth;

logSingularValues=log10(SingularValues);
logSingularValues(isinf(logSingularValues))=log10(eps); % for hard zeros
diffLogSingVal=abs(diff(logSingularValues));

% ThresholdGapSize=0; % adjust this criterion if needed
IndexSingularValuesAfterGap=nthA-NumberOfZeroSingularValues+1;

GapSize=max(diffLogSingVal);

% if (rankdydth<nthA) && (GapSize>ThresholdGapSize)
%     Iden=2;
% else
%     Iden=1;
% end

nullspace=V(:,IndexSingularValuesAfterGap:end);

if (NumberOfZeroSingularValues>0)
    CorPar=logical(abs(nullspace')>2e-2);
else
    CorPar=[]; % logical(abs(V(:,end))'>2e-2);
end

[nr,~]=size(CorPar);
if nr>1
    CorPar=any(CorPar);
end

% store symmetries when probing a network
if (NumberOfZeroSingularValues==0)
    fprintf('%s\n','No correlated parameters found');
else
    fprintf('Unidentifiable Parameters: ');
    CorParIndex=UnknownVarIndex(CorPar);
    unidentifiable_params = strings(1, numel(CorParIndex));
    for k=1:numel(CorParIndex)
        fprintf('%s  ',varNames{CorParIndex(k)});
        unidentifiable_params(k) = varNames{CorParIndex(k)};
    end
    fprintf('\n\n');
end

toc;

if StateReduction
    SingValTraj=diag(S);

    % push results up into the Workspace
    assignin('base','dydth',dydthAnalysis);
    assignin('base','TransformationMatrix',V);
    assignin('base','Sdiag',SingValTraj);
    %     [~,S,V]=svd(cat(1,dydthAnalysis,dudthAnalysis));
    % Simulate the REDUCED model
    theta=thetaNom;

    yRed=[]; yFull=[];

    x0=modelx0(theta); % in case IC depends on theta
    z0=V(:,1:kRed)'*x0;

    NtSim=400;
    timeForward=linspace(0,Tf,NtSim);
    [~,X] = integrator(@(t,x) xdotTemplate(t,x,uInput,theta,fDyn,vAlgebra),timeForward,x0,opts);
    [~,Z] = integrator(@(t,z) zdotReduceTemplate(t,z,uInput,theta,fDyn,vAlgebra,V,kRed),timeForward,z0,opts);

    yRedStep=zeros(NtSim,ny); yFullStep=zeros(NtSim,ny);
    for j=1:NtSim
        yRedStep(j,:)=yReduceTemplate(timeForward(j),Z(j,:)',uInput,theta,hObs,vAlgebra,V,kRed)';
        yFullStep(j,:)=yTemplate(timeForward,X(j,:)',uInput,theta,hObs,vAlgebra);
    end
    yRed=cat(1,yRed,yRedStep);
    yFull=cat(1,yFull,yFullStep);

    figure(2); plot(timeForward,yFull,timeForward,yRed,'ro')%,'MarkerSize',3);
    assignin('base','time',timeForward);
    assignin('base','yRed',yRed);
    assignin('base','yFull',yFull);
end

% 8.5 identifiability graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logSingularValues(isinf(logSingularValues))=log10(eps);

switch ChoiceAnalysis
    case {'1','3'}
        figure(1);
        subplot(2,1,1);
        plot(logSingularValues,'.','MarkerSize',28);
        xlabel('Singular values in descending order');
        ylabel('Log_{10}(Singular Values)');
        title('Log_{10}(Singular Values)','FontWeight','bold');
        axis([0 nthA+1 -Inf Inf]);
        xticks(0:nthA+1);
        minS=min(SingularValues(:,end));
        text(1,2+min(logSingularValues),['\sigma_{',num2str(nthA),'}=',num2str(minS,'%10.5e\n')],'Fontname','FixedWidth');

        subplot(2,1,2);
        if NumberOfZeroSingularValues>=1
            stem(V(:,IndexSingularValuesAfterGap:end));
        else
            stem(V(:,end));
        end
        axis([0 nthA+1 -1 1])
        xticks(0:nthA+1);
        labels=cell(nthA+2,1);

        if nthA<15
            for i=1:nthA
                labels{i+1}=varNames{UnknownVarIndex(i)};
            end
            xticklabels(labels);
        end

        xlabel('Parameters');
        ylabel('Components in Singular Vector(s)');
        if eq(NumberOfZeroSingularValues,0)
            title('Last Column of V','FontWeight','bold');
        elseif eq(NumberOfZeroSingularValues,1)
            title('Last Column of V','FontWeight','bold');
        else
            title(['Last ', num2str(NumberOfZeroSingularValues) ,' Columns of V'],'FontWeight','bold');
        end
end
