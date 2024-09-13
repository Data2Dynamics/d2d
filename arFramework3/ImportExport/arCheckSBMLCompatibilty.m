function [qSBMLCompatible, pProblemCond] = arCheckSBMLCompatibilty(m)
% ARCHECHECKSBMLCOMPATIBILTY checks if CONDITIONS in model.def are SBML compatible
%
% BACKGROUND:
% CONDITIONS in model.def correspond to replacements of the original model parameters
% with numerical values of mathematical expressions. In SBML this is best represented as
% "initialAssignment" rules. However, there is an important difference how CONDITIONS
% and initialAssignment rules are applied:
%   - CONDITIONS are applied exactly once and independently of each other
%   - initialAssignment rules are combined to form a set of equations that are solved
%     simultaneously
%
% This means that conditions should be encoded as initialAssignment rules*, but cannot
% be translated directly in all cases. Specifically, there must not be any dependencies
% between the conditions. This means, a parameter that is repalced must not appear in
% any mathematical expression on the RHS of the replacements.
%
% * See arExportSBML_fullmodel.m
%
% EXAMPLE OF DIFFERENT BEHAVIOR:
% consider the expression expr = "p1*x^p2"
%
% d2d:
%   CONDITIONS
%   p1   "2"
%   p2   "p1*2"
%
%   The expression after replacements becomes: "2*x^(2*p1)"
% 
% SBML:
%   <listOfInitialAssignments>
%       <initialAssignment symbol="p1">
%           <math xmlns="http://www.w3.org/1998/Math/MathML">
%               <cn> 2 </cn>
%           </math>
%       </initialAssignment>
%       <initialAssignment symbol="p2">
%           <math xmlns="http://www.w3.org/1998/Math/MathML">
%               <apply>
%                   <times/>
%                   <ci> p1 </ci>
%                   <cn> 2 </cn>
%               </apply>
%           </math>
%       </initialAssignment>
%   </listOfInitialAssignments>
%
%   The expression after replacements becomes: "2*x^(4)"
%   This is because the assignment p1<-2 is also applied within the assignment p2<-p1*2 

arguments
    m (1,1) double {mustBeInteger, mustBePositive} = 1
end


global ar

% identify parameters that are replaced in the conditions
modelP = ar.model(m).p';
modelFPsym = arSym(ar.model(m).fp);
modelFP = cellstr(string(arSym(ar.model(m).fp)));
qReplaced = ~strcmp(modelP, modelFP);

% identify which of the replaced parameters also appear in modelFP expressions
fpParams = unique(cellstr(string(symvar(modelFPsym(qReplaced)))));
qReplacedAndInFP = ismember(modelP, fpParams) & qReplaced;

% exclude init parameters that do not appear in the model equations
qStateInit = ismember(modelP, ar.model(m).px0);
qInModelEquations = ismember(modelP, union(ar.model(m).pu, ar.model(m).pv));
qPureStateInit = qStateInit & ~qInModelEquations;

% exclude compartment size parameters
% Compartment size parameters do not pose a problem, because they are not used in SBML.
% Instead compartment valume is assigned directly to the compartment ID in the SBML file.
qCompSizeParam = ismember(ar.model(m).p', ar.model(m).pc);

% final array of parameters with problematic replacements
qProblemCond = qReplacedAndInFP & ~qCompSizeParam & ~qPureStateInit;
pProblemCond = modelP(qProblemCond);

% Flag that indicates if the model is SBML compatible
qSBMLCompatible = ~any(qProblemCond);