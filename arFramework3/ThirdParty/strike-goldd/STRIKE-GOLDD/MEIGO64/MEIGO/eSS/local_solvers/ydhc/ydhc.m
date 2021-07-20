function [x, fval, exitflag, output] = ydhc(fun,x0,lb,ub,options)
% Optimization via Dynamic Hill Climbing.
%
% An implementation with slight modifications of the algorithm as described
% in [De La Maza and Yuret. Dynamic Hill Climbing].
%
% Input:
% fun     : objective function to be minimized
% x0      : initial guess for parameters
% lb, ub  : bounds for parameters, i.e. lb <= x <= ub; no value should be 
%           inf and the difference ub-lb in no component be < 0
% options : struct with options for the algorithm:
%   TolX              : tolerance of parameter
%   TolFun            : tolerance of objective function
%   MaxFunEvals       : maximum number of evaluations of fun
%   OutputFcn         : for visual output after each iteration
%   InitialStepSize   : rel. to 1 (default 0.1)
%   ExpandFactor      : factor of expansion upon success (default 2)
%   ContractFactor    : factor of contraction upon failure (default 0.45)
%   StuckSearchFactor : how far to expand again after got stuck (default 4)
%   Barrier           : use barrier on bounds (default extreme barrier)
%   Display           : off|iter|debug, text output (default off)
%
% Output:
% x   : best guess for parameters
% fval: objective function at the solution, generally fval=fun(x)
% exitflag:
%   1 : The function converged to a solution x
%   0 : Number of function evaluations exceeded options.MaxFunEvals.
%   -1: The algorithm was terminated inappropriately
% output : struct with meta information:
%   iterations  : number of iterations
%   funcCount   : number of function evaluations
%   algorithm   : name of the algorithm
%   t_cpu       : cpu time
%
% History:
% 2017/09/27 Yannik Schaelte

% number of variables
dim  = length(x0);

% interpret options
options  = f_get_options(options);
% extract often used options
tolX                = options.TolX;
tolFun              = options.TolFun;
maxFunEvals         = options.MaxFunEvals;
initialStepSize     = options.InitialStepSize;
expandFactor        = options.ExpandFactor;
contractFactor      = options.ContractFactor;
stuckSearchFactor   = options.StuckSearchFactor;
barrier             = options.Barrier;

% create column vectors
lb      = lb(:);
ub      = ub(:);
x0      = x0(:);

% normalize to unit square
normalize   = @(x) f_normalize(x,lb,ub);
denormalize = @(y) f_denormalize(y,lb,ub);
y0      = normalize(x0);
tolY    = tolX / norm(ub-lb);

% output function for in-time analysis
outputFcn = @(y,fval,funEvals,mode) f_output(y,denormalize,fval,funEvals,mode,options.OutputFcn);

% max step size
vmax = initialStepSize * ones(dim,1);

% wrap function to consider boundaries
fun = @(y,funEvals) f_wrap_fun(denormalize(y),fun,lb,ub,barrier,funEvals,maxFunEvals);

% init run variables
smax   = [vmax;-1];        % array of max step sizes, extra value for extra vector
[step,norms]   = f_init_step(vmax); % matrix of step vectors
xstep  = step;             % steps before last motion
xnorms = norms;
gradv  = zeros(dim,1);     % gradient vector
gradi  = -1;               % index of gradient vector, -1 indicates gradv is not set
prevj  = -1;               % index of last step taken
stuck  = false;            % is process stuck (in min/max)? set when step sizes are small
% then step vectors are increased
done   = false;            % is some finishing criterion fulfilled?

nVec   = 2*dim + 2;            % maximum index in step matrix
opp_j  = @(j) f_opp_j(j,nVec); % short for opposite index in step matrix

% init meta variables
exitflag  = -1;        % flag indicating exit reason
starttime = cputime;   % to measure time difference
funEvals  = 0;         % function evaluations, should be <= maxFunEvals

% init x, fval
ybst      = y0;
fbst      = fun(ybst,funEvals);
funEvals  = funEvals + 1;

outputFcn(ybst,fbst,funEvals,'init'); % create new figure and initialize
outputFcn(ybst,fbst,funEvals,'iter'); % first iteration with start point

% main loop
while ~done
    
    if stuck
        % choose the smallest step, if any is smaller than the maximum size
        [v,j] = f_min(step,norms,smax);
    else
        % choose the largest step
        [v,j] = f_max(step,norms);
    end
    
    % textual output
    f_display(options.Display,funEvals,fbst,norm(v));
    
    % j == -1 indicates minimum found
    if j ~= -1 && funEvals <= maxFunEvals
        % compute next x, fval
        ycur = ybst + v;
        fcur = fun(ycur,funEvals);
        funEvals = funEvals + 1;
        
        delta_f = fcur - fbst;
        
        % is better estimate?
        if delta_f < 0
            ybst = ycur;
            fbst = fcur;
        end
        
        % is significantly better estimate?
        if delta_f < 0 && (norm(v)>tolX/stuckSearchFactor || abs(delta_f) > tolFun)
            % we are not stuck somewhere (anymore)
            stuck = false; 
            % contract opp step to not try the previous point next
            step(:,opp_j(j)) = -contractFactor*v;
            norms(opp_j(j)) = norm(step(:,opp_j(j)));
            % if last step repeated, expand the current step
            if j == prevj
                v = expandFactor*v;
            end
            % xstep always contains the steps of the last time we moved
            xstep = step;
            xnorms = norms;
            % record the last step
            prevj = j;
            % record step
            step(:,j) = v;
            norms(j) = norm(v);
            
            % update gradient vector
            if gradi == -1
                % if gradv empty, set gradv to current step and record index
                gradv = v;
                gradi = min([j, opp_j(j)]);
            elseif gradi == min([j, opp_j(j)])
                % if gradv is parallel to current step, add the current vector
                gradv = (gradv + v);
            else
                % else update the extra vector
                step(:,dim+1) = gradv + v;
                norms(dim+1) = norm(gradv + v);
                % set extra entry in smax to max of current step smax and gradv smax
                % (bound for norm of gradient)
                smax(dim+1)   = max([smax(min([j,opp_j(j)])),smax(gradi)]);
                % set the opp step to - extra step
                step(:,dim+2) = -contractFactor*step(:,dim+1);
                norms(dim+2) = abs(contractFactor)*norms(dim+1);
                % update gradv
                gradi =-1;%        = min([j, opp_j(j)]);
            end
        else % not significantly better estimate
            if stuck
                % if already stuck, increase the current step size
                step(:,j) = expandFactor*v;
                norms(j) = norm(step(:,j));
            elseif norm(step(:,j)) >= tolY
                % if current step norm >= tolY, decrease the current step size
                step(:,j) = contractFactor*v;
                norms(j) = norm(step(:,j));
            else
                % else (norm < tolY): set the stuck flag and set all steps to
                % expandFactor times the last recorded steps
                stuck = true;
                step = expandFactor * xstep;
                norms = abs(expandFactor)*xnorms;
            end
        end
        
        % update output
        outputFcn(ybst,fbst,funEvals,'iter');
        
    else % somehow done
        done = true;
        if funEvals <= maxFunEvals
            % maybe found a local minimum
            exitflag = 1;
        else
            % needed too long
            exitflag = 0;
        end
    end
    
end

% finalize output
outputFcn(ybst,fbst,funEvals,'done');
% textual output
f_display(options.Display,funEvals,fbst,norm(v),true)

% assign return values
x                   = denormalize(ybst);
fval                = fbst;
output = struct();
output.funcCount    = funEvals;
output.iterations   = funEvals;
output.algorithm    = 'Dynamic Hill Climb';
output.t_cpu        = cputime - starttime;

end % function


%% helper functions


function y = f_normalize(x,lb,ub)
% normalize vector to [0,1]

if any(~isfinite(lb)) || any(~isfinite(ub))
    y = x;
else
    y = (x-lb)./abs(ub-lb);
end

end


function x = f_denormalize(y,lb,ub)
% denormalize vector from [0,1]

if any(~isfinite(lb)) || any(~isfinite(ub))
    x = y;
else
    x = y.*(ub-lb) + lb;
end

end


function j_opp = f_opp_j(j,nVec)
% short for opposite vector in step matrix

j_opp = nVec - (j-1);

end


function fval = f_wrap_fun(x,fun,lb,ub,barrier,funEvals,maxFunEvals)
% wrap around function to allow for a barrier function wrap

% set fun to inf whenever conditions not fulfilled
if ~isequal(barrier,'')
    fval = fun(x);
    fval = barrierFunction(fval, [], x, [lb, ub], funEvals, maxFunEvals, barrier);
else
    % extreme barrier
    if any(x>ub) || any(x<lb)
        fval = inf;
    else
        fval = fun(x);
    end
end

end


function [step,norms] = f_init_step(vmax)
% create dim x (2*dim+2)-matrix containing the step proposals as columns
% also return the norms of the steps to reduce computations

dim  = length(vmax);

nVec = 2*dim + 2;
step = zeros(dim,nVec);
norms = zeros(1,nVec);
% 2 positions in the center reserved for gradient
for j = 1:dim
    initStepSize       = abs(vmax(j));
    step(j,j)          = initStepSize;
    step(j,nVec-(j-1)) = -initStepSize;
    norms(j)           = initStepSize;
    norms(nVec-(j-1))  = initStepSize;
end

end


function [v_max,j_max] = f_max(step,norms)
% find step with maximal norm

[~,j_max] = max(norms);
v_max = step(:,j_max);

end


function [v_min,j_min] = f_min(step,norms,smax)
% find step with smallest norm, if one exists whose norm is smaller than
% the corresponding smax max step size

minnorm = -1;
j_min   = -1; % -1 used as indicator in calling function
v_min   = -1;
nVec    = size(step,2);
for j=1:nVec
    vnorm = norms(j);
    if ( vnorm <= smax(min([j, nVec-(j-1)])) && (vnorm < minnorm || minnorm < 0) && vnorm > 0) % tolerance % smax(min([j, nVec-(j-1)]))
        minnorm = vnorm;
        j_min = j;
    end
end

if (j_min ~= -1)
    v_min = step(:,j_min);
end

end


function [ options ] = f_get_options(options_in)
% fill non-existent fields with default values, and check validity

options = struct();
options.TolX                = 1e-8;
options.TolFun              = 1e-8;
options.MaxFunEvals         = Inf;
options.OutputFcn           = nan;
options.InitialStepSize     = 0.1;
options.ExpandFactor        = 2.1;
options.ContractFactor      = 0.47;
options.StuckSearchFactor   = 4;
options.Barrier             = '';
options.Display             = 'off';

% fill from input
cell_fieldnames = fieldnames(options);
cell_fieldnames_in = fieldnames(options_in);

for jf = 1:length(cell_fieldnames_in)
    fieldname = cell_fieldnames_in{jf};
    if ~any(strcmp(cell_fieldnames,fieldname))
        error(['Options field ' fieldname ' does not exist.']);
    end
    options.(fieldname) = options_in.(fieldname);
end

end


function f_output(y,f_denormalize,fval,funEvals,state,outputFcn)
% short for call to output function
% state: 'init', 'iter', or 'done'

if isa(outputFcn,'function_handle')
    x = f_denormalize(y);
    optimValues.fval = fval;
    optimValues.iteration = funEvals;
    outputFcn(x,optimValues,state);
end

end


function f_display(display,funEvals,fbst,vnorm,final)
% short for call to display on screen

if nargin < 5, final = false; end

if strcmp(display,'iter') || strcmp(display,'debug')
    if (strcmp(display,'debug'))
        show_output = true;
    else
        show_output = mod(funEvals,100) == 1;
    end
    
    if show_output && ~final
        if mod(funEvals,1000) == 1
            fprintf('fevals\t|\tfbst\t|\tstepnorm\n');
        end
        fprintf(strcat('%d\t|\t%.8e\t|\t%.8e\n'),funEvals,fbst,vnorm);
    end
    
    if final
        fprintf('final: \t funEvals: %d, \t fbst: %.8e\n',funEvals,fbst); 
    end

end

end
