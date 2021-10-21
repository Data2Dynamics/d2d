function [score,averagewidth_alphapars,alphapars] = ...
    score2d(alphapar_step,alphapar_max)
% [score,width_alphapars,alphapars] = score2d(alphapar_step,alphapar_max)
%
% Calculates the score value corresponding to the current 2D-Profile, which
% is obtained by averaging all parameter profile bounds at a specific level
% over the valdiation profile and subsequently averaging all bounds over the
% different parameter confidence levels.
%
% Setting ar.ple2d.config.bounds.vpl_mode = 2 enforces use of smoothed
% 2D-profile.
%
% Inputs:
% alphapar_step: Steps of parameter confidence levels at which parameter 
%       boundaries are calculated and averaged later on
% alphapar_max:  Maximal parameter confidence value which enters the
%       averaging procedure
%
% The optional outputs are the average parameter profile widths width_alphapars
% to the confidence levels in alphapars.
%
% See also: bounds2d

global ar

ar.ple2d.config.bounds.save_mode = 1;
% Do not save temporary bounds results in ar struct

if ar.ple2d.config.bounds.vpl_mode == 2 && ~isfield(ar.ple2d,'smooth')
    disp(['WARNING score2d: Run smooth2d to generate the smooth 2D-landscape ',...
        'or use vpl_mode == 1'])
    return
end
    
if(~exist('alphapar_max', 'var') || isempty(alphapar_max))
    alphapar_max = ar.ple2d.config.bounds.alpha_bounds+10^-4;
end
if(~exist('alphapar_step', 'var') || isempty(alphapar_step))
    alphapar_step = 0.05;
end

% Find all parameter confidence levels
alphapars = alphapar_step;
ii = 1;
while alphapars(ii) < alphapar_max
    alphapars(ii+1) = alphapars(ii) + alphapar_step;
    ii = ii+1;
end
alphapars = alphapars(1:(end-1));

% Find all validation confidence levels and corresponding profiles:
[alpha_pred,pred] = findvalidationpoints;

fprintf('Averaging over %0.3g percent of predictions... \n', ...
    100*max(alpha_pred));

% Assign a weight to each prediction:
weights_vpl = assignweights(alpha_pred);
pred_vpl = pred(~isnan(weights_vpl)); 
weights_vpl = weights_vpl(~isnan(weights_vpl));

averagewidth_alphapars = NaN(1,length(alphapars));
for ii = 1:length(alphapars)
    % Average over all parameter confidence levels:
    [~,~,bounds,par_min] = bounds2d(alphapars(ii));
    % Find parameter bounds at all predictions for a specific parameter level 
    width_tmp = NaN(length(pred_vpl),1);
    try
        for jj = 1:length(pred_vpl)
            % Average over all prediction levels to certain parameter
            % confidence level:
            bounds_tmp = bounds(jj,:);
            bounds_tmp = bounds_tmp(~isnan(bounds_tmp));
            width_tmp(jj) = combinewidth(bounds_tmp,par_min(jj));            
        end
    catch exception
        fprintf(['\n ERROR score2d: Something went wrong while calculating the ',...
         'parameter confidence intervals \n for indices ',...
         'ii=%i, jj=%i \n Error Message: %s \n Line: %s \n'],ii,jj,...
         exception.message,sprintf('%i, ',exception.stack.line));
        return
    end
    averagewidth_alphapars(ii) = weights_vpl*width_tmp;
    % Average alphapars(ii) confidence level width
end

score = sum(averagewidth_alphapars)/length(averagewidth_alphapars);
% Technically this should be a discretized averaging over a uniform
% distribution, but not having the information for large alphas and
% normalizing of the weights give it the form of a simple arithmetic mean.

ar.ple2d.config.bounds.save_mode = 2;

end


function width = combinewidth(intersects,minval)
%Sum width of possibly separated parameter confidence intervals. 

k = 1;
while intersects(k) < minval
    k = k+1;
end
%Profile minimum is between bounds k-1 and k, thus k must be an even number
%or something is wrong.
if rem(k,2) ~= 0
    disp('ERROR score2d: Profile minimum is not between determined profile boundaries.')
    return
end

width_tmp = 0;
while ~isempty(intersects)
    width_tmp = width_tmp + intersects(2)-intersects(1);
    intersects(1:2) = [];
end
width = width_tmp;

end

function weights = assignweights(alphas)
%Given a set of confidence levels (over a set of points) with exactly the
%same values in the case of repetitions, this function assigns a weight to 
%each point.

%Make repeating confidence levels really identical:
[alpha_levels,~,ind_unique_levels] = unique(round(alphas,5));
alphas = alpha_levels(ind_unique_levels);

weights = NaN(1,length(alphas));

if round(min(alphas),1) ~= 0 
    fprintf(['\n WARNING score2d: Validation averaging should comprise ',...
        'the confidence level alpha = 0, but rounded minimal value is %1.3f \n'],...
        min(round(alphas,3)));
end
if alphas(1) ~= alphas(end)
    fprintf(['\n WARNING score2d: Validation averaging reaches higher ',...
        'confidence levels to the one side than to the other. ',...
        '\n This could lead to inadequate results.',...
        'Check how much the levels differ in ar.ple2d.bounds.alpha_pred \n']);
end

for ii = 2:length(alpha_levels)
   weight_number = sum(ind_unique_levels == ii); 
   % Introduces weight loss by counting the same confidence level through 
   % multiple intersections
   if weight_number == 1
       weights(ind_unique_levels == ii) = ...
           (alpha_levels(ii)-alpha_levels(ii-1))/2;
   else
       weights(ind_unique_levels == ii) = ...
           (alpha_levels(ii)-alpha_levels(ii-1))/weight_number;
   end
end

weights = weights/sum(weights,'omitnan'); %normalize sum of weights to one

end


