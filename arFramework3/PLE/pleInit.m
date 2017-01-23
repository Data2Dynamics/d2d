% Initialize Profile Likelihood Exploit
%
% pleInit(p, q_fit, lb, ub, hb_dist, q_log10, integrate_fkt, objective_fkt, diff_fkt, ...
%         fit_fkt, setoptim_fkt, covar_fkt[, p_labels, stdalpha])
%
% p:                         parameter vector
% q_fit:                     which parameters are fitted
% lb:                        lower parameter bound
% ub:                        upper parameter bound
% q_log10:                   which parameters are logarithmic (basis 10)
% integrate_fkt(p):          set parameters and update objective function
% chi2 = objective_fkt:      get merit value
% [beta, alpha] = diff_fkt:  calculate beta=-0.5*dchi^2/dp and alpha=0.5*d^2chi^2/dp^2
% [p, covar] = fit_fkt(i):   call optimization with i'th parameter fixed
% setoptim_fkt(p):           set persistent optimal parameter values
% p_labels:                  cell array of parameter labels
% stdalpha:                  confidence level in standard deviation [1] = 68%
% 
%   Profile Likelihood Exploit works independently of the D2D framework.
%   This function can be used to apply the method for another tool.
% 
%   pleInit is called by arPLEInit

function pleInit(p, q_fit, lb, ub, q_log10, integrate_fkt, merit_fkt, ...
    diffmerit_fkt, fit_fkt, setoptim_fkt, p_labels, alpha, force)
if(~exist('force','var') || isempty(force))
    force = true;
end

global ar


% fprintf('\nProfile Likelihood Approach for Identifiability Analysis\n\n');
% fprintf('For detailed information please read:\n');
% fprintf('A. Raue, C. Kreutz, T. Maiwald, J. Bachmann, M. Schilling, U. Klingmueller and J. Timmer\n');
% fprintf('Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood.\n');
% 
% fprintf('Bioinformatics (2009), 25(15), 1923-1929.\n');
% site = '"http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btp358v1"';
% title = 'doi:10.1093/bioinformatics/btp358';
% fprintf(['<a href = ' site '>' title '</a>'])
% fprintf('\n\n');

