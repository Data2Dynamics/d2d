% chi2 = pleManual(whichp,pvals,nLHS)
% 
%   The parameter specified by whichp is fixed to values given by pvals.
%   Then, all other parameters are fitted. This function provides
%   the possibility to calculate a paraemters profile likelihood in a
%   user-defined range and step-size.
% 
%   The resuls are stored in ar.pleManual
%   
%   whichp  parameter index or parameter name (according to ar.pLabel)
%           Example: whichp = 12; 
%           Example: whichp = 'k_on'; 
% 
%   pvals   a vector of values for the specified paraemeter, 
%           Example: pvals = logspace(-3,0,20);
%
%   nLHS    number of LHS fits to do for each point, default number is 1
% 
%   chi2    the merit function after fitting

function chi2 = pleManual(whichp,pvals,nLHS)
if nargin == 2   % if the number of inputs equals 2
  nLHS = 1; % then make the third value, z, equal to my default value, 5.
end
if ischar(whichp)
    jp = strmatch(whichp,ar.pLabel,'exact');
    if length(jp)~=1
        error('Parameter %s not found or found multiple times in ar.pLabel',whichp);
    end
elseif isnumeric(whichp)
    jp = whichp;
else
    error('which p has wrong format, has to be index or parameter name.');
end

global ar


pin = ar.p+0;
qFit = ar.qFit+0;

ar.pleManual.ps{jp} = NaN(length(pvals),length(ar.p));
ar.pleManual.chi2s{jp} = NaN(length(pvals),1);
try
    for i=1:length(pvals)
        ar.qFit(jp) = 0;
        ar.p(jp) = pvals(i);
        try %continue if arFit produces an error
            if nLHS == 1
                arFit;
            else
                arFitLHS(nLHS)
            end
        catch
            disp('ERROR using arFit or arFitLHS');
        end
        ar.pleManual.ps{jp}(i,:) = ar.p;
        ar.pleManual.chi2s{jp}(i) = arGetMerit;
        disp(join(['    ',string(i),'/',string(length(pvals))]));
    end
catch ERR
    ar.p = pin;
    ar.qFit = qFit;
    rethrow(ERR)
end

ar.p = pin;
ar.qFit = qFit;

chi2 = ar.pleManual.chi2s{jp};
