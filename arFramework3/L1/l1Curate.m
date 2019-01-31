%% DEPRECATED
% l1Curate([jks], parrem, parcomp, (direction), [section], [check], [gradient], [refit])
% 
% L1 curation
% Remove linearly dependent parameters
% 
% jks        [ar.L1jks]
%            indices of the fold-factor parameters to be investigated by L1
%            regularization 
% parrem     parameter to be removed
% parcomp    []
%            parameter to compensate for removed parameter
% direction  define whether parrem and parcomp positively or negatively connected
%            Required, if parcomp is NOT empty
% section    [1:length(ar.L1linv)]
%            indices to be removed
% check      [true]
%            flag indicating whether a check for successful curation should
%            be performed
% gradient   [0]
%            specify gradient in L1
% refit      [false]
%            refit dynamics or just evaluation curation (default: false)

function l1Curate(jks, parrem, parcomp, direction, section, check, gradient, refit)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~isfield(ar,'L1linv') || isempty(ar.L1linv))
    error('Please run L1 scan first.')
end

if(~isfield(ar,'L1pcur'))
    ar.L1pcur = [];
    ar.L1seccur = {};
end

if(~exist('jks','var') || isempty(jks))
    if(~isfield(ar,'L1jks') || isempty(ar.L1jks))
        error('No L1 parameters specified and none present in ar.L1jks.')
    else
        jks = ar.L1jks;
    end
end

if(~exist('parrem','var') || isempty(parrem))
    fprintf('No parameter specified to be removed. Returning.\n')
    return
end

if(~exist('parcomp','var'))
    parcomp = [];
end

if(~exist('direction','var') || isempty(direction))
    if(~isempty(parcomp))
        error('please specify direction for compensation')
    end
end

if(~exist('section','var'))
    section = 1:length(ar.L1linv);
end

if(~exist('check','var') || isempty(check))
    check = true;
end

if(~exist('gradient','var') || isempty(gradient))
    gradient = 0;
end

if(~exist('refit','var') || isempty(refit))
    refit = false;
end

L1ps = ar.L1ps;
L1chi2s = ar.L1chi2s;
L1chi2fits = ar.L1chi2fits;
linv = ar.L1linv;

ar.type(jks) = 3;
ar.qFit(jks) = 1;

if ~isempty(parcomp)
    L1ps(section,parcomp) = L1ps(section,parcomp)+repmat(direction.*L1ps(section,parrem),1,length(parcomp));
end
L1ps(section,parrem) = 0;

for i = section
    for j = 1:length(ar.L1pcur)
        if refit
            if ismember(i,ar.L1seccur{j})
                ar.qFit(ar.L1pcur(j)) = 2;
            end
        end
    end
    ar.std(jks) = linv(i) * (1 + gradient * linspace(0,.001,length(jks)));
    ar.p = L1ps(i,:);
    arChi2
    if refit
        ar.qFit(parrem) = 2;
        arFit
        ar.qFit(jks) = 1;
        arChi2
    end
    L1ps(i,:) = ar.p;
    L1chi2s(i) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
    L1chi2fits(i) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
end

% Check if curation was successful
if check
    
    figure
    plot(L1chi2s-ar.L1chi2s)
    xlabel('Regularization index')
    ylabel('Delta log-likelihood')
    
    if ~refit
        crit = any(abs(L1chi2s-ar.L1chi2s) > 1e-3);
    else
        crit = any(diff(L1chi2fits) < -1e-3 | diff(L1chi2s) < -1e-3);
    end
    if crit 
        fprintf('Check not successful! Curation not accepted ...\n')
    else
        fprintf('Check successful, curation accepted.\n')
        ar.L1ps = L1ps;
        ar.L1chi2s = L1chi2s;
        ar.L1chi2fits = L1chi2fits;
        ar.L1pcur = [ar.L1pcur parrem];
        ar.L1seccur{end+1} = section;
    end
else
    ar.L1ps = L1ps;
    ar.L1chi2s = L1chi2s;
    ar.L1chi2fits = L1chi2fits;
    ar.L1pcur = [ar.L1pcur parrem];
    ar.L1seccur{end+1} = section;
end