% scan likelihood
%
% arSample(jks, N)
% 


function arScan(jks, N)

global ar

if(~exist('N','var'))
    N = 100;
end
if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qFit==1);
end

pReset = ar.p;

ar.scan.ps = {};
ar.scan.chi2s = {};
ar.scan.constrs = {};

ccount = 1;
tic;
arWaitbar(0);
for jk=jks
    ar.p = pReset;
    ps = linspace(ar.lb(jk), ar.ub(jk), N);
    ps = sort([ps ar.p(jk)]);
    
    ar.scan.ps{jk} = ps;
    ar.scan.chi2s{jk} = nan(size(ps));
    ar.scan.constrs{jk} = nan(size(ps));
    for j=1:length(ps);
        arWaitbar(j+((ccount-1)*(N+1)), (N+1)*length(jks), ...
            sprintf('likelihood scan for %s', strrep(ar.pLabel{jk},'_','\_')));
        ar.p(jk) = ps(j);
        try
            arCalcMerit(false, []);
            ar.scan.chi2s{jk}(j) = arGetMerit;
%             ar.scan.constrs{jk}(j) = ar.chi2constr;
        catch error_id
            fprintf('%s for %g: %s\n', ar.pLabel{jk}, ps(j), error_id.message);
        end 
    end
    ccount = ccount + 1;
end
arWaitbar(-1);
fprintf('mean evaluation time %f sec\n', toc/((N+1)*length(jks)));

ar.p = pReset;
arCalcMerit(false, []);


