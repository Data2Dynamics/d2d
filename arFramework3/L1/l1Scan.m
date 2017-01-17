% L1 scan
% jks       relative parameters to be investigated by L1 regularization
% linv      width, i.e. inverse slope of L1 penalty (Inf = no penalty; small values = large penalty)
% gradient  use a small gradient on L1 penalty ([-1 0 1]; default = 0)

function l1Scan(jks, linv, gradient)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.type == 3);
    if(isempty(jks))
        error('please initialize by l1Init')
    end
end

if(~exist('linv','var') || isempty(linv))
    linv = logspace(-4,4,49);
    linv = [linv Inf];
    linv = linv(end:-1:1);
end
ar.L1linv = linv;

if(~exist('gradient','var') || isempty(gradient))
    gradient = 0;
end

jks = sort(jks);
optim = ar.config.optimizer;
maxiter = ar.config.optim.MaxIter;

arWaitbar(0);

arFit(true)

ps = nan(length(linv),length(ar.p));
chi2s = nan(1,length(linv));
chi2fits = nan(1,length(linv));

ps(1,:) = ar.p;
chi2s(1) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
chi2fits(1) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
for i = 2:length(linv)
    arWaitbar(i, length(linv), sprintf('L1 scan'));
    ar.std(jks) = linv(i) * (1 + gradient * linspace(0,.001,length(jks)));
    try
        ar.config.optimizer = 1;
        ar.config.optim.MaxIter = 1000;
        arFit(true)
        ar.config.optimizer = 2;
        ar.config.optim.MaxIter = 20;
        arFit(true)
    catch exception
        fprintf('%s\n', exception.message);
    end
    ps(i,:) = ar.p;
    chi2s(i) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
    chi2fits(i) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
    
%     % Backward implementation
%     j = i;
%     if j > 1
%         while chi2fits(j) < max(chi2fits(1:j-1))-1e-3
%             j = j-1;
%             ar.std(jks) = linv(j) * (1 + gradient * linspace(0,.001,length(jks)));
%             try
%                 ar.config.optimizer = 1;
%                 ar.config.optim.MaxIter = 1000;
%                 arFit(true)
%                 ar.config.optimizer = 2;
%                 ar.config.optim.MaxIter = 20;
%                 arFit(true)
%             catch exception
%                 fprintf('%s\n', exception.message);
%             end
%             ps(j,:) = ar.p;
%             chi2s(j) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
%             chi2fits(j) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
%             if j == 1
%                 break
%             end
%         end
%     end
    
    if sum(abs(ps(i,jks)) > 1e-6) == 0;
        ps(i+1:end,:) = repmat(ar.p,size(ps,1)-i,1);
        chi2s(i+1:end) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
        chi2fits(i+1:end) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
        break
    end
end

arWaitbar(-1);

ar.L1ps = ps;
ar.L1chi2s = chi2s;
ar.L1chi2fits = chi2fits;

ar.config.optimizer = optim;
ar.config.optim.MaxIter = maxiter;

