function arPLErolling(jk, n)

global ar

if(~exist('jk','var') || isempty(jk))
    jk = find(ar.qFit==1);
end
if(length(jk)>1)
    for j=1:length(jk)
        arPLErolling(jk(j));
    end
    return;
end
if(~exist('n','var'))
    n = 100;
end

if(~isfield(ar, 'ple'))
    ar.ple.chi2s = {};
    ar.ple.errors = {};
    ar.ple.lambdas = {};
    ar.ple.ps = {};
    ar.ple.run = zeros(size(ar.p));
    ar.ple.alpha = 0.05;
    ar.ple.ndof = 1;
end

arCalcMerit(true);

chi2Reset = ar.chi2fit;
pReset = ar.p;
ar.ple.pStart = ar.p;

fprintf('PLE #%i for %s...\n', jk, ar.pLabel{jk});
arWaitbar(0);
tic;

[chi2sup, psup, errorsup, lambdasup] = ple_task(jk, n, 1, chi2Reset, pReset);
[chi2sdown, psdown, errorsdown, lambdasdown] = ple_task(jk, n, -1, chi2Reset, pReset);

ar.ple.chi2s{jk} = [fliplr(chi2sdown) chi2Reset chi2sup];
ar.ple.ps{jk} = [flipud(psdown); pReset; psup];
ar.ple.errors{jk} = [fliplr(errorsdown) nan errorsup];
ar.ple.lambdas{jk} = [fliplr(lambdasdown); nan(size(ar.p)); lambdasup];
ar.ple.run(jk) = 1;

arWaitbar(-1);
fprintf('PLE #%i %s elapse time\n', jk, secToHMS(toc));

ar.p = pReset;
arCalcMerit(false);


function [chi2s, ps, errors, lambdas] = ple_task(jk, n, direction, chi2Reset, pReset)

global ar

chi2s = nan(1,n);
ps = nan(n,length(pReset));
lambdas = nan(n,length(pReset));
errors = nan(1,n);

dchi2 = chi2inv(1-ar.ple.alpha, ar.ple.ndof);

dp = 1;
addres = zeros(1,length(pReset));
addres(jk) = -dp*direction;

ar.p = pReset;
pLast = pReset;

quadprog_optims = optimset('Display', 'off');
warning('off','optim:lsqlin:LinConstraints');
warning('off','optim:quadprog:SwitchToMedScale');

for j=1:n
    if(direction>0)
        arWaitbar(j,2*n, sprintf('PLE up for %s', strrep(ar.pLabel{jk},'_','\_')));
    else
        arWaitbar(j+n,2*n, sprintf('PLE down for %s', strrep(ar.pLabel{jk},'_','\_')));
    end
   
%     deltap = 1;
%     while(norm(deltap)>1e-3)
        arCalcMerit(true);
        
        qFit = ar.qFit==1;
        
        res = ar.res;
        sres = ar.sres;
        
        beta = -ar.res*ar.sres;
        beta(jk)
        
        res = [res addres(qFit)]; %#ok<AGROW>
        sres = [sres(:,qFit);eye(sum(qFit))]; %#ok<AGROW>
        
        [deltap,~,~,~,~,lambda] = lsqlin(-sres, res, [], [], [], [], ...
            ar.lb(qFit) - pLast(qFit), ar.ub(qFit) - pLast(qFit), [], quadprog_optims);
        
        %     disp([sum(lambda.lower~=0) sum(lambda.upper~=0)]);
        
        pLast(qFit) = pLast(qFit) + deltap';
        ar.p = pLast;
%     end
    
    ps(j,:) = ar.p;
    chi2s(j) = ar.chi2fit;
    lambdas(j,:) = lambda.upper~=0 | lambda.lower~=0;
        
    if(ar.chi2fit - chi2Reset > 2*dchi2)
        fprintf('PLE #%i reached confidence limit\n', jk); 
        break;
    end
    
    if(lambda.upper(jk)~=0)
        fprintf('PLE #%i reached upper parameter bound\n', jk);
        break;
    end
    if(lambda.lower(jk)~=0)
        fprintf('PLE #%i reached lower parameter bound\n', jk);
        break;
    end
    
%     if(norm(deltap) <= 1e-3)
%         dp = dp * 1.1;
%         addres(jk) = -dp*direction;
%     end
end

ar.p = pReset;



% function [chi2s, ps, errors, lambdas] = ple_task(jk, n, direction, chi2Reset, pReset)
% 
% global ar
% 
% chi2s = nan(1,n);
% ps = nan(n,length(pReset));
% lambdas = nan(n,length(pReset));
% errors = nan(1,n);
% 
% dchi2 = chi2inv(1-ar.ple.alpha, ar.ple.ndof);
% 
% dp = 1;
% addres = zeros(1,length(pReset));
% addres(jk) = -dp*direction;
% 
% ar.p = pReset;
% pLast = pReset;
% 
% quadprog_optims = optimset('Display', 'off');
% warning('off','optim:lsqlin:LinConstraints');
% warning('off','optim:quadprog:SwitchToMedScale');
% 
% for j=1:n
%     if(direction>0)
%         arWaitbar(j,2*n, sprintf('PLE up for %s', strrep(ar.pLabel{jk},'_','\_')));
%     else
%         arWaitbar(j+n,2*n, sprintf('PLE down for %s', strrep(ar.pLabel{jk},'_','\_')));
%     end
%    
% %     deltap = 1;
% %     while(norm(deltap)>1e-3)
%         arCalcMerit(true);
%         
%         qFit = ar.qFit==1;
%         
%         res = ar.res;
%         sres = ar.sres;
%         
%         beta = -ar.res*ar.sres;
%         beta(jk)
%         
%         res = [res addres(qFit)]; %#ok<AGROW>
%         sres = [sres(:,qFit);eye(sum(qFit))]; %#ok<AGROW>
%         
%         [deltap,~,~,~,~,lambda] = lsqlin(-sres, res, [], [], [], [], ...
%             ar.lb(qFit) - pLast(qFit), ar.ub(qFit) - pLast(qFit), [], quadprog_optims);
%         
%         %     disp([sum(lambda.lower~=0) sum(lambda.upper~=0)]);
%         
%         pLast(qFit) = pLast(qFit) + deltap';
%         ar.p = pLast;
% %     end
%     
%     ps(j,:) = ar.p;
%     chi2s(j) = ar.chi2fit;
%     lambdas(j,:) = lambda.upper~=0 | lambda.lower~=0;
%         
%     if(ar.chi2fit - chi2Reset > 2*dchi2)
%         fprintf('PLE #%i reached confidence limit\n', jk); 
%         break;
%     end
%     
%     if(lambda.upper(jk)~=0)
%         fprintf('PLE #%i reached upper parameter bound\n', jk);
%         break;
%     end
%     if(lambda.lower(jk)~=0)
%         fprintf('PLE #%i reached lower parameter bound\n', jk);
%         break;
%     end
%     
% %     if(norm(deltap) <= 1e-3)
% %         dp = dp * 1.1;
% %         addres(jk) = -dp*direction;
% %     end
% end
% 
% ar.p = pReset;
% 
