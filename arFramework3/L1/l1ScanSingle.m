% L1 scan, 1 parameter after the other
% jks       relative parameters to be investigated by L1 regularization
% linv      width, i.e. inverse slope of L1 penalty (Inf = no penalty; small values = large penalty)
% gradient  use a small gradient on L1 penalty ([-1 0 1]; default = 0)

function l1ScanSingle(jks, linv, gradient)

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
    linv = [1e5 1e-5];
end
linv = sort(linv([1 end]),'descend');

if(~exist('gradient','var') || isempty(gradient))
    gradient = 0;
end

jks = sort(jks);

arWaitbar(0);

ps = nan(length(linv),length(ar.p));
chi2s = nan(1,length(linv));
chi2fits = nan(1,length(linv));
for i = 1:length(linv)
    [p,chi2,chi2fit] = l1Fit(jks,linv(i),gradient);
    ps(i,:) = p;
    chi2s(i) = chi2;
    chi2fits(i) = chi2fit;
end

parsgt0 = sum(abs(ps(:,jks)) > 1e-4,2);
uni_gt0 = unique(parsgt0);

maxiter = 100;
i = 1;
while length(uni_gt0) < length(jks)+1 && i < maxiter
    i = i + 1;
    arWaitbar(length(uni_gt0), length(jks)+1, sprintf('L1 scan one by one'));
    toRefine = find(abs(diff(parsgt0))>1);
    linvtmp = logspace(log10(linv(toRefine(1))),log10(linv(toRefine(1)+1)),3);
    linv = [linv(1:toRefine(1)) linvtmp(2) linv(toRefine(1)+1:end)];
    
    ar.std(jks) = linvtmp(2) * (1 + gradient * linspace(0,.001,length(jks)));
    ar.p = ps(toRefine(1),:);
    try
        arFit(true)
    catch exception
        fprintf('%s\n', exception.message);
    end
    ps= [ps(1:toRefine(1),:); ar.p; ps(toRefine(1)+1:end,:)];
    chi2s = [chi2s(1:toRefine(1)) ar.chi2+ar.chi2err-ar.chi2prior chi2s(toRefine(1)+1:end)];
    chi2fits = [chi2fits(1:toRefine(1)) ar.chi2./ar.config.fiterrors_correction+ar.chi2err chi2fits(toRefine(1)+1:end)];
    
    % Backward implementation
    for i = 1:length(linv);
        j = i;
        if j > 1
            while chi2fits(j) < max(chi2fits(1:j-1))-1e-3
                j = j-1;
                ar.std(jks) = linv(j) * (1 + gradient * linspace(0,.001,length(jks)));
                try
                    arFit(true)
                catch exception
                    fprintf('%s\n', exception.message);
                end
                ps(j,:) = ar.p;
                chi2s(j) = ar.chi2+ar.chi2err-ar.chi2prior;
                chi2fits(j) = ar.chi2./ar.config.fiterrors_correction+ar.chi2err;
                if j == 1
                    break
                end
            end
        end

        if sum(abs(ps(i,jks)) > 1e-6) == 0;
            ps(i+1:end,:) = repmat(ar.p,size(ps,1)-i,1);
            chi2s(i+1:end) = ar.chi2+ar.chi2err-ar.chi2prior;
            chi2fits(i+1:end) = ar.chi2./ar.config.fiterrors_correction+ar.chi2err;
            break
        end
    end
    
    parsgt0 = sum(abs(ps(:,jks)) > 1e-4,2);
    uni_gt0 = unique(parsgt0);
end

arWaitbar(-1);

[~, ind_1by1] = unique(parsgt0);
ind_1by1 = ind_1by1(end:-1:1);
ps = ps(ind_1by1,:);
chi2s = chi2s(ind_1by1);
chi2fits = chi2fits(ind_1by1);
linv = linv(ind_1by1);

ar.L1ps = ps;
ar.L1chi2s = chi2s;
ar.L1chi2fits = chi2fits;
ar.L1linv = linv;

function [p,chi2,chi2fit] = l1Fit(jks,linv,gradient)

global ar
ar.std(jks) = linv(1) * (1 + gradient * linspace(0,.001,length(jks)));
try
    arFit(true)
catch exception
    fprintf('%s\n', exception.message);
end
p = ar.p;
chi2 = ar.chi2+ar.chi2err-ar.chi2prior;
chi2fit = ar.chi2./ar.config.fiterrors_correction+ar.chi2err;

