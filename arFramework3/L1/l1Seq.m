% L1 scan
% jks       relative parameters to be investigated by L1 regularization
% linv      width, i.e. inverse slope of L1 penalty (Inf = no penalty; small values = large penalty)
% gradient  use a small gradient on L1 penalty ([-1 0 1]; default = 0)

function l1Seq(jks, linv, gradient, lks,sorting)

global ar

%Check for trdog
checksum_l1   = {'67AF8E95E14AF615DFAB379CA542FD6C','f311e0c5dd8243c8e90166b03d48f17e','6BBE213BEC28A0C59A8DEDCF65CF7649'}; % Modified trdog.m
trpath = which('trdog','-all');
if sum(strcmpi(md5(trpath{1}),checksum_l1))==1
    % All good
else
    warning('Found an outdated trdog, updating...! \n');
    l1trdog
end

if(isempty(ar))
    error('please initialize by arInit')
end



if(~exist('jks','var') || isempty(jks))
    jks = find(ar.type == 3);
    if(isempty(jks))
        error('please initialize by l1Init')
    end
end

if(~exist('lks','var') || isempty(lks))
%     linv = logspace(-4,4,49);
%     linv = [linv Inf];
%     linv = linv(end:-1:1);
    lks = 1:length(ar.linv);
end

if (exist('linv','var') && ~isempty(linv))
    ar.linv = linv;
end

if(~exist('gradient','var') || isempty(gradient))
    gradient = 0;
end

if(~exist('sorting','var') || isempty(sorting))
    sorting = ones(size(lks));
else
    sorting = sorting * ones(size(lks));
end


jks = sort(jks);
optim = ar.config.optimizer;
maxiter = ar.config.optim.MaxIter;

if length(lks) > 1
    arWaitbar(0);
end

% arFit(true)

if (~isfield(ar,'L1ps') || isempty(ar.L1ps) )
    ps = nan(length(lks),length(ar.p));
else
    ps = ar.L1ps;
end

if (~isfield(ar,'L1chi2s') || isempty(ar.L1chi2s) )
    chi2s = nan(1,length(lks));
else
    chi2s = ar.L1chi2s;
end

if(~isfield(ar,'L1chi2fits') || isempty(ar.L1chi2fits) )
    chi2fits = nan(1,length(lks));
else
    chi2fits = ar.L1chi2fits;
end

chi2slam0 = ar.L1lam0chi2s;
% reference value for the discrepancy to the full model




chi2fits(:) = chi2slam0;




counter = 0;



for i = lks
    
    ar.std(jks) = ar.linv(i) * (1 + gradient * linspace(0,.001,length(jks)));
    switch ar.L1subtype(jks(1))
        case 1
            s = sprintf('L_1 scan');
        case 2
            ar.lnuweights(jks) = 1./(abs(ar.estim(jks)).^ar.gamma(i));
            s = sprintf('L_1/|OLS|^{%g} scan',ar.gamma(i));
        case 3
            ar.expo(jks) = ar.nu(i);
            s = sprintf('L_{%g} scan',ar.nu(i));
        case 4
            ar.alpha(jks) = ar.alpharange(i);
            s = sprintf('%g x L_1 + %g x L_2 scan',1-ar.alpharange(i),ar.alpharange(i));
    end
    
    if i > 1
        ar.p = ps(i-1,:);
    end
    
    arFit(true);
    
    ps(i,:) = ar.p;
    
    switch sorting(i) 
        case 1
            [~,p_sI] = sort(abs(ar.p(jks)-ar.mean(jks)));
        case 2
            p_sI = 1:length(jks);
    end
    % CASE 1 : the algorithm starts setting those values whose modulus
    % is closest to the mean value (e.g. 0)
    % CASE 2 : no sorting
    
    for j = jks(p_sI)
        
        counter = counter + 1;
        
        if ar.qFit(j) == 2
            % parameter already set to its mean value
            continue;
        end
        
        ar.type(jks) = 0;
        arWaitbar(counter, length(lks)*length(jks),...
            sprintf('Sequential %s',s));
        
                
        ar.type(j) = 3;
        
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
        
        if abs(ar.p(j)-ar.mean(j)) < ar.L1thresh
            ar.qFit(j) = 2;
        end
        
        
        ps(i,j) = ar.p(j);
  
        
    end
        
        
    
    

    
    if sum(ar.qFit(jks) == 1) == 0
        ps(i+1:end,:) = repmat(ar.p,size(ps,1)-i,1);
        chi2s(i+1:end) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
        chi2fits(i+1:end) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
        break
    end
end

arWaitbar(-1);



ar.L1ps = ps;
ar.L1chi2s_unpen = chi2s;
ar.L1chi2fits = chi2fits;

ar.config.optimizer = optim;
ar.config.optim.MaxIter = maxiter;

function md5hash = md5(filename)

mddigest   = java.security.MessageDigest.getInstance('MD5'); 
filestream = java.io.FileInputStream(java.io.File(filename)); 
digestream = java.security.DigestInputStream(filestream,mddigest);

while(digestream.read() ~= -1) end

md5hash=reshape(dec2hex(typecast(mddigest.digest(),'uint8'))',1,[]);