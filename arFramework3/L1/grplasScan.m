% L1 scan
% jks       relative parameters to be investigated by L1 regularization
% linv      width, i.e. inverse slope of L1 penalty (Inf = no penalty; small values = large penalty)
% gradient  use a small gradient on L1 penalty ([-1 0 1]; default = 0)

function grplasScan(jks, linv, gradient, lks, OptimizerSteps)

global ar

%Check for trdog
checksum_l1   = {'6BBE213BEC28A0C59A8DEDCF65CF7649'}; % Modified trdog.m
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

if(~isfield(ar,'grplas'))
    error('please initialize by grplasInit')
end


if(~exist('jks','var') || isempty(jks))
    jks = find(ar.type == 5);
    if(isempty(jks))
        error('please initialize by grplasInit')
    end
end

if(~exist('lks','var') || isempty(lks))
    lks = 1:length(ar.linv);
end

if (exist('linv','var') && ~isempty(linv))
    ar.linv = linv;
end

if(~exist('gradient','var') || isempty(gradient))
    gradient = 0;
end

if(~exist('OptimizerSteps','var') || isempty(OptimizerSteps))
    OptimizerSteps = [1000 20];
end



jks = sort(jks);
optim = ar.config.optimizer;
maxiter = ar.config.optim.MaxIter;

if length(lks) > 1
    arWaitbar(0);
end

% arFit(true)

if (~isfield(ar.grplas,'ps') || isempty(ar.grplas.ps) )
    ps = nan(length(lks),length(ar.p));
else
    ps = ar.grplas.ps;
end

if (~isfield(ar.grplas,'chi2s') || isempty(ar.grplas.chi2s) )
    chi2s = nan(1,length(lks));
else
    chi2s = ar.grplas.chi2s;
end

if(~isfield(ar.grplas,'chi2fits') || isempty(ar.grplas.chi2fits) )
    chi2fits = nan(1,length(lks));
else
    chi2fits = ar.grplas.chi2fits;
end

% ps(lks(1),:) = ar.p;
% chi2s(lks(1)) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
% chi2fits(lks(1)) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
for i = lks
    arWaitbar(find(lks == i ), length(lks), sprintf('Group Lasso Scan'));
    ar.std(jks) = ar.linv(i) * (1 + gradient * linspace(0,.001,length(jks)));

    
    if i > 1
        ar.p = ps(i-1,:);
    end

    try
        for o = 1:length(OptimizerSteps)
            if OptimizerSteps(o) > 0
                ar.config.optimizer = o;
                ar.config.optim.MaxIter = OptimizerSteps(o);
                arFit(true)
            end
        end
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
    
    if sum(abs(ps(i,jks)) > ar.grplas.thresh) == 0
        ps(i+1:end,:) = repmat(ar.p,size(ps,1)-i,1);
        chi2s(i+1:end) = arGetMerit('chi2')+arGetMerit('chi2err')-arGetMerit('chi2prior');
        chi2fits(i+1:end) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
        break
    end
end

arWaitbar(-1);

ar.grplas.ps = ps;
ar.grplas.chi2s = chi2s;
ar.grplas.chi2fits = chi2fits;

ar.config.optimizer = optim;
ar.config.optim.MaxIter = maxiter;

function md5hash = md5(filename)

mddigest   = java.security.MessageDigest.getInstance('MD5'); 
filestream = java.io.FileInputStream(java.io.File(filename)); 
digestream = java.security.DigestInputStream(filestream,mddigest);

while(digestream.read() ~= -1) end

md5hash=reshape(dec2hex(typecast(mddigest.digest(),'uint8'))',1,[]);