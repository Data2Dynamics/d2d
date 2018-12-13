% L1 scan
% jks    relative parameters to be investigated by L1 regularization
% linv   width, i.e. inverse slope of L1 penalty (Inf = no penalty; small values = large penalty)

function l1Unpen(jks)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.type == 3);
    if(isempty(jks))
        error('please initialize by l1Init and run l1scan')
    end
end

if(~isfield(ar,'linv') || isempty(ar.linv))
    error('please initialize by l1Init and run l1scan')
end

linv = ar.linv;

ps = ar.L1ps;

ps_unpen = nan(length(linv),length(ar.p));
chi2s_unpen = nan(1,length(linv));

arWaitbar(0);

for i = 1:length(linv)
   arWaitbar(i, length(linv), sprintf('Unpenalized solution'));
   ar.p = ps(i,:);
   ar.type(jks) = 0;
   ar.qFit(jks) = 1;
   excl = jks(abs(ps(i,jks)) <= ar.L1thresh);
   ar.qFit(excl) = 2;
   try
       arFit(true)
       ps_unpen(i,:) = ar.p;
       chi2s_unpen(i) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
   catch exception
       fprintf('%s\n', exception.message);
   end

   j = i;
   if j == 1

       if (chi2s_unpen(j) < ar.L1lam0chi2s && isempty(excl))
           ar.L1lam0chi2s = chi2s_unpen(j);
       end
           
   else
       while chi2s_unpen(j) < max(chi2s_unpen(1:j-1))-1e-3
           j = j-1;
           ar.type(jks) = 0;
           ar.qFit(jks) = 1;
           ar.qFit(jks(abs(ps(j,jks)) <= ar.L1thresh)) = 2;

           try
               arFit(true)
               ps_unpen(j,:) = ar.p;
               chi2s_unpen(j) = arGetMerit('chi2')./ar.config.fiterrors_correction+arGetMerit('chi2err');
           catch exception
               fprintf('%s\n', exception.message);
           end

           if j == 1
               if sum(abs(ps(j,jks)) <= ar.L1thresh) == 0
                    ar.L1lam0chi2s = chi2s_unpen(j);
               end
               break
           end
       end
   end

end


arWaitbar(-1);

ar.L1jks = jks;
ar.L1ps_unpen = ps_unpen;
ar.L1chi2s_unpen = chi2s_unpen;
