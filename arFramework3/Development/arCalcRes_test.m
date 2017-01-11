%  This function can be used to check the new function arCalcRes.m
%  The following fields are recalulated with arCalcRes and compared to the
%  valued computed in arSimuCalc.c (in the deprecated version):
%  data.res, reserr, sres, sreserr, chi2, chi2err
% 
%  It works only whith the old version of arSimuCalc.c and if arCalcRes is
%  NOT called in arSimu (or elsewhere) befor this fuction is called.
% 
%  If it should be used later, then
%   1) exchange arSimuCalc.c with Deprecated/arSimuCalc_2016.c
%   2) comment the call of arCalcRes in arSimu.m
%   3) Call arChi2 (or arSimu with proper arguments)
%   4) Call arCalcRes_test

function [maxdiff,maxd] = arCalcRes_test

global ar
arOld = arDeepCopy(ar);
old = copy_fields(arOld);
fn = fieldnames(old);

arNew = cell(0);

try
    arCalcRes(0);
    arNew{1} = arDeepCopy(ar);
    arCalcRes(1);
    arNew{end+1} = arDeepCopy(ar);
    
    for i=1:length(arNew)
        new = copy_fields(arNew{i});
        
        maxdiff = 0;
        maxd = zeros(length(fn),length(ar.model),length(ar.model(1).data));
        for f=1:length(fn)
            for m=1:length(ar.model)
                for d=1:length(ar.model(m).data)
                    dd = old.(fn{f}){m}{d} - new.(fn{f}){m}{d};
                    maxd(f,m,d) = max(abs(dd(:)));
                end
            end
            tmp = squeeze(maxd(f,:,:));
            maxdiff = max(tmp(:),maxdiff);
        
            if max(tmp(:))>1e-10
                save error
                fprintf('%e  ',max(maxd));
                error('maxd >1e-10 for %s',fn{f})
            end
        end
        
%         fprintf('%e  ',max(maxd));
    end

catch
    ar = arOld;    
    rethrow(lasterror)
end

ar = arOld;

end

function old = copy_fields(ar)
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)

        old.res{m}{d} = ar.model(m).data(d).res + 0.0;
        old.reserr{m}{d} = ar.model(m).data(d).reserr + 0.0;
        old.chi2{m}{d} = ar.model(m).data(d).chi2 + 0.0;
        old.chi2err{m}{d} = ar.model(m).data(d).chi2err + 0.0;
        old.sres{m}{d} = ar.model(m).data(d).sres + 0.0;
        old.sreserr{m}{d} = ar.model(m).data(d).sreserr + 0.0;
    end
end


end