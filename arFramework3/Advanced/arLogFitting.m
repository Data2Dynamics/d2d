% arLogFitting(logfitting)
%
% Set log-fitting for plotting.
%
%   logfitting      Enable logfitting?
function arLogFitting(logfitting)
global ar

for jm = 1:length(ar.model)
    if(isfield(ar.model(jm), 'data'))
        for jd = 1:length(ar.model(jm).data)
            if(~logfitting)
                ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==1) = ...
                    10.^ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==1);
            else
                ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==0) = ...
                    log10(ar.model(jm).data(jd).yExp(:,ar.model(jm).data(jd).logfitting==0));
            end
            ar.model(jm).data(jd).logfitting(:) = logfitting;
        end
    end
end