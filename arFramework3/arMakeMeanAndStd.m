% converts data from repeated measurements to mean and mean error data
%
% function arMakeMeanAndStd(m, d)

function arMakeMeanAndStd(m, d)

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('m','var'))
    for jm=1:length(ar.model)
        for jd=1:length(ar.model(jm).data)
            arMakeMeanAndStd(jm, jd)
        end
    end
    return
end
if(~exist('d','var'))
    for jd=1:length(ar.model(m).data)
        arMakeMeanAndStd(m, jd)
    end
    return
end

t = ar.model(m).data(d).tExp;
y = ar.model(m).data(d).yExp;

tnew = unique(t);
ynew = nan(length(tnew), size(y,2));
ystdnew = nan(length(tnew), size(y,2));

for j=1:length(tnew)
	q = t==tnew(j);	
	for jj=1:size(y,2)
		ytmp = y(q,jj);
		qnonnan = ~isnan(ytmp);
		if(sum(qnonnan)>1)
			ynew(j,jj) = mean(ytmp(qnonnan));
			ystdnew(j,jj) = std(ytmp(qnonnan)) / sqrt(sum(qnonnan)); % SEM standard error of mean not STD standard deviation!
        elseif(sum(qnonnan)==1)
            ynew(j,jj) = ytmp(qnonnan);
			ystdnew(j,jj) = ytmp(qnonnan);
            warning('arMakeMeanAndStd: only single replicates!');
		end
	end
end

ar.model(m).data(d).tExp = tnew;
ar.model(m).data(d).yExp = ynew;
ar.model(m).data(d).yExpStd = ystdnew;

if(isfield(ar.model(m), 'condition'))
    arLink(true);
end