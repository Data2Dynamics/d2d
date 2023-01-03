function convertt = arMagnitudeConversion

% Magnitude_Conversion.m:
%   Switch order of magnitude of time in range 0-1000 
%   (if range condition allows, set dtmin > 1)

%   Switch order of magnitude of data in Range >10^-3 
%   or <-10^-3 (if negative sign)

global ar

tFine = ar.model(1).data.tFine;
y = ar.model(1).data.yFineSimu;

converty = ones(size(y,2),1);
for i=1:size(y,2)
        if (max(y(:,i)) > 10^(4) || max(y(:,i)) < 10^(-4))
            converty(i) = 10.^(floor(log10(range(y(:,i)))));
            y(:,i) = y(:,i)./converty(i); 
        end
end

% delete points of steady-state tail
% Estimate overall tRange and convertt between 50 and 600
convertt = ones(size(y,2),1);
% dt_end = ones(size(y,2),1).*max(tFine);
t = repmat(tFine,1,size(y,2));
for i=1:size(y,2)
    for j=size(y,1)-1:-1:2
        if abs(y(j+1,i)-y(j-1,i))/max(y(:,i)) < 10^(-6)
            y(j+1,i) = nan;
            t(j+1,i) = nan;
        else
            break
        end
    end
  % dt_end(i)=tFine(j);
    tRange = max(t(:,i));
    if tRange < 10
        convertt(i) = round(100/tRange,1,'significant');
    end
    if tRange > 100
        convertt(i) = round(100/tRange,1,'significant');   % so tRange is shorten to range(t) = 100;
    end
    t(:,i) = t(:,i) .* convertt(i);
end

ar.model(1).data.tFine = t;
ar.model(1).data.yFineSimu = y;

end

