function ar = L1penalty(ar,L1sd,dolog,L1mean)

indP = [strmatch('k_',ar.pLabel);strmatch('khomo',ar.pLabel)];
ar.type(indP) = 3; % L1
ar.std(indP)  = L1sd;

% [stot,fak]=Abundance;
% fn = fieldnames(fak);

if(dolog==0)
    
    ar.lb(indP) = L1mean;
    ar.mean(indP) = L1mean; % L1
    
    for i=1:length(indP)
        if(ar.qLog10(indP(i))==1)
            ar.qLog10(indP(i)) = 0;
            ar.p(indP(i)) = 10.^ar.p(indP(i));
            ar.ub(indP(i)) = 10.^ar.ub(indP(i));
        end
    end
    
    
else  % dolog==1
        
    ar.lb(indP) = L1mean;
    ar.mean(indP) = L1mean; % L1

    for i=1:length(indP)
        if(ar.qLog10(indP(i))==0)
            ar.qLog10(indP(i)) = 1;
            ar.p(indP(i)) = log10(ar.p(indP(i)));
            ar.ub(indP(i)) = log10(ar.ub(indP(i)));
        end
    end
    
%     for f=1:length(fn)
%         ind = strmatch(fn{f},ar.pLabel,'exact');
%         ar.lb(ind) = -fak.(fn{f})+ar.lb(ind);
%         ar.ub(ind) = -fak.(fn{f})+ar.ub(ind);
%         ar.mean(ind) = -fak.(fn{f})+ar.mean(ind);
%     end
    
end

