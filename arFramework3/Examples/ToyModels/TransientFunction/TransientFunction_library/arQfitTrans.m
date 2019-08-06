% arQfitTrans(option,varargin)
% 
% Examples:
% 
% arQfitTrans(1data',m,d)   % for one m/d, data.qFit=1, for all other data.qFit=0
% 
% arQfitTrans('1model',m);  % only data and parameters of model m are fitted
% 
function arQfitTrans(option,varargin)
global ar
    

switch option
    case '1data'  % for one m/d, data.qFit=1, for all other data.qFit=0
        mIn = varargin{1};
        dIn = varargin{2};

        psym = arSym(ar.pLabel);
        vars = symvar(arSym(ar.model(mIn).data(dIn).fp)); % using fp is a bit critical
        [~,found] = intersect(psym,vars);
        
        qFitIn = ar.qFit;
        ar.qFit(ar.qFit==1) = 0;
        ar.qFit(intersect(find(qFitIn==1 | qFitIn==0),found)) = 1; % only parameters in fy are fitted
        
        for m=1:length(ar.model)
            for d=1:length(ar.model(m).data)
                if m==mIn && d==dIn
                    ar.model(m).data(d).qFit(:)=1;                    
                else
                    ar.model(m).data(d).qFit(:)=0;
                end
            end
        end
    case '1model'  % for one m/d, data.qFit=1, for all other data.qFit=0
        mIn = varargin{1};

        
        for m=1:length(ar.model)
            for d=1:length(ar.model(m).data)
                if m==mIn 
                    if d==1
                        vars = symvar(arSym(ar.model(m).data(d).fp)); % using fp is a bit critical
                    else
                        vars = union(vars,symvar(arSym(ar.model(m).data(d).fp))); % using fp is a bit critical
                    end                        
                    ar.model(m).data(d).qFit(:)=1;                    
                else
                    ar.model(m).data(d).qFit(:)=0;
                end
            end
        end
        
        psym = arSym(ar.pLabel);
        [~,found] = intersect(psym,vars);
        
        qFitIn = ar.qFit;
        ar.qFit(ar.qFit==1) = 0;
        ar.qFit(intersect(find(qFitIn==1 | qFitIn==0),found)) = 1; % only parameters in fy are fitted

end