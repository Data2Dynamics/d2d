% arHighlightY(jm,jd,jy,jt)
%   sets the highlight flag ar.model(jm).data(jd).highlight(jt,jy) = 1
%   This flag is used by arPlotY to highlight data points.
% 
% arHighlightY      sets all highlight flags to zero.
% arHighlightY(0)   sets all highlight flags to zero.

function arHighlightY(jm,jd,jy,jt)
global ar
if(~exist('jm','var') || isempty(jm))
    jm = [];
end
if(~exist('jd','var') || isempty(jd))
    jd = [];
end
if(~exist('jy','var') || isempty(jy))
    jy = [];
end
if(~exist('jt','var') || isempty(jt))
    jt = [];
end

if(nargin>0)
    if(jm==0)
        arHighlightY
        return
    end
    
    if(length(jm)==1 & isempty(jd))
        jd = 1:length(ar.model(jm).data);
    end
    if(isempty(jy) & length(jm)==1 & length(jd)==1)
        jy = size(ar.model(jm).data(jd).yExp,2);
    end
    if(isempty(jt) & length(jm)==1 & length(jd)==1)
        jt = size(ar.model(jm).data(jd).yExp,1);
    end
end

if(nargin==0)
    for m=1:length(ar.model)
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).highlight = zeros(size(ar.model(m).data(d).yExp));
        end
    end
else
    for m=jm
        for d=jd
            if(~isfield(ar.model(m).data(d),'highlight'))
                ar.model(m).data(d).highlight = zeros(size(ar.model(m).data(d).yExp));
            end
            if(sum(abs(size(ar.model(m).data(d).highlight)-size(ar.model(m).data(d).yExp)))~=0)
                ar.model(m).data(d).highlight = zeros(size(ar.model(m).data(d).yExp));
            end
            
            for y=jy
                for t=jt                        
                    if(t>size(ar.model(m).data(d).yExp,1))
                        [m,d,t,y]
                        ar.model(m).data(d).yExp
                        error('t>size(ar.model(m).data(d).yExp,1)')
                    end
                    if(y>size(ar.model(m).data(d).yExp,2))
                        [m,d,t,y]
                        ar.model(m).data(d).yExp
                        error('y>size(ar.model(m).data(d).yExp,2)')
                    end
                    ar.model(m).data(d).highlight(t,y) = 1;
                end
            end
        end
    end
end
