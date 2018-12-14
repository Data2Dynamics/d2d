% [X,Xall,annoX] = arGetDesignMatrix([xtype])
% 
% Sensitivity-Matrix of the observations as e.g. used in the definition of
% sloppiness extracted from ar.model.data.syExpSimu
%  
% Dimension: ndata x n_qFit
% Corresponds to the design matrix X in the linear case H == X'*X;
%
%   xtype   type of sensitivities, eg 'sres', 'sreslog', 'sy' ['sy']
% 
%   X       only data points with ar.model.data.qFit==1
%   Xall    all data 
%   annoX   assignment of the rows
%           annoX.m     model index
%           annoX.d     data index
%           annoX.it    time index
%           annoX.iy    observable index
%           annoX.dop   find(ar.qFit==1)
% 
% Attention: The sensitivities have to be calculated first by arSimu for
% the current parameter vector.

function [X,Xall,annoX] = arGetDesignMatrix(xtype)
if(~exist('xtype','var') || isempty(xtype))
    xtype = 'sy';
end

switch lower(xtype)
    case {'sres',''}
        [X,Xall,annoX] = arGetDesignMatrixSres(false);
    case {'sreslog'}
        [X,Xall,annoX] = arGetDesignMatrixSres(true);
    case {'sy','','syexpsimu'}
        [X,Xall,annoX] = arGetDesignMatrixSy;
    otherwise
        error('arGetDesignMatrix.m: xtype=%s unknown.',xtype);
end



function [X,Xall,annoX] = arGetDesignMatrixSres(forcelog)
if ~exist('forcelog','var') || isempty(forcelog)
    forcelog = false;
end

arCheckCache(true);
arCalcMerit(true);

global ar
nonlog = find(ar.qLog10~=1 & ar.qFit==1);

Xall = ar.sres;
% dres/dplog = dres/dp * dp/dplog
% dp/dplog = 10^plog * ln(10) = p * ln(10)
if forcelog && ~isempty(nonlog)
    punlog = ar.p;
    Xall = Xall.*(ones(size(Xall,1),1)*punlog) * log(10);
end


dop = find(ar.qFit==1);
X = zeros(size(ar.sres));
X(:,dop) = ar.sres(:,dop);

% dres/dplog = dres/dp * dp/dplog
% dp/dplog = 10^plog * ln(10) = p * ln(10)
if forcelog && ~isempty(nonlog)
    punlog = ar.p(nonlog);
    X(:,nonlog) = X(:,nonlog).*(ones(size(X,1),1)*punlog) * log(10);
end

annoX.m = [];
annoX.d = [];
annoX.it = [];
annoX.iy = [];
annoX.dop = dop;
annoX.stype = 'sy';





% Design matrix taken from ar.model.data.syExpSimu
function [X,Xall,annoX] = arGetDesignMatrixSy
global ar

dop = find(ar.qFit==1);

X = [];
Xall = [];
annoX.m = [];
annoX.d = [];
annoX.it = [];
annoX.iy = [];
for m=1:length(ar.model)
    for d = 1:length(ar.model(m).data)
        [~,ia,ib] = intersect(ar.model(m).data(d).p,ar.pLabel);
        dop_d = find(ar.qFit(ib)==1);
        if(sum(ar.model(m).data(d).qFit)>0)            
            s = ar.model(m).data(d).syExpSimu(:,ar.model(m).data(d).qFit,dop_d);
            tmp = reshape(s,size(s,2)*size(s,1),size(s,3));
            indnan = isnan(ar.model(m).data(d).yExp);
            tmp(indnan(:),:) = NaN;

            Xnew = zeros(size(tmp,1),length(ar.pLabel));
            
            Xnew(:,dop_d) = tmp;
            
            X = [X;Xnew];
            annoX.m = [annoX.m;m*ones(size(Xnew))];
            annoX.d = [annoX.d;d*ones(size(Xnew))];
            it = (1:size(s,1))'* ones(1,size(s,2));
            annoX.it = [annoX.it;it(:)];
            iy = ones(size(s,1),1)* (1:size(s,2));
            annoX.iy = [annoX.iy;iy(:)];
        end
        s = ar.model(m).data(d).syExpSimu(:,:,dop_d);
        s = reshape(s,size(s,2)*size(s,1),size(s,3));
        Xtmp = zeros(size(s,1),length(ar.pLabel));
        Xtmp(:,dop_d) = s;
        Xall = [Xall;Xtmp];
    end
end

drin = nansum(abs(X),2)>1e-10;
X = X(drin,:);
annoX.m = annoX.m(drin);
annoX.d = annoX.d(drin);
annoX.it = annoX.it(drin);
annoX.iy = annoX.iy(drin);
annoX.dop = dop;
annoX.stype = 'sy';

Xall = Xall(sum(abs(Xall),2)>1e-10,:);


