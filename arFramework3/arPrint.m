%   arPrint
%   arPrint(3:4)
%   arPrint('sd')
%   arPrint({'para1','para2'})
% 
% print parameter values
% 
%   js      Indices of the parameters to be displayed
%           (see Example below)
% 
% 
% Examples:
% arPrint('turn')
%            name                      lb       value       ub          10^value        fitted   prior
% #  22|DI | geneA_turn              |       -5      -0.41         +3 | 1      +0.39 |       1 | uniform(-5,3) 
% #  31|DI | geneB_turn              |       -5       -2.7         +3 | 1    +0.0022 |       1 | uniform(-5,3) 
% #  40|DI | geneC_turn              |       -5      -0.91         +3 | 1      +0.12 |       1 | uniform(-5,3) 
%      |   |                         |                                |              |         |      
% #  49|DI | geneD_turn              |       -5       -2.1         +3 | 1     +0.008 |       1 | uniform(-5,3) 
% #  58|DI | geneE_turn              |       -5       -2.5         +3 | 1     +0.003 |       1 | uniform(-5,3) 
% 
% arPrint(15:17)
% Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model
% 
% Example:
%            name                      lb       value       ub          10^value        fitted   prior
% #  15|D  | geneA_deg1              |       -2         -2         +3 | 1      +0.01 |       1 | uniform(-2,3) 
% #  16|D  | geneA_deg2              |       -2      +0.36         +3 | 1       +2.3 |       1 | uniform(-2,3) 
% #  17|D  | geneA_deg3              |       -2       -1.2         +3 | 1     +0.063 |       1 | uniform(-2,3) 

function varargout = arPrint(js)
global ar

if(~exist('js','var') | isempty(js))
    js = 1:length(ar.p);
elseif(islogical(js))
    js = find(js);
elseif(isnumeric(js))
    if(size(js,1)>1)
        js = js'; %should not be a row
    end
    if(sum(js==1 | js==0) ==length(js) && length(js)>1)
        js = find(js);
    end
elseif(ischar(js))
    js = find(~cellfun(@isempty,regexp(ar.pLabel,js)));
    if isempty(js)
        disp('Pattern not found in ar.pLabel');
        return;
    end
elseif(iscell(js)) % cell of pNames
    [~,js] = intersect(ar.pLabel,js);
else
    error('Argument has to be a string or an array of indices.')
end

if(sum(isnan(js))>0 || sum(isinf(js))>0 || min(js)<1 || max(js-round(js))>eps)
    js
    warning('arPrint.m: argument js is not plausible (should be an array of indices).')
else
    if(size(js,1)~=1)
        js = js';
    end
end

if nargout>0
    varargout{1} = js;
end

if(isempty(ar))
    error('please initialize by arInit')
end

pTrans = ar.p;
pTrans(ar.qLog10==1) = 10.^pTrans(ar.qLog10==1);

qLog10 = ar.qLog10 == 1;
ar.qCloseToBound(qLog10) = ar.p(qLog10) - ar.lb(qLog10) < ar.config.par_close_to_bound | ...
    ar.ub(qLog10) - ar.p(qLog10) < ar.config.par_close_to_bound;
qPos = ar.p>0;
ar.qCloseToBound(~qLog10 & ~qPos) = ar.p(~qLog10 & ~qPos) - ar.lb(~qLog10 & ~qPos) < ar.config.par_close_to_bound | ...
    ar.ub(~qLog10 & ~qPos) - ar.p(~qLog10 & ~qPos) < ar.config.par_close_to_bound;
ar.qCloseToBound(~qLog10 & qPos) = (ar.p(~qLog10 & qPos)) - (ar.lb(~qLog10 & qPos)) < ar.config.par_close_to_bound | ...
    (ar.ub(~qLog10 & qPos)) - (ar.p(~qLog10 & qPos)) < ar.config.par_close_to_bound;

maxlabellength = max(cellfun(@length, ar.pLabel));

fprintf('Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model\n\n');
printHead;
for j=js
    printPar(j, ar.qCloseToBound(j));
	if(mod(j,10)==0 && j<max(js))
		fprintf(['     |   | ' arExtendStr('', maxlabellength) ' |                                |              |         |      \n']);
	end
end

    function printHead
        strhead = ['     |   | ' arExtendStr('name', maxlabellength) ' | lb       value       ub        | 10^value      | fitted | prior\n'];
        strhead = strrep(strhead, '|', ' ');
        fprintf(strhead);
    end

    function printPar(j, qclosetobound)
        strdyn = ' ';
        if(ar.qDynamic(j))
            strdyn = 'D';
        end
        strerr = ' ';
        if(ar.qError(j))
            strerr = 'E';
        end
        strinit = ' ';
        if(ar.qInitial(j))
            strinit = 'I';
        end
        
        if(qclosetobound)
            outstream = 2;
        else
            outstream = 1;
        end
        if(ar.qFit(j)==2)
            fit_flag = 'C';
        else
            fit_flag = '#';
        end
        fprintf(outstream, '%s%4i|%s%s%s| %s | %+8.2g   %+8.2g   %+8.2g | %i   %+8.2g | %7i | %s \n', ...
            fit_flag, j, strdyn, strinit, strerr, arExtendStr(ar.pLabel{j}, maxlabellength), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), pTrans(j), ar.qFit(j), priorStr(j));
        
    end

    function str = priorStr(j)
        if(ar.type(j) == 0)
            str = sprintf('uniform(%g,%g)', ar.lb(j), ar.ub(j));
        elseif(ar.type(j) == 1)
            str = sprintf('normal(%g,%g^2)', ar.mean(j), ar.std(j));
        elseif(ar.type(j) == 2)
            str = sprintf('uniform(%g,%g) with soft bounds', ar.lb(j), ar.ub(j));
        elseif(ar.type(j) == 3)
            str = sprintf('L1(%g,%g)', ar.mean(j), ar.std(j));
        end
    end
end



