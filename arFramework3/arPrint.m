% print parameter values

function arPrint

global ar

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
ar.qCloseToBound(~qLog10 & qPos) = log10(ar.p(~qLog10 & qPos)) - log10(ar.lb(~qLog10 & qPos)) < ar.config.par_close_to_bound | ...
    log10(ar.ub(~qLog10 & qPos)) - log10(ar.p(~qLog10 & qPos)) < ar.config.par_close_to_bound;

fprintf('Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model\n\n');
printHead;
for j=1:length(ar.p)
    printPar(j, ar.qCloseToBound(j));
	if(mod(j,10)==0 && j<length(ar.p))
		fprintf('     |   |                                |                                |              |         |      \n');
	end
end

    function printHead
        strhead = '     |   | name                           | lb       value       ub        | 10^value      | fitted | prior\n';
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
        fprintf(outstream, '%s%4i|%s%s%s| %30s | %+8.2g   %+8.2g   %+8.2g | %i   %+8.2g | %7i | %s \n', ...
            fit_flag, j, strdyn, strinit, strerr, ar.pLabel{j}, ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), pTrans(j), ar.qFit(j), priorStr(j));
        
    end

    function str = priorStr(j)
        if(ar.type(j) == 0)
            str = sprintf('uniform(%g,%g)', ar.lb(j), ar.ub(j));
        elseif(ar.type(j) == 1)
            str = sprintf('normal(%g,%g^2)', ar.mean(j), ar.std(j));
        elseif(ar.type(j) == 2)
            str = sprintf('uniform(%g,%g) with soft bounds', ar.lb(j), ar.ub(j));
        end
    end
end