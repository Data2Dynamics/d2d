% load model parameters and parameter setting from .mat
% and reconcile with current parameters
%
% arLoadPars(filename, fixAssigned)

function arLoadPars(filename, fixAssigned)

global ar

if(nargin<1)
    [~, filename] = fileChooser('./Results', 1, true);
	filename = ['./Results/' filename '/workspace.mat'];
else
    if(strcmp(filename,'end'))
        filelist = fileList('./Results');
        filename = filelist{end};
    end
    filename = ['./Results/' filename '/workspace.mat'];
end

if(~exist('fixAssigned', 'var'))
    fixAssigned = false;  
end

S = load(filename);
fprintf('parameters loaded from file %s:\n', filename);

for j=1:length(ar.p)
    qi = strmatch(ar.pLabel{j}, strvcat(S.ar.pLabel), 'exact'); %#ok<MATCH3,REMFF1>
    
    if(isempty(qi))
        fprintf('                      %s\n', ar.pLabel{j});
    else
        ar.p(j) = S.ar.p(qi);
        ar.qLog10(j) = S.ar.qLog10(qi);
        ar.qFit(j) = S.ar.qFit(qi);
        ar.lb(j) = S.ar.lb(qi);        
        ar.ub(j) = S.ar.ub(qi);        
        
        if(fixAssigned)
            ar.qFit(j) = 0;
            fprintf('fixed and assigned -> %s\n', ar.pLabel{j});
        else
            fprintf('          assigned -> %s\n', ar.pLabel{j});
        end
    end
end
