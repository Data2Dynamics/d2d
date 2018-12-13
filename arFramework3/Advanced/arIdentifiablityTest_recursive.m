% arIdentifiablityTest_recursive(varargin)
% 
% arIdentifiablityTest_recursive iteratively applies arIdentifiabilityTest while fixing
% the non-identifiability with largest dp (most flat direction). 
% 
%   varargin      these arguments are passed to arIdentifiablityTest.m
% 
% For each non-identifiability, the related parameters are detemined via
% L1-penalization. This step, however, currently only works for parameters without priors
% (i.e. ar.type=0)
% 
% See also arIdentifiablityTest
% 
% References: C. Kreutz, An easy and efficient approach for testing
% identifiability of parameters. Bioinformatics, 2018, 34(11), 1913-1921. 

function arIdentifiablityTest_recursive(varargin)

global ar
ar.NI = cell(0);
tic
starttime = cputime;
ar.NI_timing = [];
qfit0 = ar.qFit;

itest = 1;
mess = sprintf('IdentifiabilityTest iteration %i ...\n',itest);
fprintf(mess);
ar.NI{itest}.messages = {mess};

arIdentifiablityTest(varargin{:})
qfit = qfit0;

while ~isempty(ar.IdentifiabilityTest.p_nonIdentifiable)
    ar.NI{itest}.IdentifiabilityTest = ar.IdentifiabilityTest;
    ar.NI{itest}.pLabel = ar.IdentifiabilityTest.p_nonIdentifiable{1};
    ar.NI{itest}.p_val = ar.IdentifiabilityTest.p(ar.IdentifiabilityTest.dp_order(1));
    ar.NI{itest}.qFit = ar.p;
    ar.NI{itest}.qFit = ar.qFit;
    ar.NI{itest}.whichP = arIdentifiablityTest_WhichP(ar.NI{itest}.pLabel,ar.NI{itest}.p_val);

    %% next one ...
    itest = itest+1;    
    mess = sprintf('Fix the non-identifiable parameter %s ...\n',ar.IdentifiabilityTest.p_nonIdentifiable{1});
    fprintf(mess);
    ar.NI{itest}.messages = {mess};
    
    [~,ia] = intersect(ar.pLabel,ar.IdentifiabilityTest.p_nonIdentifiable{1});
    try
        qfit(ia) = 0;
        ar.qFit = qfit;
        
        mess = sprintf('IdentifiabilityTest iteration %i ...\n',itest);
        fprintf(mess);
        ar.NI{itest}.messages{end+1} = mess;
        arIdentifiablityTest(varargin{:})
    catch ERR
        ar.qFit = qfit0;
        rethrow(ERR);
    end
        
end
ar.qFit = qfit0;

if itest>1
    mess = sprintf('\narIdentifiablityTest_recursive: %i non-identifiabilities detected.\n\n',itest-1);
    fprintf(mess);
    ar.NI{itest-1}.messages{end+1} = mess;
    mess = sprintf('Fixing the following parameters yields an identifiable model:\n\n');
    fprintf(mess);
    ar.NI{itest-1}.messages{end+1} = mess;    
else
    ar.NI{1}.IdentifiabilityTest = ar.IdentifiabilityTest;
    mess = sprintf('\narIdentifiablityTest_recursive: No non-identifiabilities detected.\n\n');
    fprintf(mess);
    ar.NI{1}.messages{end+1} = mess;
end

ar.NI_timing = cputime - starttime;

