% Simulate for current parameter settings
%
% arSimu(sensi,fine)
%   sensi:          calculate sensitivities         [true]
%   fine:           fine grid for plotting          [false]

function arSimu(sensi,fine)

global ar

if(~exist('sensi', 'var'))
    sensi = true;
end
if(~exist('fine', 'var'))
    fine = false;
end

if(~isfield(ar,'p'))
    fprintf('WARNING: forgot linking\n');
    arLink;
end

if(~isfield(ar.config,'useParallel'))
    ar.config.useParallel = true;
end
if(~isfield(ar.config,'fiterrors_correction'))
    ar.config.fiterrors_correction = 1;
end

ar.stop = 0;

% propagate parameters
for m=1:length(ar.model)  
    for c=1:length(ar.model(m).condition)       
        ar.model(m).condition(c).status = 0;
        ar.model(m).condition(c).pNum = ar.p(ar.model(m).condition(c).pLink);
        ar.model(m).condition(c).qLog10 = ar.qLog10(ar.model(m).condition(c).pLink);
        ar.model(m).condition(c).pNum(ar.model(m).condition(c).qLog10 == 1) = ...
            10.^ar.model(m).condition(c).pNum(ar.model(m).condition(c).qLog10 == 1);
        ar.model(m).condition(c).start = 0;
        ar.model(m).condition(c).stop = 0;
        ar.model(m).condition(c).stop_data = 0;
    end
    
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).pNum = ar.p(ar.model(m).data(d).pLink);
            ar.model(m).data(d).qLog10 = ar.qLog10(ar.model(m).data(d).pLink);
            ar.model(m).data(d).pNum(ar.model(m).data(d).qLog10 == 1) = ...
                10.^ar.model(m).data(d).pNum(ar.model(m).data(d).qLog10 == 1);
        end
    end
end

% if(~exist(ar.fkt, 'file'))
%     fprintf('WARNING: file %s does not exist, recomplining...\n', ar.fkt);
%     arWriteCFiles;
% end

% initialize fine senistivities
if(fine && sensi)
    initFineSensis;
end

eval([ar.fkt '(ar, fine, ar.config.useSensis && sensi);'])

% integration error ?
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(ar.model(m).condition(c).status>0)
            error('arSimuCalc failed at %s', ar.info.arsimucalc_flags{ar.model(m).condition(c).status});
        elseif(ar.model(m).condition(c).status<0)
            error('cvodes failed at %s for model %i, condition %i', ...
                ar.info.cvodes_flags{abs(ar.model(m).condition(c).status)}, m, c);
        end
    end
end

% manually transform sensitivities from normal to log10 for fine time points
if(fine && sensi && ar.config.useSensis)
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).condition)
            for j=find(ar.qLog10(ar.model(m).condition(c).pLink)==1)
                ar.model(m).condition(c).suFineSimu(:,:,j) = ar.model(m).condition(c).suFineSimu(:,:,j) * ...
                    ar.model(m).condition(c).pNum(j) * log(10);
                ar.model(m).condition(c).sxFineSimu(:,:,j) = ar.model(m).condition(c).sxFineSimu(:,:,j) * ...
                    ar.model(m).condition(c).pNum(j) * log(10);
            end
        end
        if(isfield(ar.model(m), 'data'))
            for d=1:length(ar.model(m).data)
                for j=find(ar.qLog10(ar.model(m).data(d).pLink)==1)
                    ar.model(m).data(d).syFineSimu(:,:,j) = ar.model(m).data(d).syFineSimu(:,:,j) * ...
                        ar.model(m).data(d).pNum(j) * log(10);
                    ar.model(m).data(d).systdFineSimu(:,:,j) = ar.model(m).data(d).systdFineSimu(:,:,j) * ...
                        ar.model(m).data(d).pNum(j) * log(10);
                end
            end
        end
    end
end


% Initialize arrays for fine sensitivities with zeros
function initFineSensis

global ar

for m = 1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d = 1:length(ar.model(m).data)
            ar.model(m).data(d).syFineSimu = zeros(length(ar.model(m).data(d).tFine), length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
            ar.model(m).data(d).systdFineSimu = zeros(length(ar.model(m).data(d).tFine), length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
        end
    end
    for c = 1:length(ar.model(m).condition)
        ar.model(m).condition(c).suFineSimu = zeros(length(ar.model(m).condition(c).tFine), length(ar.model(m).u), length(ar.model(m).condition(c).p));
        ar.model(m).condition(c).svFineSimu = zeros(length(ar.model(m).condition(c).tFine), length(ar.model(m).vs), length(ar.model(m).condition(c).p));
        ar.model(m).condition(c).sxFineSimu = zeros(length(ar.model(m).condition(c).tFine), length(ar.model(m).x), length(ar.model(m).condition(c).p));
    end
end