% Link models, inputs and data sets
%
% arLink(silent, tExpAdd)


function arLink(silent, tExpAdd)

matVer = ver('MATLAB');

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(nargin==0)
    silent = false;
end

if(~silent)
    fprintf('\nlinking time points...\n');
end

for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        % clear condition time points
        for c=1:length(ar.model(m).condition)
            if(exist('tExpAdd','var'))
                ar.model(m).condition(c).tExp = tExpAdd;
            else
                ar.model(m).condition(c).tExp = [];
            end
            ar.model(m).condition(c).tFine = ...
                linspace(ar.model(m).tLim(1), ar.model(m).tLim(2), ar.config.nFinePoints)';
        end
        
        % collect time points
        for d=1:length(ar.model(m).data)
            if(isfield(ar.model(m).data(d), 'tExp'))
                ar.model(m).condition(ar.model(m).data(d).cLink).tExp = union(union( ... %R2013a compatible
                    ar.model(m).condition(ar.model(m).data(d).cLink).tExp, ...
                    ar.model(m).data(d).tExp), 0);
                ar.model(m).condition(ar.model(m).data(d).cLink).tExp = ...
                    ar.model(m).condition(ar.model(m).data(d).cLink).tExp(:);
            end
            
            ar.model(m).data(d).tFine = ...
                linspace(ar.model(m).data(d).tLim(1), ar.model(m).data(d).tLim(2), ar.config.nFinePoints)';
            
			ar.model(m).condition(ar.model(m).data(d).cLink).tFine = union( ... %R2013a compatible
			   ar.model(m).condition(ar.model(m).data(d).cLink).tFine, ...
			   ar.model(m).data(d).tFine);
            
            ar.model(m).condition(ar.model(m).data(d).cLink).tstart = ...
                min(ar.model(m).condition(ar.model(m).data(d).cLink).tFine);
        end
        
        % collect time points for multiple shooting
        if(isfield(ar, 'ms_count_snips') && ar.ms_count_snips>0)
            fprintf('\n');
            for jms=1:ar.model(m).ms_count
                for c=1:length(ar.model(m).condition)
                    for c2=1:length(ar.model(m).condition)
                        if(~isempty(ar.model(m).condition(c).ms_index) && ...
                                ~isempty(ar.model(m).condition(c2).ms_index))
                            qc = ar.model(m).condition(c).ms_index == jms;
                            qc2 = ar.model(m).condition(c2).ms_index == jms;
                            
                            if(sum(qc)>1 || sum(qc2)>1)
                                error('wrong multiple shooting indexing');
                            end
                            if(sum(qc)==1 && sum(qc2)==1 && ...
                                    ar.model(m).condition(c).ms_snip_index(qc)+1 == ar.model(m).condition(c2).ms_snip_index(qc2))
                                
                                tlink = ar.model(m).condition(c2).ms_snip_start;
                                fprintf('linking condition %i and %i for multiple shooting at t = %f\n', c, c2, tlink);
                                
                                ar.model(m).condition(c).tExp = ...
                                    union(ar.model(m).condition(c).tExp, tlink); %R2013a compatible
								ar.model(m).condition(c2).tExp = ...
                                    union(ar.model(m).condition(c2).tExp, tlink); %R2013a compatible 

                                
                                if(~isfield(ar.model(m), 'ms_link'))
                                    ar.model(m).ms_link = [c c2 tlink];
                                else
                                    ar.model(m).ms_link(end+1,1) = c;
                                    ar.model(m).ms_link(end,2) = c2;
                                    ar.model(m).ms_link(end,3) = tlink;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % link back time points
        for d=1:length(ar.model(m).data)
            if(isfield(ar.model(m).data(d), 'tExp'))
                if(str2double(matVer.Version)>=8.1)
                    [qtime, itime] = ismember(ar.model(m).data(d).tExp, ...
                        ar.model(m).condition(ar.model(m).data(d).cLink).tExp,'legacy'); %#ok<ASGLU>
                else
                    [qtime, itime] = ismember(ar.model(m).data(d).tExp, ...
                        ar.model(m).condition(ar.model(m).data(d).cLink).tExp); %#ok<ASGLU>
                end
                ar.model(m).data(d).tLinkExp = itime;
            end
            if(str2double(matVer.Version)>=8.1)
                [qtime, itime] = ismember(ar.model(m).data(d).tFine, ...
                    ar.model(m).condition(ar.model(m).data(d).cLink).tFine,'legacy'); %#ok<ASGLU>
            else
                [qtime, itime] = ismember(ar.model(m).data(d).tFine, ...
                    ar.model(m).condition(ar.model(m).data(d).cLink).tFine); %#ok<ASGLU>
            end
            ar.model(m).data(d).tLinkFine = itime;
        end
        
        % link back for multiple shooting
        if(isfield(ar.model(m), 'ms_link') && ~isempty(ar.model(m).ms_link))
            for jms = 1:size(ar.model(m).ms_link, 1)
                ar.model(m).ms_link(jms,4) = ...
                    find(ar.model(m).condition(ar.model(m).ms_link(jms,1)).tExp == ar.model(m).ms_link(jms,3));
                ar.model(m).ms_link(jms,5) = ...
                    find(ar.model(m).condition(ar.model(m).ms_link(jms,2)).tExp == ar.model(m).ms_link(jms,3));
            end
        end
        
        % statistics
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).qFit = true(size(ar.model(m).data(d).y));
            if(isfield(ar.model(m).data(d), 'tExp') && isfield(ar.model(m).data(d), 'yExp'))
                ar.model(m).data(d).ndata = sum(~isnan(ar.model(m).data(d).yExp),1);
            else
                ar.model(m).data(d).ndata = zeros(size(ar.model(m).data(d).y));
            end
            ar.model(m).data(d).chi2 = zeros(size(ar.model(m).data(d).ndata));
            ar.model(m).data(d).chi2err = zeros(size(ar.model(m).data(d).ndata));
        end
    else
        for c = 1:length(ar.model(m).condition)
            ar.model(m).condition(c).tFine = linspace(ar.model(m).tLim(1), ar.model(m).tLim(2), ar.config.nFinePoints);
            ar.model(m).condition(c).tstart = min(ar.model(m).condition(c).tFine);
        end
    end
end

% loading array for x, z and y
for m = 1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d = 1:length(ar.model(m).data)
            ntf = length(ar.model(m).data(d).tFine);
            ny = length(ar.model(m).data(d).y);
            np = length(ar.model(m).data(d).p);
            
            if(isfield(ar.model(m).data(d), 'tExp'))
                nt = length(ar.model(m).data(d).tExp);
                if(nt>0)
                    ar.model(m).data(d).yExpSimu = zeros(nt, ny);
                    ar.model(m).data(d).syExpSimu = zeros(nt, ny, np);
                    ar.model(m).data(d).ystdExpSimu = zeros(nt, ny);
                    ar.model(m).data(d).systdExpSimu = zeros(nt, ny, np);
                    if(isfield(ar.model(m).data(d), 'yExp') && ~isempty(ar.model(m).data(d).yExp))
                        ar.model(m).data(d).res = zeros(nt, ny);
                        ar.model(m).data(d).reserr = zeros(nt, ny);
                        ar.model(m).data(d).sres = zeros(nt, ny, np);
                        ar.model(m).data(d).sreserr = zeros(nt, ny, np);
                        ar.model(m).data(d).has_yExp = true;
                    else
                        ar.model(m).data(d).has_yExp = false;
                    end
                    ar.model(m).data(d).has_tExp = true;
                end
            else
                ar.model(m).data(d).has_tExp = false;
                ar.model(m).data(d).has_yExp = false;
            end
            
            ar.model(m).data(d).yFineSimu = zeros(ntf, ny);
            ar.model(m).data(d).ystdFineSimu = zeros(ntf, ny);
        end
    end
    for c = 1:length(ar.model(m).condition)
        ntf = length(ar.model(m).condition(c).tFine);
        nu = length(ar.model(m).u);
        nv = length(ar.model(m).vs);
        nx = length(ar.model(m).x);
        nz = length(ar.model(m).z);
        np = length(ar.model(m).condition(c).p);
        
        ar.model(m).condition(c).uNum = zeros(1, nu);
        ar.model(m).condition(c).vNum = zeros(1, nv);
        ar.model(m).condition(c).dvdxNum = zeros(nv, nx); 
        ar.model(m).condition(c).dvduNum = zeros(nv, nu); 
        ar.model(m).condition(c).dvdpNum = zeros(nv, np); 
        ar.model(m).condition(c).suNum = zeros(nu, np);
        ar.model(m).condition(c).svNum = zeros(1, nv);

        if(isfield(ar.model(m).condition(c), 'tExp'))
            nt = length(ar.model(m).condition(c).tExp);
            
            ar.model(m).condition(c).uExpSimu = zeros(nt, nu);
            ar.model(m).condition(c).suExpSimu = zeros(nt, nu, np);
            ar.model(m).condition(c).vExpSimu = zeros(nt, nv);
            ar.model(m).condition(c).svExpSimu = zeros(nt, nv, np);
            ar.model(m).condition(c).xExpSimu = zeros(nt, nx);
            ar.model(m).condition(c).sxExpSimu = zeros(nt, nx, np);
            ar.model(m).condition(c).zExpSimu = zeros(nt, nz);
            ar.model(m).condition(c).szExpSimu = zeros(nt, nz, np);
            ar.model(m).condition(c).has_tExp = true;
        else
            ar.model(m).condition(c).has_tExp = false;
        end
        
        ar.model(m).condition(c).uFineSimu = zeros(ntf, nu);
        ar.model(m).condition(c).vFineSimu = zeros(ntf, nv);
        ar.model(m).condition(c).xFineSimu = zeros(ntf, nx);
        ar.model(m).condition(c).zFineSimu = zeros(ntf, nz);
        
        % steady state sensitivities
        ar.model(m).condition(c).qSteadyState = false(1,nx);
        ar.model(m).condition(c).dxdt = zeros(1, nx);
        ar.model(m).condition(c).ddxdtdp = zeros(nx, np);
        ar.model(m).condition(c).stdSteadyState = zeros(1,nx) + ar.config.steady_state_constraint;
        
        ar.model(m).condition(c).start = 0;
        ar.model(m).condition(c).stop = 0;
        ar.model(m).condition(c).stop_data = 0;
    end
end

if(~silent)
    fprintf('linking parameters...\n');
end

% remember existing values
if(isfield(ar, 'pLabel'))
    plabel = ar.pLabel;
    p = ar.p;
    qfit = ar.qFit;
    qlog10 = ar.qLog10;
    lb = ar.lb;
    ub = ar.ub;
    type = ar.type;
    meanp = ar.mean;
    stdp = ar.std;
end

% collecting parameters
ar.pLabel = {};
for m = 1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d = 1:length(ar.model(m).data)
            if(str2double(matVer.Version)>=8.1)
                ar.pLabel = union(ar.pLabel, ar.model(m).data(d).p,'legacy'); %R2013a compatible
            else
                ar.pLabel = union(ar.pLabel, ar.model(m).data(d).p);
            end
        end
    end
    for c = 1:length(ar.model(m).condition)
        if(str2double(matVer.Version)>=8.1)
            ar.pLabel = union(ar.pLabel, ar.model(m).condition(c).p,'legacy'); %R2013a compatible                  
        else
            ar.pLabel = union(ar.pLabel, ar.model(m).condition(c).p);
        end
    end
end
ar.qFit = ones(size(ar.pLabel));

% determine parameters influencing model dynamics
ar.qDynamic = zeros(size(ar.pLabel));
for m = 1:length(ar.model)
    for c = 1:length(ar.model(m).condition)
        if(~isempty(ar.model(m).condition(c).p))
            qdyn = ismember(ar.pLabel, ar.model(m).condition(c).p); %R2013a compatible
            ar.qDynamic(qdyn) = 1;
        end
    end
end

% determine parameters influencing initial values of model dynamics
ar.qInitial = zeros(size(ar.pLabel));
for m = 1:length(ar.model)
    for c = 1:length(ar.model(m).condition)
        if(~isempty(ar.model(m).condition(c).p) && ~isempty(ar.model(m).condition(c).px0))
            qinit = ismember(ar.pLabel, ar.model(m).condition(c).px0); %R2013a compatible
            ar.qInitial(qinit) = 1;
        end
    end
end

% determine parameters influencing the error model
ar.qError = zeros(size(ar.pLabel));
for m = 1:length(ar.model)
    if(isfield(ar.model(m),'data'))
        for d = 1:length(ar.model(m).data)
            if(~isempty(ar.model(m).data(d).pystd))
                qerr = ismember(ar.pLabel, ar.model(m).data(d).pystd); %R2013a compatible
                ar.qError(qerr) = 1;
            end
        end
    end
end

% fix volumen parameters
for m = 1:length(ar.model)
    qvolpara = ismember(ar.pLabel, ar.model(m).pc); %R2013a compatible
    ar.qFit(qvolpara) = 2;
end

ar.qLog10 = ones(size(ar.pLabel));
ar.p = ones(size(ar.pLabel)) * -1;

ar.ub = ones(size(ar.pLabel)) * +3;
ar.lb = ones(size(ar.pLabel)) * -5;

ar.type = zeros(size(ar.pLabel));
ar.mean = zeros(size(ar.pLabel));
ar.std = zeros(size(ar.pLabel));

% link back parameters
for m = 1:length(ar.model)
    for c = 1:length(ar.model(m).condition)
        if(~isempty(ar.model(m).condition(c).p))
            ar.model(m).condition(c).pLink = ismember(ar.pLabel, ar.model(m).condition(c).p); %R2013a compatible
        else
            ar.model(m).condition(c).pLink = [];
        end
        ar.model(m).condition(c).pNum = 10.^ar.p(ar.model(m).condition(c).pLink);
    end
    if(isfield(ar.model(m), 'data'))
        for d = 1:length(ar.model(m).data)
            if(~isempty(ar.model(m).data(d).p))
                ar.model(m).data(d).pLink = ismember(ar.pLabel, ar.model(m).data(d).p); %R2013a compatible
            else
                ar.model(m).data(d).pLink = [];
            end
            ar.model(m).data(d).pNum = 10.^ar.p(ar.model(m).data(d).pLink);
        end
    end
end

% populate threads
ar.config.threads = [];
ar.config.threads(1).id = 0;
ar.config.threads(1).n = 0;
ar.config.threads(1).nd = 0;
ar.config.threads(1).ms = int32([]);
ar.config.threads(1).cs = int32([]);
ithread = 1;
ar.config.nTasks = 0;
for m = 1:length(ar.model)
    for c = 1:length(ar.model(m).condition)
        if(length(ar.config.threads)<ithread)
            ar.config.threads(ithread).id = ithread-1;
            ar.config.threads(ithread).n = 0;
            ar.config.threads(ithread).nd = 0;
            ar.config.threads(ithread).ms = int32([]);
            ar.config.threads(ithread).cs = int32([]);
        end
        ar.config.nTasks = ar.config.nTasks + 1;
        ar.config.threads(ithread).n = ar.config.threads(ithread).n + 1;
        ar.config.threads(ithread).nd = ar.config.threads(ithread).nd + ...
            length(ar.model(m).condition(c).dLink);        
        ar.config.threads(ithread).ms(end+1) = int32(m-1);
        ar.config.threads(ithread).cs(end+1) = int32(c-1);
        ithread = ithread + 1;
        if(ithread>ar.config.nParallel)
            ithread = 1;
        end
    end
end
ar.config.nThreads = length(ar.config.threads);

% reset values
if(exist('plabel','var'))
    arSetPars(plabel, p, qfit, qlog10, lb, ub, type, meanp, stdp);
end

% plotting
for jm=1:length(ar.model)
    if(~isfield(ar.model(jm), 'qPlotYs') || isempty(ar.model(jm).qPlotYs))
        if(length(ar.model(jm).plot) > 10)
            if(~silent)
                fprintf('Automatic plotting disabled for model %i. Please use arPlotter for plotting.\n', jm);
            end
            ar.model(jm).qPlotYs = false(1,length(ar.model(jm).plot));
            ar.model(jm).qPlotXs = false(1,length(ar.model(jm).plot));
            ar.model(jm).qPlotVs = false(1,length(ar.model(jm).plot));
        else
            ar.model(jm).qPlotYs = true(1,length(ar.model(jm).plot));
            if(isfield(ar.model(jm),'data'))
                ar.model(jm).qPlotXs = false(1,length(ar.model(jm).plot));
            else
                ar.model(jm).qPlotXs = true(1,length(ar.model(jm).plot));
            end
            ar.model(jm).qPlotVs = false(1,length(ar.model(jm).plot));
        end
    end
end

% set external parameters
if(isfield(ar, 'pExternLabels'))
    arSetPars(ar.pExternLabels, ar.pExtern, ar.qFitExtern, ar.qLog10Extern, ...
        ar.lbExtern, ar.ubExtern);
end

