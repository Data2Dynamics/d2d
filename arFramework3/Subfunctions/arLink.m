% Link models, inputs and data sets
%
% arLink(silent, tExpAdd)


function arLink(silent, tExpAdd, dataAdd, ix, id, im, newData, yStd, add_sec)

matVer = ver('MATLAB');

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(nargin==0)
    silent = false;
end

if(~silent)
    arFprintf(1, '\nlinking model...');
    arFprintf(2, '\nlinking time points...\n');
end

if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

useMS 		= 0;
useEvents 	= 0;

%Set variables to temporarily add data points
if(~exist('add_sec','var'))
    add_sec = false;
end

if(~exist('dataAdd','var'))
    dataAdd = false;
    add_sec = false;
    newData = [];
    ix = [];
    id = [];
    im = [];
    yStd = [];
end

if(~exist('yStd','var'))
    if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==1) || ...
            (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(im,id)==1) )
        yStd = NaN;
    else
        yStd = 0.1;

    end
end

for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(~isfield(ar.model(m).condition(c), 'tEvents'))
            ar.model(m).condition(c).tEvents = [];
        end
    end
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
            % Initialize data fields needed in arCalcRes if they don't
            % exist
            if(~isfield(ar.model(m).data(d), 'yExp'))
                ar.model(m).data(d).yExp = [];
            end
            if(~isfield(ar.model(m).data(d), 'yExpSimu'))
                ar.model(m).data(d).yExpSimu = [];
            end             
            if(~isfield(ar.model(m).data(d), 'yExpStd'))
                ar.model(m).data(d).yExpStd = [];
            end
            if(~isfield(ar.model(m).data(d), 'systdExpSimu'))
                ar.model(m).data(d).systdExpSimu = [];
            end
            if(~isfield(ar.model(m).data(d), 'ystdExpSimu'))
                ar.model(m).data(d).ystdExpSimu = [];
            end
            
            if(isfield(ar.model(m).data(d), 'tExp'))
                %delete data point
                if(dataAdd && d==id && m==im && isfield(ar.model(m).data(d),'ppl'))
                    if((isfield(ar.model(m).data(d).ppl,'Added') && ~isempty(ar.model(m).data(d).ppl.Added)) ...
                        || (isfield(ar.model(m).data(d).ppl,'Added2') && ~isempty(ar.model(m).data(d).ppl.Added2)))
                        if(~add_sec && ~isempty(ar.model(m).data(d).ppl.Added))
                            ar.model(m).data(d).tExp(ar.model(m).data(d).ppl.Added)=[];
                            ar.model(m).data(d).yExp(ar.model(m).data(d).ppl.Added,:)=[];
                            ar.model(m).data(d).yExpStd(ar.model(m).data(d).ppl.Added,:)=[];
                            ar.model(m).data(d).ppl.Added=[];                       
                        end
                    end
                end
                ar.model(m).condition(ar.model(m).data(d).cLink).tExp = union(union( ... %R2013a compatible
                    ar.model(m).condition(ar.model(m).data(d).cLink).tExp, ...
                    ar.model(m).data(d).tExp), 0);
                ar.model(m).condition(ar.model(m).data(d).cLink).tExp = ...
                    ar.model(m).condition(ar.model(m).data(d) ...
                                          .cLink).tExp(:);
                %add new data point in id, im
                if(dataAdd && d==id && m==im && ~isnan(newData))
                    if(add_sec && isempty(ar.model(m).data(d).ppl.Added))
                        error('\n You have to add a first data point before adding a second \n')                       
                    end
                    infl_tExp = [ ar.model(m).data(d).tExp; tExpAdd ];
                    %[ar.model(m).data(d).tExp, iAdd, ic]=unique( infl_tExp );
                    ar.model(m).data(d).tExp=sort( infl_tExp );
                    new_data=NaN(1,size(ar.model(m).data(d).yExp,2));
                    new_data(ix) = newData;
                    [~,it_tmp] = min(abs(ar.model(m).data(d).tExp-tExpAdd));
                    if(length(find(ar.model(m).data(d).tExp==tExpAdd))==1)
                        ar.model(m).data(d).yExp =  [ar.model(m).data(d).yExp(1:(it_tmp-1),:); new_data; ar.model(m).data(d).yExp(it_tmp:end,:)];
                        new_data(ix) = yStd;
                        ar.model(m).data(d).yExpStd =  [ar.model(m).data(d).yExpStd(1:(it_tmp-1),:); new_data; ar.model(m).data(d).yExpStd(it_tmp:end,:)];
                    else
                        ar.model(m).data(d).yExp =  [ar.model(m).data(d).yExp(1:it_tmp,:); new_data; ar.model(m).data(d).yExp(it_tmp+1:end,:)];
                        new_data(ix) = yStd;
                        ar.model(m).data(d).yExpStd =  [ar.model(m).data(d).yExpStd(1:it_tmp,:); new_data; ar.model(m).data(d).yExpStd(it_tmp+1:end,:)];
                    end
                    %possibility to add a second data point
                    if(add_sec)
                        ar.model(m).data(d).ppl.Added2=find(ar.model(m).data(d).tExp==tExpAdd,1,'last');
                        if(~isempty(ar.model(m).data(d).ppl.Added) && ar.model(m).data(d).ppl.Added2<ar.model(m).data(d).ppl.Added)
                            ar.model(m).data(d).ppl.Added = ar.model(m).data(d).ppl.Added + 1;
                        end
                    else
                        ar.model(m).data(d).ppl.Added=find(ar.model(m).data(d).tExp==tExpAdd,1,'last');               
                    end
                    
                    %set variable to retrieve data point in ar.res
                    resi = 0;
                    for m_res=1:im-1
                        for d_res=1:length(ar.model(m_res).data)
                            resi = resi + sum(sum(~isnan(ar.model(m_res).data(d_res).yExp)));
                        end
                    end
                    for d_res = 1:id-1
                        resi = resi + sum(sum(~isnan(ar.model(im).data(d_res).yExp)));
                    end
                    for x_res = 1:ix-1
                        resi = resi + sum(~isnan(ar.model(im).data(id).yExp(:,x_res)));
                    end
                    ar.ppl.resi_tmp = resi + sum(~isnan(ar.model(im).data(id).yExp(1:ar.model(m).data(d).ppl.Added,ix)));
                    
                    if(add_sec)
                        ar.ppl.resi2_tmp = resi + sum(~isnan(ar.model(im).data(id).yExp(1:ar.model(m).data(d).ppl.Added2,ix)));
                    end
                end
            end
            
            ar.model(m).data(d).tFine = ...
                linspace(ar.model(m).data(d).tLim(1), ar.model(m).data(d).tLim(2), ar.config.nFinePoints)';
            
            % Add extra time points if desired
            if isfield(ar.model(m).data(d), 'tExtra')
                ar.model(m).data(d).tFine = union(ar.model(m).data(d).tFine, ar.model(m).data(d).tExtra);
            end
            
            ar.model(m).condition(ar.model(m).data(d).cLink).tFine = union( ... %R2013a compatible
                ar.model(m).condition(ar.model(m).data(d).cLink).tFine, ...
                ar.model(m).data(d).tFine);
            
            ar.model(m).condition(ar.model(m).data(d).cLink).tstart = ...
                min(ar.model(m).condition(ar.model(m).data(d).cLink).tFine);
            
            % Events defined in the data
            if(isfield(ar.model(m).data(ar.model(m).data(d).cLink), 'tEvents'))
                ar.model(m).condition(ar.model(m).data(d).cLink).tEvents = ...
                    union( ar.model(m).condition(ar.model(m).data(d).cLink).tEvents, ...
                    ar.model(m).data(d).tEvents );
            end
        end
        
        % collect time points for multiple shooting
        if(isfield(ar, 'ms_count_snips') && ar.ms_count_snips>0)
            arFprintf(2, '\n');
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
                                arFprintf(2, 'linking condition %i and %i for multiple shooting at t = %f\n', c, c2, tlink);
                                
                                % Add multiple shooting points to the event list
                                ar.model(m).condition(c).tEvents = ...
                                    union(ar.model(m).condition(c).tEvents, tlink); %R2013a compatible
                                ar.model(m).condition(c2).tEvents = ...
                                    union(ar.model(m).condition(c2).tEvents, tlink); %R2013a compatible
                                
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
        
        % Remove events before tstart
        for c = 1 : length( ar.model(m).condition )
            if ~isempty(ar.model(m).condition(c).tEvents)
                ar.model(m).condition(c).tEvents( ar.model(m).condition(c).tEvents < ar.model(m).condition(c).tstart ) = [];
            end
        end
        
        % Add events to tFine and tExp (if it exists)
        for c = 1 : length( ar.model(m).condition )
            if isfield(ar.model(m).condition(c), 'tExp')
                ar.model(m).condition(c).tExp = ...
                    union(ar.model(m).condition(c).tExp, ar.model(m).condition(c).tEvents);
            end
            
            ar.model(m).condition(c).tFine = ...
                union(ar.model(m).condition(c).tFine, ar.model(m).condition(c).tEvents);
        end
        
        % Add extra time points if desired
        if isfield(ar.model(m).condition(c), 'tExtra')
            ar.model(m).conditions(c).tFine = union(ar.model(m).condition(c).tFine, ar.model(m).condition(c).tExtra);
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
            % Only overwrite qFit if it doesn't exist yet or if it is empty
            if ((~isfield( ar.model(m).data(d), 'qFit' ))||(isempty(ar.model(m).data(d).qFit)))
                ar.model(m).data(d).qFit = true(size(ar.model(m).data(d).y));
            end
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
        
        % Add events to tFine and tExp (if it exists)
        for c = 1 : length( ar.model(m).condition )
            if isfield(ar.model(m).condition(c), 'tExp')
                ar.model(m).condition(c).tExp = ...
                    union(ar.model(m).condition(c).tExp, ar.model(m).condition(c).tEvents);
            end
            
            ar.model(m).condition(c).tFine = ...
                union(ar.model(m).condition(c).tFine, ar.model(m).condition(c).tEvents);
        end
        
        % Add extra time points if desired
        if isfield(ar.model(m).condition(c), 'tExtra')
            ar.model(m).conditions(c).tFine = union(ar.model(m).condition(c).tFine, ar.model(m).condition(c).tExtra);
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
                    if(length(ar.model(m).x)>0) %#ok
                        ar.model(m).data(d).y_scale = zeros(nt, ny, length(ar.model(m).x));                       
                    else
                        ar.model(m).data(d).y_scale = zeros(nt, ny, 1);                        
                    end
                    ar.model(m).data(d).dydt = zeros(nt, ny);
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
        ar.model(m).condition(c).y_atol = zeros(nx,1);
        
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
            ar.model(m).condition(c).dzdx = zeros(nt, nz, nx);
            ar.model(m).condition(c).dxdts = zeros(nt, nx);
            ar.model(m).condition(c).has_tExp = true;
        else
            ar.model(m).condition(c).has_tExp = false;
        end
        
        ar.model(m).condition(c).uFineSimu = zeros(ntf, nu);
        ar.model(m).condition(c).vFineSimu = zeros(ntf, nv);
        ar.model(m).condition(c).xFineSimu = zeros(ntf, nx);
        ar.model(m).condition(c).zFineSimu = zeros(ntf, nz);
        
        % event assignment/override operations
        nte = length( ar.model(m).condition(c).tEvents );
        
        modx_A = ones(nte, nx);
        modx_B = zeros(nte, nx);
        modsx_A = ones(nte, nx, np);
        modsx_B = zeros(nte, nx, np);
        
        % which events were already in the list upon last link?
        if ( isfield( ar.model(m).condition(c), 'modt' ) )
            % preserve the old values for the modification parameters
            mod_id = ismember( ar.model(m).condition(c).tEvents, ...
                ar.model(m).condition(c).modt );
            
            if ( ~isempty( mod_id ) && sum( mod_id ) > 0 )
                modx_A(mod_id,:) = ar.model(m).condition(c).modx_A;
                modx_B(mod_id,:) = ar.model(m).condition(c).modx_B;
                modsx_A(mod_id,:,:) = ar.model(m).condition(c).modsx_A;
                modsx_B(mod_id,:,:) = ar.model(m).condition(c).modsx_B;
            end
        end
        
        ar.model(m).condition(c).modt = ar.model(m).condition(c).tEvents;
        ar.model(m).condition(c).modx_A = modx_A;
        ar.model(m).condition(c).modx_B = modx_B;
        ar.model(m).condition(c).modsx_A = modsx_A;
        ar.model(m).condition(c).modsx_B = modsx_B;
        
        % Equilibration time
        ar.model(m).condition(c).tEq = NaN;
        
        % steady state sensitivities
        ar.model(m).condition(c).qSteadyState = false(1,nx);
        ar.model(m).condition(c).ssRelative = false(1,nx);
        ar.model(m).condition(c).dxdt = zeros(1, nx);
        ar.model(m).condition(c).ddxdtdp = zeros(nx, np);
        ar.model(m).condition(c).stdSteadyState = zeros(1,nx) + ar.config.steady_state_constraint;
        
        ar.model(m).condition(c).start = 0;
        ar.model(m).condition(c).stop = 0;
        ar.model(m).condition(c).stop_data = 0;
        
        % conditions with events
        if(~isempty(ar.model(m).condition(c).tEvents))
            useEvents = 1;
            ar.model(m).condition(c).qEvents = 1;
        else
            ar.model(m).condition(c).qEvents = 0;
        end
    end
end

if(~silent)
    arFprintf(2, 'linking parameters...\n');
end

% copy old values for event and multiple shooting settings if they exist
if (isfield(ar, 'config'))
    if (isfield(ar.config, 'useMS'))
        useMS = ar.config.useMS;
    end
    if (isfield(ar.config, 'useEvents'))
        useEvents = ar.config.useEvents;
    end
end
ar.config.useEvents = useEvents;
ar.config.useMS 	= useMS;

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
else 
    p = [];
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
                %fystd_regexp= regexp(ar.model(m).data(d).fystd,'sd_[a-zA-Z_0-9]+','match');
                %qerr =  ismember(ar.pLabel,[fystd_regexp{:}]); %R2013a compatible
                qerr = ismember(ar.pLabel, ar.model(m).data(d).pystd);
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

% Link time points intercondition constraints
if ( isfield( ar, 'conditionconstraints' ) )
    for jC = 1 : length( ar.conditionconstraints )
        m1 = ar.conditionconstraints(jC).m1;
        m2 = ar.conditionconstraints(jC).m2;
        c1 = ar.conditionconstraints(jC).c1;
        c2 = ar.conditionconstraints(jC).c2;
        t  = ar.conditionconstraints(jC).t;
        ar.conditionconstraints(jC).tLink1 = find( ismember( ar.model(m1).condition(c1).tExp, t ) );
        ar.conditionconstraints(jC).tLink2 = find( ismember( ar.model(m2).condition(c2).tExp, t ) );
    end
end

% populate threads
populate_threads( 'threads', 'condition', 'nTasks');
populate_threads( 'ss_threads', 'ss_condition', 'nssTasks');
ar.config.nThreads = length(ar.config.threads);
ar.config.nssThreads = length(ar.config.ss_threads);

% reset values
if(exist('plabel','var'))
    arSetPars(plabel, p, qfit, qlog10, lb, ub, type, meanp, stdp);
end

% plotting
for jm=1:length(ar.model)
    if(~isfield(ar.model(jm), 'qPlotYs') || isempty(ar.model(jm).qPlotYs))
        if(length(ar.model(jm).plot) > 10)
            if(~silent)
                arFprintf(1, 'Automatic plotting disabled for model %i. Please use arPlotter for plotting.\n', jm);
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
    
    if ~isempty(p) && sum(  abs( p-ar.p )>1e-10 & ar.p~=-1)>0
        warning('arLink overwrites parameter values by the default ar.pExtern defined in the def file.');
    end
end



% assigning first to cells 'tmp' and then in an array of SIMILAR structs is
% required for serveral models/several data sets/ several plots:
ar = orderfields(ar);
tmp = cell(size(ar.model));
for m=1:length(ar.model)
    tmp{m} = orderfields(ar.model(m));
end
ar.model = [tmp{:}];

for m=1:length(ar.model)
    if(isfield(ar.model(m),'data'))
        tmp = cell(size(ar.model(m).data));
        for d=1:length(ar.model(m).data)
            tmp{d} = orderfields(ar.model(m).data(d));
        end
        ar.model(m).data = [tmp{:}];
    end

    if(isfield(ar.model(m),'plot'))
        tmp = cell(size(ar.model(m).plot));
        for p=1:length(ar.model(m).plot)
            tmp{p} = orderfields(ar.model(m).plot(p));
        end
        ar.model(m).plot = [tmp{:}];
    end

end


function populate_threads( thread_fieldname, condition_fieldname, ntask_fieldname)

global ar;

% populate threads
ar.config.(thread_fieldname) = [];
ar.config.(thread_fieldname)(1).id = 0;
ar.config.(thread_fieldname)(1).n = 0;
ar.config.(thread_fieldname)(1).nd = 0;
ar.config.(thread_fieldname)(1).ms = int32([]);
ar.config.(thread_fieldname)(1).cs = int32([]);
ithread = 1;
ar.config.(ntask_fieldname) = 0;
for m = 1:length(ar.model)
    if (isfield(ar.model(m), condition_fieldname))
        for c = 1:length(ar.model(m).(condition_fieldname) )
            if(length(ar.config.(thread_fieldname))<ithread)
                ar.config.(thread_fieldname)(ithread).id = ithread-1;
                ar.config.(thread_fieldname)(ithread).n = 0;
                ar.config.(thread_fieldname)(ithread).nd = 0;
                ar.config.(thread_fieldname)(ithread).ms = int32([]);
                ar.config.(thread_fieldname)(ithread).cs = int32([]);
            end
            ar.config.(ntask_fieldname) = ar.config.(ntask_fieldname) + 1;
            ar.config.(thread_fieldname)(ithread).n = ...
                ar.config.(thread_fieldname)(ithread).n + 1;
            ar.config.(thread_fieldname)(ithread).nd = ...
                ar.config.(thread_fieldname)(ithread).nd + ...
                length(ar.model(m).(condition_fieldname)(c).dLink);
            ar.config.(thread_fieldname)(ithread).ms(end+1) = int32(m-1);
            ar.config.(thread_fieldname)(ithread).cs(end+1) = int32(c-1);
            ithread = ithread + 1;
            if(ithread>ar.config.nParallel)
                ithread = 1;
            end
        end
    end
end

% Invalidate cache so simulations do not get skipped
arCheckCache(1);