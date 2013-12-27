% Parse Model and write function files
%
% arParseModel(forceParsing)
%   forceParsing:                                   [false]
%
% Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)

function arParseModel(forceParsing)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('forceParsing','var'))
    forceParsing = false;
end

checksum_global = addToCheckSum(ar.info.c_version_code);
for m=1:length(ar.model)
    fprintf('\n');
    
    % parse model
    arParseODE(m);
    
    % extract & parse conditions
    ar.model(m).condition = [];
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            
            % conditions checksum
            qdynparas = ismember(ar.model(m).data(d).p, ar.model(m).px) | ... %R2013a compatible
                ismember(ar.model(m).data(d).p, ar.model(m).data(d).pu); %R2013a compatible
            
            checksum_cond = addToCheckSum(ar.model(m).data(d).fu);
            checksum_cond = addToCheckSum(ar.model(m).p, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).fv, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).N, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).cLink, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).data(d).fp(qdynparas), checksum_cond);
            checkstr_cond = getCheckStr(checksum_cond);
            
            % data checksum
            checksum_data = addToCheckSum(ar.model(m).data(d).fu);
            checksum_data = addToCheckSum(ar.model(m).data(d).p, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fy, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fystd, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fp, checksum_data);
            checkstr_data = getCheckStr(checksum_data);
            
            ar.model(m).data(d).checkstr = checkstr_data;
            ar.model(m).data(d).fkt = [ar.model(m).data(d).name '_' checkstr_data];
            
            cindex = -1;
            for c=1:length(ar.model(m).condition)
                if(strcmp(checkstr_cond, ar.model(m).condition(c).checkstr))
                    cindex = c;
                end
            end
            
            % global checksum
            if(isempty(checksum_global))
                checksum_global = addToCheckSum(ar.model(m).data(d).fkt);
            else
                checksum_global = addToCheckSum(ar.model(m).data(d).fkt, checksum_global);
            end
            
            if(cindex == -1) % append new condition
                cindex = length(ar.model(m).condition) + 1;
                
                ar.model(m).condition(cindex).status = 0;
                
                ar.model(m).condition(cindex).fu = ar.model(m).data(d).fu;
                ar.model(m).condition(cindex).fp = ar.model(m).data(d).fp(qdynparas);
                ar.model(m).condition(cindex).p = ar.model(m).data(d).p(qdynparas);
                
                ar.model(m).condition(cindex).checkstr = checkstr_cond;
                ar.model(m).condition(cindex).fkt = [ar.model(m).name '_' checkstr_cond];
                
                ar.model(m).condition(cindex).dLink = d;
                
                % global checksum
                checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt, checksum_global);
                
                % link data to condition
                ar.model(m).data(d).cLink = length(ar.model(m).condition);
                
                % for multiple shooting
                if(isfield(ar.model(m).data(d), 'ms_index') && ~isempty(ar.model(m).data(d).ms_index))
                    ar.model(m).condition(cindex).ms_index = ...
                        ar.model(m).data(d).ms_index;
                    ar.model(m).condition(cindex).ms_snip_index = ...
                        ar.model(m).data(d).ms_snip_index;
                    ar.model(m).condition(cindex).ms_snip_start = ar.model(m).data(d).tLim(1);
                end
            else
                % link data to condition
                ar.model(m).condition(cindex).dLink(end+1) = d;
                ar.model(m).data(d).cLink = cindex;
                
                % for multiple shooting
                if(isfield(ar.model(m).data(d), 'ms_index') && ~isempty(ar.model(m).data(d).ms_index))
                    ar.model(m).condition(cindex).ms_index(end+1) = ...
                        ar.model(m).data(d).ms_index;
                    ar.model(m).condition(cindex).ms_snip_index(end+1) = ...
                        ar.model(m).data(d).ms_snip_index;
                    ar.model(m).condition(cindex).ms_snip_start(end+1) = ar.model(m).data(d).tLim(1);
                end
            end
        end
        
        % skip parse conditions
        doskip = nan(1,length(ar.model(m).condition));
        for c=1:length(ar.model(m).condition)
            doskip(c) = ~forceParsing && exist(['./Compiled/' ar.info.c_version_code '/' ar.model(m).condition(c).fkt '.c'],'file');
        end
        
        % parse conditions
        config = ar.config;
        model.name = ar.model(m).name;
        model.fv = ar.model(m).fv;
        model.px0 = ar.model(m).px0;
        model.sym = ar.model(m).sym;
        model.t = ar.model(m).t;
        model.x = ar.model(m).x;
        model.u = ar.model(m).u;
        model.us = ar.model(m).us;
        model.xs = ar.model(m).xs;
        model.vs = ar.model(m).vs;
        model.N = ar.model(m).N;
        condition = ar.model(m).condition;
        parfor c=1:length(ar.model(m).condition)
             condition_new(c) = arParseCondition(config, model, condition(c), m, c, doskip(c));
        end
        ar.model(m).condition = condition_new;
        
        % skip parse data
        doskip = nan(1,length(ar.model(m).data));
        for d=1:length(ar.model(m).data)            
            doskip(d) = ~forceParsing && exist(['./Compiled/' ar.info.c_version_code '/' ar.model(m).data(d).fkt '.c'],'file');
        end
        
        % parse data
        data = ar.model(m).data;
        parfor d=1:length(ar.model(m).data)
            c = data(d).cLink;
            data_new(d) = arParseOBS(config, model, condition_new(c), data(d), m, c, d, doskip(d));
        end
        ar.model(m).data = data_new;
    else
        qdynparas = ismember(ar.model(m).p, ar.model(m).px) | ... %R2013a compatible
            ismember(ar.model(m).p, ar.model(m).pu); %R2013a compatible
        
        % conditions checksum
        checksum_cond = addToCheckSum(ar.model(m).fu);
        checksum_cond = addToCheckSum(ar.model(m).p(qdynparas), checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fv, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).N, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).cLink, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fp, checksum_cond);
        
        % append condition
        cindex = 1;
        
        ar.model(m).condition(cindex).status = 0;
        ar.model(m).condition(cindex).fu = ar.model(m).fu;
        ar.model(m).condition(cindex).fp = ar.model(m).fp(qdynparas);
        ar.model(m).condition(cindex).p = ar.model(m).p(qdynparas);
        ar.model(m).condition(cindex).checkstr = getCheckStr(checksum_cond);
        ar.model(m).condition(cindex).fkt = [ar.model(m).name '_' ar.model(m).condition(cindex).checkstr];
        ar.model(m).condition(cindex).dLink = [];
        
        % global checksum
        if(isempty(checksum_global))
            checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt);
        else
            checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt, checksum_global);
        end
        
        % skip parse conditions
        doskip = nan(1,length(ar.model(m).condition));
        for c=1:length(ar.model(m).condition)
            doskip(c) = ~forceParsing && exist(['./Compiled/' ar.info.c_version_code '/' ar.model(m).condition(c).fkt '.c'],'file');
        end
        
        % parse conditions
        config = ar.config;
        model.name = ar.model(m).name;
        model.fv = ar.model(m).fv;
        model.px0 = ar.model(m).px0;
        model.sym = ar.model(m).sym;
        model.t = ar.model(m).t;
        model.x = ar.model(m).x;
        model.u = ar.model(m).u;
        model.us = ar.model(m).us;
        model.xs = ar.model(m).xs;
        model.vs = ar.model(m).vs;
        model.N = ar.model(m).N;
        condition = ar.model(m).condition;
        parfor c=1:length(ar.model(m).condition)
             condition_new(c) = arParseCondition(config, model, condition(c), m, c, doskip(c));
        end
        ar.model(m).condition = condition_new;
        
        % plot setup
        if(~isfield(ar.model(m), 'plot'))
            ar.model(m).plot(1).name = ar.model(m).name;
        else
            ar.model(m).plot(end+1).name = ar.model(m).name;
        end
        ar.model(m).plot(end).doseresponse = false;
        ar.model(m).plot(end).dLink = 0;
        ar.model(m).plot(end).ny = 0;
        ar.model(m).plot(end).condition = {};
    end
end

ar.checkstr = getCheckStr(checksum_global);
ar.fkt = ['arSimuCalcFun_' ar.checkstr];



% ODE
function arParseODE(m)
global ar

fprintf('parsing model m%i, %s...', m, ar.model(m).name);

% make short strings
ar.model(m).xs = {};
ar.model(m).us = {};
ar.model(m).vs = {};

for j=1:length(ar.model(m).x)
    ar.model(m).xs{j} = sprintf('x[%i]',j);
end
fprintf('x=%i, ', length(ar.model(m).xs));
for j=1:length(ar.model(m).u)
    ar.model(m).us{j} = sprintf('u[%i]',j);
end
fprintf('u=%i, ', length(ar.model(m).u));
for j=1:length(ar.model(m).fv)
    ar.model(m).vs{j} = sprintf('v[%i]',j);
end
fprintf('v=%i, ', length(ar.model(m).fv));

% make syms
ar.model(m).sym.x = sym(ar.model(m).x);
ar.model(m).sym.xs = sym(ar.model(m).xs);
ar.model(m).sym.px0 = sym(ar.model(m).px0);
ar.model(m).sym.u = sym(ar.model(m).u);
ar.model(m).sym.us = sym(ar.model(m).us);
ar.model(m).sym.vs = sym(ar.model(m).vs);
ar.model(m).sym.fv = sym(ar.model(m).fv);

% compartment volumes
if(~isempty(ar.model(m).pc)) 
    % make syms
    ar.model(m).sym.pc = sym(ar.model(m).pc);
    ar.model(m).sym.C = sym(ones(size(ar.model(m).N)));
    
    if(~isfield(ar.model(m),'isAmountBased') || ~ar.model(m).isAmountBased)
        for j=1:size(ar.model(m).N,1) % for every species j
            qinfluxwitheducts = ar.model(m).N(j,:) > 0 & sum(ar.model(m).N < 0,1) > 0;
            eductcompartment = zeros(size(qinfluxwitheducts));
            for jj=find(qinfluxwitheducts)
				eductcompartment(jj) = unique(ar.model(m).cLink(ar.model(m).N(:,jj)<0)); %R2013a compatible
            end
            
            cfaktor = sym(ones(size(qinfluxwitheducts)));
            for jj=find(qinfluxwitheducts & eductcompartment~=ar.model(m).cLink(j))
                cfaktor(jj) = ar.model(m).sym.pc(eductcompartment(jj)) / ...
                    ar.model(m).sym.pc(ar.model(m).cLink(j));
            end
            ar.model(m).sym.C(j,:) = transpose(cfaktor);
        end
    else
        for j=1:size(ar.model(m).N,1) % for every species j
            ar.model(m).sym.C(j,:) = ar.model(m).sym.C(j,:) / ar.model(m).sym.pc(ar.model(m).cLink(j));
        end
    end
else
    ar.model(m).sym.C = sym(ones(size(ar.model(m).N)));
end

% derivatives
if(~isempty(ar.model(m).sym.fv))
    ar.model(m).sym.dfvdx = jacobian(ar.model(m).sym.fv, ar.model(m).sym.x);
    if(~isempty(ar.model(m).sym.us))
        ar.model(m).sym.dfvdu = jacobian(ar.model(m).sym.fv, ar.model(m).sym.u);
    else
        ar.model(m).sym.dfvdu = sym(ones(length(ar.model(m).sym.fv), 0));
    end
else
    ar.model(m).sym.dfvdx = sym(ones(0, length(ar.model(m).sym.x)));
    ar.model(m).sym.dfvdu = sym(ones(0, length(ar.model(m).sym.u)));
end

ar.model(m).qdvdx_nonzero = logical(ar.model(m).sym.dfvdx~=0);
ar.model(m).qdvdu_nonzero = logical(ar.model(m).sym.dfvdu~=0);

tmpsym = ar.model(m).sym.dfvdx;
tmpsym = mysubs(tmpsym, ar.model(m).sym.x, ones(size(ar.model(m).sym.x))/2);
tmpsym = mysubs(tmpsym, ar.model(m).sym.u, ones(size(ar.model(m).sym.u))/2);
tmpsym = mysubs(tmpsym, sym(ar.model(m).p), ones(size(ar.model(m).p))/2);

ar.model(m).qdvdx_negative = double(tmpsym) < 0;

tmpsym = ar.model(m).sym.dfvdu;
tmpsym = mysubs(tmpsym, ar.model(m).sym.x, ones(size(ar.model(m).sym.x))/2);
tmpsym = mysubs(tmpsym, ar.model(m).sym.u, ones(size(ar.model(m).sym.u))/2);
tmpsym = mysubs(tmpsym, sym(ar.model(m).p), ones(size(ar.model(m).p))/2);

ar.model(m).qdvdu_negative = double(tmpsym) < 0;

fprintf('done\n');



% Condition
function condition = arParseCondition(config, model, condition, m, c, doskip)

fprintf('parsing condition m%i c%i, %s (%s)...', m, c, model.name, condition.checkstr);

% hard code conditions
condition.sym.p = sym(condition.p);
condition.sym.fp = sym(condition.fp);
condition.sym.fpx0 = sym(model.px0);
condition.sym.fpx0 = mysubs(condition.sym.fpx0, condition.sym.p, condition.sym.fp);
condition.sym.fv = sym(model.fv);
condition.sym.fv = mysubs(condition.sym.fv, condition.sym.p, condition.sym.fp);
condition.sym.fu = sym(condition.fu);
condition.sym.fu = mysubs(condition.sym.fu, condition.sym.p, condition.sym.fp);
condition.sym.C = mysubs(model.sym.C, condition.sym.p, condition.sym.fp);

% predictor
condition.sym.fv = mysubs(condition.sym.fv, sym(model.t), sym('t'));
condition.sym.fu = mysubs(condition.sym.fu, sym(model.t), sym('t'));

% remaining initial conditions
qinitial = ismember(condition.p, model.px0); %R2013a compatible

varlist = cellfun(@symvar, condition.fp(qinitial), 'UniformOutput', false);
condition.px0 = union(vertcat(varlist{:}), [])'; %R2013a compatible

% remaining parameters
varlist = cellfun(@symvar, condition.fp, 'UniformOutput', false);
condition.pold = condition.p;
condition.p = setdiff(setdiff(union(vertcat(varlist{:}), [])', model.x), model.u); %R2013a compatible

if(doskip)
    fprintf('skipped\n');
    
    condition.ps = {};
    condition.qfu_nonzero = [];
    condition.qdvdx_nonzero = [];
    condition.qdvdu_nonzero = [];
    condition.qdvdp_nonzero = [];
    condition.dvdx = {};
    condition.dvdu = {};
    condition.dvdp = {};
    condition.qdfxdx_nonzero = [];
    condition.dfxdx = {};
    condition.su = {};
    condition.sx = {};
    condition.qfsv_nonzero = [];
    condition.sv = {};
    
    return;
end

% make short strings
condition.ps = {};
for j=1:length(condition.p)
    condition.ps{j} = sprintf('p[%i]',j);
end
fprintf('p=%i, ', length(condition.p));

% make syms
condition.sym.p = sym(condition.p);
condition.sym.ps = sym(condition.ps);
condition.sym.px0s = mysubs(sym(condition.px0), ...
    condition.sym.p, condition.sym.ps);

% make syms
condition.sym.fv = mysubs(condition.sym.fv, model.sym.x, model.sym.xs);
condition.sym.fv = mysubs(condition.sym.fv, model.sym.u, model.sym.us);

condition.sym.fv = mysubs(condition.sym.fv, condition.sym.p, condition.sym.ps);
condition.sym.fu = mysubs(condition.sym.fu, condition.sym.p, condition.sym.ps);
condition.sym.fpx0 = mysubs(condition.sym.fpx0, condition.sym.p, condition.sym.ps);

% remove zero inputs
condition.qfu_nonzero = logical(condition.sym.fu ~= 0);
if(~isempty(model.sym.us))
    condition.sym.fv = mysubs(condition.sym.fv, model.sym.us(~condition.qfu_nonzero), ...
        sym(zeros(1,sum(~condition.qfu_nonzero))));
end

% derivatives
if(~isempty(condition.sym.fv))
    condition.sym.dfvdx = jacobian(condition.sym.fv, model.sym.xs);
    if(~isempty(model.sym.us))
        condition.sym.dfvdu = jacobian(condition.sym.fv, model.sym.us);
    else
        condition.sym.dfvdu = sym(ones(length(condition.sym.fv), 0));
    end
    condition.sym.dfvdp = jacobian(condition.sym.fv, condition.sym.ps);
else
    condition.sym.dfvdx = sym(ones(0, length(model.sym.xs)));
    condition.sym.dfvdu = sym(ones(0, length(model.sym.us)));
    condition.sym.dfvdp = sym(ones(0, length(condition.sym.ps)));
end

% flux signs
condition.qdvdx_nonzero = logical(condition.sym.dfvdx~=0);
condition.qdvdu_nonzero = logical(condition.sym.dfvdu~=0);
condition.qdvdp_nonzero = logical(condition.sym.dfvdp~=0);

% short terms
condition.dvdx = cell(length(model.vs), length(model.xs));
for j=1:length(model.vs)
    for i=1:length(model.xs)
        if(condition.qdvdx_nonzero(j,i))
            condition.dvdx{j,i} = sprintf('dvdx[%i]', j + (i-1)*length(model.vs));
        else
            condition.dvdx{j,i} = '0';
        end
    end
end
condition.sym.dvdx = sym(condition.dvdx);
fprintf('dvdx=%i, ', sum(condition.qdvdx_nonzero(:)));

condition.dvdu = cell(length(model.vs), length(model.us));
for j=1:length(model.vs)
    for i=1:length(model.us)
        if(condition.qdvdu_nonzero(j,i))
            condition.dvdu{j,i} = sprintf('dvdu[%i]', j + (i-1)*length(model.vs));
        else
            condition.dvdu{j,i} = '0';
        end
    end
end
condition.sym.dvdu = sym(condition.dvdu);
fprintf('dvdu=%i, ', sum(condition.qdvdu_nonzero(:)));

condition.dvdp = cell(length(model.vs), length(condition.ps));
for j=1:length(model.vs)
    for i=1:length(condition.ps)
        if(condition.qdvdp_nonzero(j,i))
            condition.dvdp{j,i} = sprintf('dvdp[%i]', j + (i-1)*length(model.vs));
        else
            condition.dvdp{j,i} = '0';
        end
    end
end
condition.sym.dvdp = sym(condition.dvdp);
fprintf('dvdp=%i, ', sum(condition.qdvdp_nonzero(:)));

% make equations
condition.sym.C = mysubs(condition.sym.C, condition.sym.p, condition.sym.ps);
condition.sym.fx = (model.N .* condition.sym.C) * transpose(model.sym.vs);

% Jacobian dfxdx
if(config.useJacobian)
    condition.sym.dfxdx = (model.N .* condition.sym.C) * condition.sym.dvdx;
    condition.qdfxdx_nonzero = logical(condition.sym.dfxdx~=0);
    for j=1:length(model.xs)
        for i=1:length(model.xs)
            if(condition.qdfxdx_nonzero(j,i))
                condition.dfxdx{j,i} = sprintf('dfxdx[%i]', j + (i-1)*length(model.xs));
            else
                condition.dfxdx{j,i} = '0';
            end
        end
    end
    fprintf('dfxdx=%i, ', sum(condition.qdfxdx_nonzero(:)));
end

% sx sensitivities
if(config.useSensis)
	% su
    condition.su = cell(length(model.us), length(condition.ps));
    for j=1:length(model.us)
        for i=1:length(condition.ps)
            if(condition.qfu_nonzero(j))
                condition.su{j,i} = sprintf('su[%i]', j + (i-1)*length(model.us));
            else
                condition.su{j,i} = '0';
            end
        end
    end
    condition.sym.su = sym(condition.su);
    fprintf('su=%i, ', length(condition.ps)*sum(condition.qfu_nonzero(:)));
    
    % input derivatives 
    if(~isempty(condition.sym.ps))
        if(~isempty(condition.sym.fu))
        condition.sym.dfudp = ...
            jacobian(condition.sym.fu, condition.sym.ps);
        else
            condition.sym.dfudp = sym(ones(0,length(condition.sym.ps)));
        end
        % derivatives of step1 (DISABLED)
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'step1('))
                condition.sym.dfudp(j,:) = 0;
            end
        end
        
        % derivatives of step2 (DISABLED)
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'step2('))
                condition.sym.dfudp(j,:) = 0;
            end
        end
        
        % derivatives of spline3
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline3('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline3)('))
                        ustr = strrep(ustr, 'D([3], spline3)(', 'Dspline3(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline3)('))
                        ustr = strrep(ustr, 'D([5], spline3)(', 'Dspline3(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline3)('))
                        ustr = strrep(ustr, 'D([7], spline3)(', 'Dspline3(');
                        ustr = strrep(ustr, ')', ', 3)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end
        
        % derivatives of spline_pos3
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline_pos3('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline_pos3)('))
                        ustr = strrep(ustr, 'D([3], spline_pos3)(', 'Dspline_pos3(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline_pos3)('))
                        ustr = strrep(ustr, 'D([5], spline_pos3)(', 'Dspline_pos3(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline_pos3)('))
                        ustr = strrep(ustr, 'D([7], spline_pos3)(', 'Dspline_pos3(');
                        ustr = strrep(ustr, ')', ', 3)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end
        
        % derivatives of spline4
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline4('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline4)('))
                        ustr = strrep(ustr, 'D([3], spline4)(', 'Dspline4(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline4)('))
                        ustr = strrep(ustr, 'D([5], spline4)(', 'Dspline4(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline4)('))
                        ustr = strrep(ustr, 'D([7], spline4)(', 'Dspline4(');
                        ustr = strrep(ustr, ')', ', 3)');
                    elseif(strfind(ustr, 'D([9], spline4)('))
                        ustr = strrep(ustr, 'D([9], spline4)(', 'Dspline4(');
                        ustr = strrep(ustr, ')', ', 4)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end
        
        % derivatives of spline_pos4
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline_pos4('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline_pos4)('))
                        ustr = strrep(ustr, 'D([3], spline_pos4)(', 'Dspline_pos4(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline_pos4)('))
                        ustr = strrep(ustr, 'D([5], spline_pos4)(', 'Dspline_pos4(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline_pos4)('))
                        ustr = strrep(ustr, 'D([7], spline_pos4)(', 'Dspline_pos4(');
                        ustr = strrep(ustr, ')', ', 3)');
                    elseif(strfind(ustr, 'D([9], spline_pos4)('))
                        ustr = strrep(ustr, 'D([9], spline_pos4)(', 'Dspline_pos4(');
                        ustr = strrep(ustr, ')', ', 4)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end
        
        % derivatives of spline5
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline5('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline5)('))
                        ustr = strrep(ustr, 'D([3], spline5)(', 'Dspline5(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline5)('))
                        ustr = strrep(ustr, 'D([5], spline5)(', 'Dspline5(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline5)('))
                        ustr = strrep(ustr, 'D([7], spline5)(', 'Dspline5(');
                        ustr = strrep(ustr, ')', ', 3)');
                    elseif(strfind(ustr, 'D([9], spline5)('))
                        ustr = strrep(ustr, 'D([9], spline5)(', 'Dspline5(');
                        ustr = strrep(ustr, ')', ', 4)');
                    elseif(strfind(ustr, 'D([11], spline5)('))
                        ustr = strrep(ustr, 'D([11], spline5)(', 'Dspline5(');
                        ustr = strrep(ustr, ')', ', 5)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end
        
        % derivatives of spline_pos5
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline_pos5('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline_pos5)('))
                        ustr = strrep(ustr, 'D([3], spline_pos5)(', 'Dspline_pos5(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline_pos5)('))
                        ustr = strrep(ustr, 'D([5], spline_pos5)(', 'Dspline_pos5(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline_pos5)('))
                        ustr = strrep(ustr, 'D([7], spline_pos5)(', 'Dspline_pos5(');
                        ustr = strrep(ustr, ')', ', 3)');
                    elseif(strfind(ustr, 'D([9], spline_pos5)('))
                        ustr = strrep(ustr, 'D([9], spline_pos5)(', 'Dspline_pos5(');
                        ustr = strrep(ustr, ')', ', 4)');
                    elseif(strfind(ustr, 'D([11], spline_pos5)('))
                        ustr = strrep(ustr, 'D([11], spline_pos5)(', 'Dspline_pos5(');
                        ustr = strrep(ustr, ')', ', 5)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end
        
          % derivatives of spline10
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline10('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline10)('))
                        ustr = strrep(ustr, 'D([3], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline10)('))
                        ustr = strrep(ustr, 'D([5], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline10)('))
                        ustr = strrep(ustr, 'D([7], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 3)');
                    elseif(strfind(ustr, 'D([9], spline10)('))
                        ustr = strrep(ustr, 'D([9], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 4)');
                    elseif(strfind(ustr, 'D([11], spline10)('))
                        ustr = strrep(ustr, 'D([11], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 5)');
                    elseif(strfind(ustr, 'D([13], spline10)('))
                        ustr = strrep(ustr, 'D([13], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 6)');
                    elseif(strfind(ustr, 'D([15], spline10)('))
                        ustr = strrep(ustr, 'D([15], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 7)');
                    elseif(strfind(ustr, 'D([17], spline10)('))
                        ustr = strrep(ustr, 'D([17], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 8)');
                    elseif(strfind(ustr, 'D([19], spline10)('))
                        ustr = strrep(ustr, 'D([19], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 9)');
                    elseif(strfind(ustr, 'D([21], spline10)('))
                        ustr = strrep(ustr, 'D([21], spline10)(', 'Dspline10(');
                        ustr = strrep(ustr, ')', ', 10)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end
        
        % derivatives of spline_pos10
        for j=1:length(model.u)
            if(strfind(condition.fu{j}, 'spline_pos10('))
                for j2=1:length(condition.sym.dfudp(j,:))
                    ustr = char(condition.sym.dfudp(j,j2));
                    if(strfind(ustr, 'D([3], spline_pos10)('))
                        ustr = strrep(ustr, 'D([3], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 1)');
                    elseif(strfind(ustr, 'D([5], spline_pos10)('))
                        ustr = strrep(ustr, 'D([5], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 2)');
                    elseif(strfind(ustr, 'D([7], spline_pos10)('))
                        ustr = strrep(ustr, 'D([7], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 3)');
                    elseif(strfind(ustr, 'D([9], spline_pos10)('))
                        ustr = strrep(ustr, 'D([9], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 4)');
                    elseif(strfind(ustr, 'D([11], spline_pos10)('))
                        ustr = strrep(ustr, 'D([11], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 5)');
                    elseif(strfind(ustr, 'D([13], spline_pos10)('))
                        ustr = strrep(ustr, 'D([13], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 6)');
                    elseif(strfind(ustr, 'D([15], spline_pos10)('))
                        ustr = strrep(ustr, 'D([15], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 7)');
                    elseif(strfind(ustr, 'D([17], spline_pos10)('))
                        ustr = strrep(ustr, 'D([17], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 8)');
                    elseif(strfind(ustr, 'D([19], spline_pos10)('))
                        ustr = strrep(ustr, 'D([19], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 9)');
                    elseif(strfind(ustr, 'D([21], spline_pos10)('))
                        ustr = strrep(ustr, 'D([21], spline_pos10)(', 'Dspline_pos10(');
                        ustr = strrep(ustr, ')', ', 10)');
                    end
                    condition.sym.dfudp(j,j2) = sym(ustr);
                end
            end
        end 
    end
    
	% sx
    condition.sx = cell(length(model.xs), length(condition.ps));
    for j=1:length(model.xs)
        for i=1:length(condition.ps)
            condition.sx{j,i} = sprintf('sx[%i]', j);
        end
    end
	condition.sym.sx = sym(condition.sx);
    
    condition.sym.fsv = condition.sym.dvdx * condition.sym.sx + ...
        condition.sym.dvdu * condition.sym.su + condition.sym.dvdp;
    
	% sv
    condition.qfsv_nonzero = logical(condition.sym.fsv ~= 0);
    condition.sv = cell(length(model.vs), length(condition.ps));
    for j=1:length(model.vs)
        for i=1:length(condition.ps)
            if(condition.qfsv_nonzero(j,i))
                condition.sv{j,i} = sprintf('sv[%i]', j);
            else
                condition.sv{j,i} = '0';
            end
        end
    end
    condition.sym.sv = sym(condition.sv);
    fprintf('sv=%i, ', sum(condition.qfsv_nonzero(:)));
    
    if(config.useSensiRHS)
        condition.sym.fsx = (model.N .* condition.sym.C) * condition.sym.sv;
        fprintf('sx=%i... ', numel(condition.sym.fsx));
    else
        fprintf('sx=skipped ');
    end
    
    % sx initials
    if(~isempty(condition.sym.fpx0))
        condition.sym.fsx0 = jacobian(condition.sym.fpx0, condition.sym.ps);
    else
        condition.sym.fsx0 = sym(ones(0, length(condition.sym.ps)));
    end
    
    % steady state sensitivities
    condition.sym.dfxdp = (model.N .* condition.sym.C) * (condition.sym.dvdp + ...
        condition.sym.dvdx*condition.sym.fsx0 + ...
        condition.sym.dvdu * condition.sym.dfudp);
end

fprintf('done\n');



% OBS
function data = arParseOBS(config, model, condition, data, m, c, d, doskip)

fprintf('parsing data m%i d%i -> c%i, %s (%s)...', m, d, c, data.name, data.checkstr);

% hard code conditions
data.sym.p = sym(data.p);
data.sym.fp = sym(data.fp);
data.sym.fy = sym(data.fy);
data.sym.fy = mysubs(data.sym.fy, data.sym.p, data.sym.fp);
data.sym.fystd = sym(data.fystd);
data.sym.fystd = mysubs(data.sym.fystd, data.sym.p, data.sym.fp);

data.sym.fu = sym(condition.fu);
data.sym.fu = mysubs(data.sym.fu, data.sym.p, data.sym.fp);
data.qfu_nonzero = logical(data.sym.fu ~= 0);

% predictor
data.sym.fu = mysubs(data.sym.fu, sym(model.t), sym('t'));
data.sym.fy = mysubs(data.sym.fy, sym(model.t), sym('t'));
data.sym.fystd = mysubs(data.sym.fystd, sym(model.t), sym('t'));

% remaining parameters
varlist = cellfun(@symvar, data.fp, 'UniformOutput', false);
data.pold = data.p;
data.p = setdiff(setdiff(union(vertcat(varlist{:}), [])', model.x), model.u); %R2013a compatible

if(doskip)
    fprintf('skipped\n');
    
    data.ps = {};
    data.ys = [];
    data.qu_measured = [];
    data.qx_measured = [];
    data.dfydxnon0 = {};
    data.sx = {};
    data.sy = {};
    
    return;
end

% make short strings
for j=1:length(data.p)
    data.ps{j} = sprintf('p[%i]',j);
end
data.ys = {};
for j=1:length(data.y)
    data.ys{j} = sprintf('y[%i]',j);
end
fprintf('y=%i, p=%i, ', length(data.y), length(data.p));

% make syms
data.sym.p = sym(data.p);
data.sym.ps = sym(data.ps);
data.sym.y = sym(data.y);
data.sym.ys = sym(data.ys);

% substitute
data.sym.fy = mysubs(data.sym.fy, ...
    model.sym.x, model.sym.xs);
data.sym.fy = mysubs(data.sym.fy, ...
    model.sym.u, model.sym.us);
data.sym.fy = mysubs(data.sym.fy, ...
    data.sym.p, data.sym.ps);

data.sym.fystd = mysubs(data.sym.fystd, ...
    model.sym.x, model.sym.xs);
data.sym.fystd = mysubs(data.sym.fystd, ...
    model.sym.u, model.sym.us);
data.sym.fystd = mysubs(data.sym.fystd, ...
    data.sym.y, data.sym.ys);
data.sym.fystd = mysubs(data.sym.fystd, ...
    data.sym.p, data.sym.ps);

% remove zero inputs
% data.sym.fy = mysubs(data.sym.fy, model.sym.us(~condition.qfu_nonzero), ...
%     sym(zeros(1,sum(~condition.qfu_nonzero))));
% data.sym.fystd = mysubs(data.sym.fystd, model.sym.us(~condition.qfu_nonzero), ...
%     sym(zeros(1,sum(~condition.qfu_nonzero))));

% derivatives fy
if(~isempty(data.sym.fy))
    if(~isempty(model.sym.us))
        data.sym.dfydu = jacobian(data.sym.fy, model.sym.us);
    else
        data.sym.dfydu = sym(ones(length(data.y), 0));
    end
    if(~isempty(model.x))
        data.sym.dfydx = jacobian(data.sym.fy, model.sym.xs);
    else
        data.sym.dfydx = [];
    end
	data.sym.dfydp = jacobian(data.sym.fy, data.sym.ps);
else
	data.sym.dfydu = [];
	data.sym.dfydx = [];
	data.sym.dfydp = [];
end

% what is measured ?
data.qu_measured = sum(logical(data.sym.dfydu~=0),1)>0;
data.qx_measured = sum(logical(data.sym.dfydx~=0),1)>0;

% derivatives fystd
if(~isempty(data.sym.fystd))
    if(~isempty(model.sym.us))
        data.sym.dfystddu = jacobian(data.sym.fystd, model.sym.us);
    else
        data.sym.dfystddu = sym(ones(length(data.y), 0));
    end
    if(~isempty(model.x))
        data.sym.dfystddx = jacobian(data.sym.fystd, model.sym.xs);
    else
        data.sym.dfystddx = [];
    end
    data.sym.dfystddp = jacobian(data.sym.fystd, data.sym.ps);
    data.sym.dfystddy = jacobian(data.sym.fystd, data.sym.ys);
else
	data.sym.dfystddp = [];
	data.sym.dfystddy = [];
    data.sym.dfystddx = [];
end

% observed directly and exclusively
data.dfydxnon0 = logical(data.sym.dfydx ~= 0);

if(config.useSensis)
    % sx sensitivities
    data.sx = {};
    for j=1:length(model.xs)
        for i=1:length(condition.p)
            data.sx{j,i} = sprintf('sx[%i]', j + (i-1)*length(model.xs));
        end
    end
    data.sym.sx = sym(data.sx);
    
    % sy sensitivities
    data.sy = {};
    for j=1:length(data.sym.fy)
        for i=1:length(data.sym.ps)
            data.sy{j,i} = sprintf('sy[%i]', j + (i-1)*length(data.sym.fy));
        end
    end
	data.sym.sy = sym(data.sy);
    
    % calculate sy
    if(~isempty(data.sym.sy))
        data.sym.fsy = data.sym.dfydp;
        if(~isempty(condition.p))
            qdynpara = ismember(data.p, condition.p); %R2013a compatible
        else
            qdynpara = false(size(data.p));
        end
        
        if(~isempty(condition.p))
            tmpfsx = data.sym.dfydx * ...
                data.sym.sx;
            if(~isfield(condition, 'sym') || ~isfield(condition.sym, 'su'))
                % su
                condition.su = cell(length(model.us), length(condition.p));
                for j=1:length(model.us)
                    for i=1:length(condition.p)
                        if(data.qfu_nonzero(j))
                            condition.su{j,i} = sprintf('su[%i]', j + (i-1)*length(model.us));
                        else
                            condition.su{j,i} = '0';
                        end
                    end
                end
                condition.sym.su = sym(condition.su);
                fprintf('su=%i, ', length(condition.p)*sum(data.qfu_nonzero(:)));
            end
            tmpfsu = data.sym.dfydu * ...
                condition.sym.su;
            if(~isempty(model.x))
                data.sym.fsy(:,qdynpara) = data.sym.fsy(:,qdynpara) + tmpfsx + tmpfsu;
            else
                data.sym.fsy(:,qdynpara) = data.sym.fsy(:,qdynpara) + tmpfsu;
            end
        end
    else
        data.sym.fsy = [];
    end
    fprintf('sy=%i, ', sum(logical(data.sym.fsy(:)~=0)));
    
    % calculate systd
    if(~isempty(data.sym.sy))
        data.sym.fsystd = data.sym.dfystddp + ...
            data.sym.dfystddy * data.sym.sy;
        if(~isempty(condition.p))
            qdynpara = ismember(data.p, condition.p); %R2013a compatible
        else
            qdynpara = false(size(data.p));
        end
        
        if(~isempty(condition.p))
            tmpfsx = data.sym.dfystddx * ...
                data.sym.sx;
            if(~isfield(condition, 'sym') || ~isfield(condition.sym, 'su'))
                % su
                condition.su = cell(length(model.us), length(condition.p));
                for j=1:length(model.us)
                    for i=1:length(condition.p)
                        if(data.qfu_nonzero(j))
                            condition.su{j,i} = sprintf('su[%i]', j + (i-1)*length(model.us));
                        else
                            condition.su{j,i} = '0';
                        end
                    end
                end
                condition.sym.su = sym(condition.su);
                fprintf('su=%i, ', length(condition.p)*sum(data.qfu_nonzero(:)));
            end
            tmpfsu = data.sym.dfystddu * ...
                condition.sym.su;
            if(~isempty(model.x))
                data.sym.fsystd(:,qdynpara) = data.sym.fsystd(:,qdynpara) + tmpfsx + tmpfsu;
            else
                data.sym.fsystd(:,qdynpara) = data.sym.fsystd(:,qdynpara) + tmpfsu;
            end
        end
    else
        data.sym.fsystd = [];
    end
    fprintf('systd=%i, ', sum(logical(data.sym.fsystd(:)~=0)));
end

fprintf('done\n');



% better subs
function out = mysubs(in, old, new)
if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
    matVer = ver('MATLAB');
    if(str2double(matVer.Version)>=8.1)
        out = subs(in, old(:), new(:));
    else
        out = subs(in, old(:), new(:), 0);
    end
else
    out = in;
end

function checksum = addToCheckSum(str, checksum)
algs = {'MD2','MD5','SHA-1','SHA-256','SHA-384','SHA-512'};
if(nargin<2)
    checksum = java.security.MessageDigest.getInstance(algs{2});
end
if(iscell(str))
    for j=1:length(str)
        checksum = addToCheckSum(str{j}, checksum);
    end
else
    if(~isempty(str))
        checksum.update(uint8(str(:)));
    end
end

function checkstr = getCheckStr(checksum)
h = typecast(checksum.digest,'uint8');
checkstr = dec2hex(h)';
checkstr = checkstr(:)';

clear checksum
