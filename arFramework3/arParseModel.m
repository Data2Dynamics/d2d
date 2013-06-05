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
    % parse model
    arParseODE(m);
    
    % extract & parse conditions
    ar.model(m).condition = [];
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            
            % conditions checksum
            qdynparas = ismember(ar.model(m).data(d).p, ar.model(m).p) | ...
                ismember(ar.model(m).data(d).p, ar.model(m).data(d).pu);
            
            checksum_cond = addToCheckSum(ar.model(m).data(d).fu);
            checksum_cond = addToCheckSum(ar.model(m).p, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).fv, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).fxeq, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).N, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).cLink, checksum_cond);
            checksum_cond = addToCheckSum(ar.model(m).data(d).fp(qdynparas), checksum_cond);
            checkstr_cond = getCheckStr(checksum_cond);
            
            % data checksum
            checksum_data = addToCheckSum(ar.model(m).data(d).fu);
            checksum_data = addToCheckSum(ar.model(m).data(d).p, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fy, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fystd, checksum_data);
            checksum_data = addToCheckSum(ar.model(m).data(d).fxeq, checksum_data);
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
                
                ar.model(m).condition(cindex).fxeq = ar.model(m).data(d).fxeq;
                ar.model(m).condition(cindex).fu = ar.model(m).data(d).fu;
                ar.model(m).condition(cindex).fp = ar.model(m).data(d).fp(qdynparas);
                ar.model(m).condition(cindex).p = ar.model(m).data(d).p(qdynparas);
                
                ar.model(m).condition(cindex).checkstr = checkstr_cond;
                ar.model(m).condition(cindex).fkt = [ar.model(m).name '_' checkstr_cond];
                
                ar.model(m).condition(cindex).dLink = d;
                
                % global checksum
                checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt, checksum_global);
                
                doskip = ~forceParsing && exist(['./Compiled/' ar.model(m).condition(cindex).fkt '.c'],'file');
                arParseCondition(m, cindex, doskip);
                
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
        
        % parse data
        for d=1:length(ar.model(m).data)
            doskip = ~forceParsing && exist(['./Compiled/' ar.model(m).data(d).fkt '.c'],'file');
            arParseOBS(m, d, doskip);
        end
    else
        % conditions checksum
        checksum_cond = addToCheckSum(ar.model(m).fu);
        checksum_cond = addToCheckSum(ar.model(m).p, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fv, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fxeq, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).N, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).cLink, checksum_cond);
        checksum_cond = addToCheckSum(ar.model(m).fp, checksum_cond);
        
        % append new condition
        cindex = length(ar.model(m).condition) + 1;
        
        ar.model(m).condition(cindex).status = 0;
        
        ar.model(m).condition(cindex).fxeq = ar.model(m).fxeq;
        ar.model(m).condition(cindex).fu = ar.model(m).fu;
        ar.model(m).condition(cindex).fp = ar.model(m).fp;
        ar.model(m).condition(cindex).p = ar.model(m).p;
        
        ar.model(m).condition(cindex).checkstr = getCheckStr(checksum_cond);
        ar.model(m).condition(cindex).fkt = [ar.model(m).name '_' ar.model(m).condition(cindex).checkstr];
        
        ar.model(m).condition(cindex).dLink = [];
        
        % global checksum
        if(isempty(checksum_global))
            checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt);
        else
            checksum_global = addToCheckSum(ar.model(m).condition(cindex).fkt, checksum_global);
        end
        
        doskip = ~forceParsing && exist(['./Compiled/' ar.model(m).condition(cindex).fkt '.c'],'file');
        arParseCondition(m, cindex, doskip);
        
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
                eductcompartment(jj) = unique(ar.model(m).cLink(ar.model(m).N(:,jj)<0));
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
ar.model(m).sym.dfvdx = jacobian(ar.model(m).sym.fv, ar.model(m).sym.x);
if(~isempty(ar.model(m).sym.us))
    ar.model(m).sym.dfvdu = jacobian(ar.model(m).sym.fv, ar.model(m).sym.u);
else
    ar.model(m).sym.dfvdu = sym(ones(length(ar.model(m).sym.fv), 0));
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
function arParseCondition(m, c, doskip)
global ar

fprintf('parsing condition m%i c%i, %s (%s)...', m, c, ar.model(m).name, ar.model(m).condition(c).checkstr);

% hard code conditions
ar.model(m).condition(c).sym.p = sym(ar.model(m).condition(c).p);
ar.model(m).condition(c).sym.fp = sym(ar.model(m).condition(c).fp);
ar.model(m).condition(c).sym.fpx0 = sym(ar.model(m).px0);
ar.model(m).condition(c).sym.fpx0 = mysubs(ar.model(m).condition(c).sym.fpx0, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.fp);
ar.model(m).condition(c).sym.fv = sym(ar.model(m).fv);
ar.model(m).condition(c).sym.fv = mysubs(ar.model(m).condition(c).sym.fv, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.fp);
if(~isempty(ar.model(m).condition(c).fxeq))
    for j=1:length(ar.model(m).condition(c).fxeq)
        ar.model(m).condition(c).sym.fxeq(j) = sym(ar.model(m).condition(c).fxeq{j});
    end
else
    ar.model(m).condition(c).sym.fxeq = sym([]);
end
ar.model(m).condition(c).sym.fxeq = mysubs(ar.model(m).condition(c).sym.fxeq, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.fp);
ar.model(m).condition(c).sym.fu = sym(ar.model(m).condition(c).fu);
ar.model(m).condition(c).sym.fu = mysubs(ar.model(m).condition(c).sym.fu, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.fp);
ar.model(m).condition(c).sym.C = mysubs(ar.model(m).sym.C, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.fp);

% predictor
ar.model(m).condition(c).sym.fv = mysubs(ar.model(m).condition(c).sym.fv, sym(ar.model(m).t), sym('t'));
ar.model(m).condition(c).sym.fu = mysubs(ar.model(m).condition(c).sym.fu, sym(ar.model(m).t), sym('t'));

% remaining initial conditions
qinitial = ismember(ar.model(m).condition(c).p, ar.model(m).px0);
varlist = cellfun(@symvar, ar.model(m).condition(c).fp(qinitial), 'UniformOutput', false);
ar.model(m).condition(c).px0 = union(vertcat(varlist{:}), [])';

% remaining parameters
varlist = cellfun(@symvar, ar.model(m).condition(c).fp, 'UniformOutput', false);
ar.model(m).condition(c).pold = ar.model(m).condition(c).p;
ar.model(m).condition(c).p = union(vertcat(varlist{:}), [])';

if(doskip)
    fprintf('skipping\n');
    return;
end

% make short strings
ar.model(m).condition(c).ps = {};
for j=1:length(ar.model(m).condition(c).p)
    ar.model(m).condition(c).ps{j} = sprintf('p[%i]',j);
end
fprintf('p=%i, ', length(ar.model(m).condition(c).p));

% make syms
ar.model(m).condition(c).sym.p = sym(ar.model(m).condition(c).p);
ar.model(m).condition(c).sym.ps = sym(ar.model(m).condition(c).ps);
ar.model(m).condition(c).sym.px0s = mysubs(sym(ar.model(m).condition(c).px0), ...
    ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.ps);

% make syms
ar.model(m).condition(c).sym.fv = mysubs(ar.model(m).condition(c).sym.fv, ar.model(m).sym.x, ar.model(m).sym.xs);
ar.model(m).condition(c).sym.fv = mysubs(ar.model(m).condition(c).sym.fv, ar.model(m).sym.u, ar.model(m).sym.us);
ar.model(m).condition(c).sym.fv = mysubs(ar.model(m).condition(c).sym.fv, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.ps);

ar.model(m).condition(c).sym.fu = mysubs(ar.model(m).condition(c).sym.fu, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.ps);
ar.model(m).condition(c).sym.fpx0 = mysubs(ar.model(m).condition(c).sym.fpx0, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.ps);
ar.model(m).condition(c).sym.fxeq = mysubs(ar.model(m).condition(c).sym.fxeq, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.ps);

% remove zero inputs
ar.model(m).condition(c).qfu_nonzero = logical(ar.model(m).condition(c).sym.fu ~= 0);
if(~isempty(ar.model(m).sym.us))
    ar.model(m).condition(c).sym.fv = mysubs(ar.model(m).condition(c).sym.fv, ar.model(m).sym.us(~ar.model(m).condition(c).qfu_nonzero), ...
        sym(zeros(1,sum(~ar.model(m).condition(c).qfu_nonzero))));
end

% derivatives
ar.model(m).condition(c).sym.dfvdx = jacobian(ar.model(m).condition(c).sym.fv, ar.model(m).sym.xs);
if(~isempty(ar.model(m).sym.us))
    ar.model(m).condition(c).sym.dfvdu = jacobian(ar.model(m).condition(c).sym.fv, ar.model(m).sym.us);
else
    ar.model(m).condition(c).sym.dfvdu = sym(ones(length(ar.model(m).condition(c).sym.fv), 0));
end
ar.model(m).condition(c).sym.dfvdp = jacobian(ar.model(m).condition(c).sym.fv, ar.model(m).condition(c).sym.ps);

% flux signs
ar.model(m).condition(c).qdvdx_nonzero = logical(ar.model(m).condition(c).sym.dfvdx~=0);
ar.model(m).condition(c).qdvdu_nonzero = logical(ar.model(m).condition(c).sym.dfvdu~=0);
ar.model(m).condition(c).qdvdp_nonzero = logical(ar.model(m).condition(c).sym.dfvdp~=0);

% short terms
ar.model(m).condition(c).dvdx = cell(length(ar.model(m).vs), length(ar.model(m).xs));
for j=1:length(ar.model(m).vs)
    for i=1:length(ar.model(m).xs)
        if(ar.model(m).condition(c).qdvdx_nonzero(j,i))
            ar.model(m).condition(c).dvdx{j,i} = sprintf('dvdx[%i]', j + (i-1)*length(ar.model(m).vs));
        else
            ar.model(m).condition(c).dvdx{j,i} = '0';
        end
    end
end
ar.model(m).condition(c).sym.dvdx = sym(ar.model(m).condition(c).dvdx);
fprintf('dvdx=%i, ', sum(ar.model(m).condition(c).qdvdx_nonzero(:)));

ar.model(m).condition(c).dvdu = cell(length(ar.model(m).vs), length(ar.model(m).us));
for j=1:length(ar.model(m).vs)
    for i=1:length(ar.model(m).us)
        if(ar.model(m).condition(c).qdvdu_nonzero(j,i))
            ar.model(m).condition(c).dvdu{j,i} = sprintf('dvdu[%i]', j + (i-1)*length(ar.model(m).vs));
        else
            ar.model(m).condition(c).dvdu{j,i} = '0';
        end
    end
end
ar.model(m).condition(c).sym.dvdu = sym(ar.model(m).condition(c).dvdu);
fprintf('dvdu=%i, ', sum(ar.model(m).condition(c).qdvdu_nonzero(:)));

ar.model(m).condition(c).dvdp = cell(length(ar.model(m).vs), length(ar.model(m).condition(c).ps));
for j=1:length(ar.model(m).vs)
    for i=1:length(ar.model(m).condition(c).ps)
        if(ar.model(m).condition(c).qdvdp_nonzero(j,i))
            ar.model(m).condition(c).dvdp{j,i} = sprintf('dvdp[%i]', j + (i-1)*length(ar.model(m).vs));
        else
            ar.model(m).condition(c).dvdp{j,i} = '0';
        end
    end
end
ar.model(m).condition(c).sym.dvdp = sym(ar.model(m).condition(c).dvdp);
fprintf('dvdp=%i, ', sum(ar.model(m).condition(c).qdvdp_nonzero(:)));

% make equations
ar.model(m).condition(c).sym.C = mysubs(ar.model(m).condition(c).sym.C, ar.model(m).condition(c).sym.p, ar.model(m).condition(c).sym.ps);
ar.model(m).condition(c).sym.fx = (ar.model(m).N .* ar.model(m).condition(c).sym.C) * transpose(ar.model(m).sym.vs);

% Jacobian dfxdx
if(ar.config.useJacobian)
    ar.model(m).condition(c).sym.dfxdx = (ar.model(m).N .* ar.model(m).condition(c).sym.C) * ar.model(m).condition(c).sym.dvdx;
    ar.model(m).condition(c).qdfxdx_nonzero = logical(ar.model(m).condition(c).sym.dfxdx~=0);
    for j=1:length(ar.model(m).xs)
        for i=1:length(ar.model(m).xs)
            if(ar.model(m).condition(c).qdfxdx_nonzero(j,i))
                ar.model(m).condition(c).dfxdx{j,i} = sprintf('dfxdx[%i]', j + (i-1)*length(ar.model(m).xs));
            else
                ar.model(m).condition(c).dfxdx{j,i} = '0';
            end
        end
    end
    fprintf('dfxdx=%i, ', sum(ar.model(m).condition(c).qdfxdx_nonzero(:)));
end

% sx sensitivities
if(ar.config.useSensis)
	% su
    ar.model(m).condition(c).su = cell(length(ar.model(m).us), length(ar.model(m).condition(c).ps));
    for j=1:length(ar.model(m).us)
        for i=1:length(ar.model(m).condition(c).ps)
            if(ar.model(m).condition(c).qfu_nonzero(j))
                ar.model(m).condition(c).su{j,i} = sprintf('su[%i]', j + (i-1)*length(ar.model(m).us));
            else
                ar.model(m).condition(c).su{j,i} = '0';
            end
        end
    end
    ar.model(m).condition(c).sym.su = sym(ar.model(m).condition(c).su);
    fprintf('su=%i, ', length(ar.model(m).condition(c).ps)*sum(ar.model(m).condition(c).qfu_nonzero(:)));
    
    ar.model(m).condition(c).sym.dfudp = [];
    if(~isempty(ar.model(m).condition(c).sym.ps))
        ar.model(m).condition(c).sym.dfudp = ...
            jacobian(ar.model(m).condition(c).sym.fu, ar.model(m).condition(c).sym.ps);
    end
    
	% sx
    ar.model(m).condition(c).sx = cell(length(ar.model(m).xs), length(ar.model(m).condition(c).ps));
    for j=1:length(ar.model(m).xs)
        for i=1:length(ar.model(m).condition(c).ps)
            ar.model(m).condition(c).sx{j,i} = sprintf('sx[%i]', j);
        end
    end
	ar.model(m).condition(c).sym.sx = sym(ar.model(m).condition(c).sx);
    
    ar.model(m).condition(c).sym.fsv = ar.model(m).condition(c).sym.dvdx * ar.model(m).condition(c).sym.sx + ...
        ar.model(m).condition(c).sym.dvdu * ar.model(m).condition(c).sym.su + ar.model(m).condition(c).sym.dvdp;
    
	% sv
    ar.model(m).condition(c).qfsv_nonzero = logical(ar.model(m).condition(c).sym.fsv ~= 0);
    ar.model(m).condition(c).sv = cell(length(ar.model(m).vs), length(ar.model(m).condition(c).ps));
    for j=1:length(ar.model(m).vs)
        for i=1:length(ar.model(m).condition(c).ps)
            if(ar.model(m).condition(c).qfsv_nonzero(j,i))
                ar.model(m).condition(c).sv{j,i} = sprintf('sv[%i]', j);
            else
                ar.model(m).condition(c).sv{j,i} = '0';
            end
        end
    end
    ar.model(m).condition(c).sym.sv = sym(ar.model(m).condition(c).sv);
    fprintf('sv=%i, ', sum(ar.model(m).condition(c).qfsv_nonzero(:)));
    
    ar.model(m).condition(c).sym.fsx = (ar.model(m).N .* ar.model(m).condition(c).sym.C) * ar.model(m).condition(c).sym.sv;
    fprintf('sx=%i... ', numel(ar.model(m).condition(c).sym.fsx));
    
    % sx initials
    ar.model(m).condition(c).sym.fsx0 = jacobian(ar.model(m).condition(c).sym.fpx0, ar.model(m).condition(c).sym.ps);
    
    % steady state sensitivities
    ar.model(m).condition(c).sym.dfxdp = (ar.model(m).N .* ar.model(m).condition(c).sym.C) * (ar.model(m).condition(c).sym.dvdp + ...
        ar.model(m).condition(c).sym.dvdx*ar.model(m).condition(c).sym.fsx0 + ...
        ar.model(m).condition(c).sym.dvdu * ar.model(m).condition(c).sym.dfudp);
end

fprintf('done\n');



% OBS
function arParseOBS(m, d, doskip)
global ar

c = ar.model(m).data(d).cLink;
fprintf('parsing data m%i d%i, %s (%s -> %s)...', m, d, ar.model(m).data(d).name, ar.model(m).data(d).checkstr, ...
    ar.model(m).condition(c).checkstr);

% hard code conditions
ar.model(m).data(d).sym.p = sym(ar.model(m).data(d).p);
ar.model(m).data(d).sym.fp = sym(ar.model(m).data(d).fp);
ar.model(m).data(d).sym.fy = sym(ar.model(m).data(d).fy);
ar.model(m).data(d).sym.fy = mysubs(ar.model(m).data(d).sym.fy, ar.model(m).data(d).sym.p, ar.model(m).data(d).sym.fp);
ar.model(m).data(d).sym.fystd = sym(ar.model(m).data(d).fystd);
ar.model(m).data(d).sym.fystd = mysubs(ar.model(m).data(d).sym.fystd, ar.model(m).data(d).sym.p, ar.model(m).data(d).sym.fp);

ar.model(m).data(d).sym.fu = sym(ar.model(m).condition(c).fu);
ar.model(m).data(d).sym.fu = mysubs(ar.model(m).data(d).sym.fu, ar.model(m).data(d).sym.p, ar.model(m).data(d).sym.fp);
ar.model(m).data(d).qfu_nonzero = logical(ar.model(m).data(d).sym.fu ~= 0);

% predictor
ar.model(m).data(d).sym.fu = mysubs(ar.model(m).data(d).sym.fu, sym(ar.model(m).t), sym('t'));
ar.model(m).data(d).sym.fy = mysubs(ar.model(m).data(d).sym.fy, sym(ar.model(m).t), sym('t'));
ar.model(m).data(d).sym.fystd = mysubs(ar.model(m).data(d).sym.fystd, sym(ar.model(m).t), sym('t'));

% remaining parameters
varlist = cellfun(@symvar, ar.model(m).data(d).fp, 'UniformOutput', false);
ar.model(m).data(d).pold = ar.model(m).data(d).p;
ar.model(m).data(d).p = union(vertcat(varlist{:}), [])';

if(doskip)
    fprintf('skipping\n');
    return;
end

% make short strings
for j=1:length(ar.model(m).data(d).p)
    ar.model(m).data(d).ps{j} = sprintf('p[%i]',j);
end
ar.model(m).data(d).ys = {};
for j=1:length(ar.model(m).data(d).y)
    ar.model(m).data(d).ys{j} = sprintf('y[%i]',j);
end
fprintf('y=%i, p=%i, ', length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));

% make syms
ar.model(m).data(d).sym.p = sym(ar.model(m).data(d).p);
ar.model(m).data(d).sym.ps = sym(ar.model(m).data(d).ps);
ar.model(m).data(d).sym.y = sym(ar.model(m).data(d).y);
ar.model(m).data(d).sym.ys = sym(ar.model(m).data(d).ys);

% substitute
ar.model(m).data(d).sym.fy = mysubs(ar.model(m).data(d).sym.fy, ...
    ar.model(m).sym.x, ar.model(m).sym.xs);
ar.model(m).data(d).sym.fy = mysubs(ar.model(m).data(d).sym.fy, ...
    ar.model(m).sym.u, ar.model(m).sym.us);
ar.model(m).data(d).sym.fy = mysubs(ar.model(m).data(d).sym.fy, ...
    ar.model(m).data(d).sym.p, ar.model(m).data(d).sym.ps);

ar.model(m).data(d).sym.fystd = mysubs(ar.model(m).data(d).sym.fystd, ...
    ar.model(m).sym.x, ar.model(m).sym.xs);
ar.model(m).data(d).sym.fystd = mysubs(ar.model(m).data(d).sym.fystd, ...
    ar.model(m).sym.u, ar.model(m).sym.us);
ar.model(m).data(d).sym.fystd = mysubs(ar.model(m).data(d).sym.fystd, ...
    ar.model(m).data(d).sym.y, ar.model(m).data(d).sym.ys);
ar.model(m).data(d).sym.fystd = mysubs(ar.model(m).data(d).sym.fystd, ...
    ar.model(m).data(d).sym.p, ar.model(m).data(d).sym.ps);

% remove zero inputs
% ar.model(m).data(d).sym.fy = mysubs(ar.model(m).data(d).sym.fy, ar.model(m).sym.us(~ar.model(m).condition(c).qfu_nonzero), ...
%     sym(zeros(1,sum(~ar.model(m).condition(c).qfu_nonzero))));
% ar.model(m).data(d).sym.fystd = mysubs(ar.model(m).data(d).sym.fystd, ar.model(m).sym.us(~ar.model(m).condition(c).qfu_nonzero), ...
%     sym(zeros(1,sum(~ar.model(m).condition(c).qfu_nonzero))));

% derivatives fy
if(~isempty(ar.model(m).data(d).sym.fy))
    if(~isempty(ar.model(m).sym.us))
        ar.model(m).data(d).sym.dfydu = jacobian(ar.model(m).data(d).sym.fy, ar.model(m).sym.us);
    else
        ar.model(m).data(d).sym.dfydu = sym(ones(length(ar.model(m).data(d).y), 0));
    end
    if(~isempty(ar.model(m).x))
        ar.model(m).data(d).sym.dfydx = jacobian(ar.model(m).data(d).sym.fy, ar.model(m).sym.xs);
    else
        ar.model(m).data(d).sym.dfydx = [];
    end
	ar.model(m).data(d).sym.dfydp = jacobian(ar.model(m).data(d).sym.fy, ar.model(m).data(d).sym.ps);
else
	ar.model(m).data(d).sym.dfydu = [];
	ar.model(m).data(d).sym.dfydx = [];
	ar.model(m).data(d).sym.dfydp = [];
end

% what is measured ?
ar.model(m).data(d).qu_measured = sum(logical(ar.model(m).data(d).sym.dfydu~=0),1)>0;
ar.model(m).data(d).qx_measured = sum(logical(ar.model(m).data(d).sym.dfydx~=0),1)>0;

% derivatives fystd
if(~isempty(ar.model(m).data(d).sym.fystd))
    if(~isempty(ar.model(m).sym.us))
        ar.model(m).data(d).sym.dfystddu = jacobian(ar.model(m).data(d).sym.fystd, ar.model(m).sym.us);
    else
        ar.model(m).data(d).sym.dfystddu = sym(ones(length(ar.model(m).data(d).y), 0));
    end
    if(~isempty(ar.model(m).x))
        ar.model(m).data(d).sym.dfystddx = jacobian(ar.model(m).data(d).sym.fystd, ar.model(m).sym.xs);
    else
        ar.model(m).data(d).sym.dfystddx = [];
    end
    ar.model(m).data(d).sym.dfystddp = jacobian(ar.model(m).data(d).sym.fystd, ar.model(m).data(d).sym.ps);
    ar.model(m).data(d).sym.dfystddy = jacobian(ar.model(m).data(d).sym.fystd, ar.model(m).data(d).sym.ys);
else
	ar.model(m).data(d).sym.dfystddp = [];
	ar.model(m).data(d).sym.dfystddy = [];
    ar.model(m).data(d).sym.dfystddx = [];
end

% observed directly and exclusively
ar.model(m).data(d).dfydxnon0 = logical(ar.model(m).data(d).sym.dfydx ~= 0);

if(ar.config.useSensis)
    % sx sensitivities
    ar.model(m).data(d).sx = {};
    for j=1:length(ar.model(m).xs)
        for i=1:length(ar.model(m).condition(c).p)
            ar.model(m).data(d).sx{j,i} = sprintf('sx[%i]', j + (i-1)*length(ar.model(m).xs));
        end
    end
    ar.model(m).data(d).sym.sx = sym(ar.model(m).data(d).sx);
    
    % sy sensitivities
    ar.model(m).data(d).sy = {};
    for j=1:length(ar.model(m).data(d).sym.fy)
        for i=1:length(ar.model(m).data(d).sym.ps)
            ar.model(m).data(d).sy{j,i} = sprintf('sy[%i]', j + (i-1)*length(ar.model(m).data(d).sym.fy));
        end
    end
	ar.model(m).data(d).sym.sy = sym(ar.model(m).data(d).sy);
    
    % calculate sy
    if(~isempty(ar.model(m).data(d).sym.sy))
        ar.model(m).data(d).sym.fsy = ar.model(m).data(d).sym.dfydp;
        if(~isempty(ar.model(m).condition(c).p))
            qdynpara = ismember(ar.model(m).data(d).p, ar.model(m).condition(c).p);
        else
            qdynpara = false(size(ar.model(m).data(d).p));
        end
        
        if(~isempty(ar.model(m).condition(c).p))
            tmpfsx = ar.model(m).data(d).sym.dfydx * ...
                ar.model(m).data(d).sym.sx;
            if(~isfield(ar.model(m).condition(c), 'sym') || ~isfield(ar.model(m).condition(c).sym, 'su'))
                % su
                ar.model(m).condition(c).su = cell(length(ar.model(m).us), length(ar.model(m).condition(c).p));
                for j=1:length(ar.model(m).us)
                    for i=1:length(ar.model(m).condition(c).p)
                        if(ar.model(m).data(d).qfu_nonzero(j))
                            ar.model(m).condition(c).su{j,i} = sprintf('su[%i]', j + (i-1)*length(ar.model(m).us));
                        else
                            ar.model(m).condition(c).su{j,i} = '0';
                        end
                    end
                end
                ar.model(m).condition(c).sym.su = sym(ar.model(m).condition(c).su);
                fprintf('su=%i, ', length(ar.model(m).condition(c).p)*sum(ar.model(m).data(d).qfu_nonzero(:)));
            end
            tmpfsu = ar.model(m).data(d).sym.dfydu * ...
                ar.model(m).condition(c).sym.su;
            if(~isempty(ar.model(m).x))
                ar.model(m).data(d).sym.fsy(:,qdynpara) = ar.model(m).data(d).sym.fsy(:,qdynpara) + tmpfsx + tmpfsu;
            else
                ar.model(m).data(d).sym.fsy(:,qdynpara) = ar.model(m).data(d).sym.fsy(:,qdynpara) + tmpfsu;
            end
        end
    else
        ar.model(m).data(d).sym.fsy = [];
    end
    fprintf('sy=%i, ', sum(logical(ar.model(m).data(d).sym.fsy(:)~=0)));
    
    % calculate systd
    if(~isempty(ar.model(m).data(d).sym.sy))
        ar.model(m).data(d).sym.fsystd = ar.model(m).data(d).sym.dfystddp + ...
            ar.model(m).data(d).sym.dfystddy * ar.model(m).data(d).sym.sy;
        if(~isempty(ar.model(m).condition(c).p))
            qdynpara = ismember(ar.model(m).data(d).p, ar.model(m).condition(c).p);
        else
            qdynpara = false(size(ar.model(m).data(d).p));
        end
        
        if(~isempty(ar.model(m).condition(c).p))
            tmpfsx = ar.model(m).data(d).sym.dfystddx * ...
                ar.model(m).data(d).sym.sx;
            if(~isfield(ar.model(m).condition(c), 'sym') || ~isfield(ar.model(m).condition(c).sym, 'su'))
                % su
                ar.model(m).condition(c).su = cell(length(ar.model(m).us), length(ar.model(m).condition(c).p));
                for j=1:length(ar.model(m).us)
                    for i=1:length(ar.model(m).condition(c).p)
                        if(ar.model(m).data(d).qfu_nonzero(j))
                            ar.model(m).condition(c).su{j,i} = sprintf('su[%i]', j + (i-1)*length(ar.model(m).us));
                        else
                            ar.model(m).condition(c).su{j,i} = '0';
                        end
                    end
                end
                ar.model(m).condition(c).sym.su = sym(ar.model(m).condition(c).su);
                fprintf('su=%i, ', length(ar.model(m).condition(c).p)*sum(ar.model(m).data(d).qfu_nonzero(:)));
            end
            tmpfsu = ar.model(m).data(d).sym.dfystddu * ...
                ar.model(m).condition(c).sym.su;
            if(~isempty(ar.model(m).x))
                ar.model(m).data(d).sym.fsystd(:,qdynpara) = ar.model(m).data(d).sym.fsystd(:,qdynpara) + tmpfsx + tmpfsu;
            else
                ar.model(m).data(d).sym.fsystd(:,qdynpara) = ar.model(m).data(d).sym.fsystd(:,qdynpara) + tmpfsu;
            end
        end
    else
        ar.model(m).data(d).sym.fsystd = [];
    end
    fprintf('systd=%i, ', sum(logical(ar.model(m).data(d).sym.fsystd(:)~=0)));
end

fprintf('done\n');



% better subs
function out = mysubs(in, old, new, flag)
if(~exist('flag','var'))
    flag = 0;
end
if(~isnumeric(in) && ~isempty(old) && ~isempty(findsym(in)))
    out = subs(in, old(:), new(:), flag);
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
