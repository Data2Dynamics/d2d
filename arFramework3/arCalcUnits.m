function arCalcUnits(m)


clc

global ar

if(~exist('m','var'))
    m = 1;
end

%% reset all entries
for jx = 1:length(ar.model(m).x)
    ar.model(m).xUnits{jx,1} = '';
end

%% species units
for jv = 1:length(ar.model(m).fv)
    % [ar.model(m).x' ar.model(m).xUnits]
    xs = union(ar.model(m).fv_source{jv}, ar.model(m).fv_target{jv});
    q = ismember(ar.model(m).x, xs);
    
    xu = ar.model(m).xUnits(q,1);
    xul = ar.model(m).xUnits(q,2);
    xuu = unique(xu);
    xuul = unique(xul);
    if(length(xuul)>1)
        fprintf(2, 'Inconsistent species units:\n');
        for jj=find(q)
            fprintf(2, '\t%s [%s]\n', ar.model(m).x{jj}, ar.model(m).xUnits{jj,2});
        end
        error('Inconsistent species units');
    end
    xuu = setdiff(xuu, {''});
    if(isempty(xuu))
        ar.model(m).xUnits(q,1) = repmat(xs(1), sum(q),1);
    elseif(length(xuu)==1)
        ar.model(m).xUnits(q,1) = repmat(xuu(1), sum(q),1);
    else
        ar.model(m).xUnits(q,1) = repmat(xuu(1), sum(q),1);
        for jxuu=2:length(xuu)
            ar.model(m).xUnits(:,1) = strrep(ar.model(m).xUnits(:,1), ...
                xuu(jxuu), xuu(1));
        end
    end
end

known_vars = ar.model(m).x;
known_units = ar.model(m).xUnits(:,1)';

%% parameter units for volumes and initial conditions
ar.model(m).pUnits = cell(size(ar.model(m).p));
for jp = 1:length(ar.model(m).p)
    % volumes
    q = ismember(ar.model(m).pc, ar.model(m).p{jp});
    if(sum(q) == 1)
        ar.model(m).pUnits(jp) = ar.model(m).c(q);
        known_vars{end+1} = ar.model(m).p{jp}; %#ok<AGROW>
        known_units(end+1) = ar.model(m).pUnits(jp); %#ok<AGROW>
        continue
    elseif(sum(q) > 1)
        error('double matches');
    end
    
    % initial conditions
    q = ismember(ar.model(m).px0, ar.model(m).p{jp});
    if(sum(q) == 1)
        ar.model(m).pUnits(jp) = ar.model(m).xUnits(q,1);
        known_vars{end+1} = ar.model(m).p{jp}; %#ok<AGROW>
        known_units(end+1) = ar.model(m).pUnits(jp); %#ok<AGROW>
        continue
    elseif(sum(q) > 1)
        error('double matches');
    end
end

%% remaining dynamic parameters
px = {};
for jp = 1:length(ar.model(m).p)
    q = ismember(ar.model(m).pv, ar.model(m).p{jp});
    if(sum(q) == 1)
        px{end+1} = ar.model(m).p{jp}; %#ok<AGROW>
    elseif(sum(q) > 1)
        error('double matches');
    end
end

%% setup flux equations for unit solver
fv = ar.model(m).fv;
fv = sym(fv);
fv = subs(fv, sym(known_vars), sym(known_units));
for jv = 1:length(fv)
    xs = union(ar.model(m).fv_source{jv}, ar.model(m).fv_target{jv});
    fv(jv) = fv(jv) * sym(ar.model(m).tUnits{1}) / sym(xs{1}) - 1;
end

punits = unit_solver(fv, px);

% %% remaining dynamic parameters
% 
% for jp = 1:length(ar.model(m).p)
%     q = ismember(ar.model(m).pv, ar.model(m).p{jp});
%     if(sum(q) == 1)
%         ar.model(m).p{jp}
%         q = ismembern(ar.model(m).pvs,ar.model(m).p{jp});
%         vs = ar.model(m).fv(q);
%         vsi = find(q);
%         xs = union(ar.model(m).fv_source{vsi(1)}, ar.model(m).fv_target{vsi(1)});
%         xsu = ar.model(m).xUnits{ismember(ar.model(m).x, xs{1}),1};
%         vs = sym(vs);
%         vs = subs(vs, sym(ar.model(m).x), sym(ar.model(m).xUnits(:,1)'))
%         vs = vs * sym(ar.model(m).tUnits{1}) / sym(xsu) - 1;
%         solve(vs(:), ar.model(m).p{jp})
%         
%         ptmp = '';
%         
%         ar.model(m).pUnits{jp} = ptmp;
%         continue
%     elseif(sum(q) > 1)
%         error('double matches');
%     end
% end


function q = ismembern(c,l)
q = logical(size(c));
for j=1:length(c)
    q(j) = sum(ismember(c{j},l))>0;
end
    
function punits = unit_solver(fv, px)

psolved = {};
punits = {};
px = sym(px);

for j=1:length(px)
    fvtmp = findasym(fv, px(j));
    
    px(j)
    fvtmp
    S = solve(fvtmp, px(j))
end


function fxout = findasym(fx, x)
fxout = sym([]);
for j=1:length(fx)
    q = sum(ismember(symvar(fx(j)),x))>0;
    if(q)
        fxout(end+1) = fx(j); %#ok<AGROW>
    end
end



