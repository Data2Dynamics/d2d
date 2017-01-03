% Print results of PLE
%
% plePrint(fid)
%
% fid:     file ID     [1=Standard Output to Console]

function plePrint(fid)

global pleGlobals;

if(isempty(pleGlobals))
    error('perform ple before usage');
end
if(isempty(pleGlobals.ps))
    return
end

if(~exist('fid', 'var'))
    fid = 1;
end

if(isfield(pleGlobals, 'p_true'))
    spacer = '-----------------------------------------------------------------------------------------------------------------------------------\n';
else
    spacer = '---------------------------------------------------------------------------------------------------------------\n';
end

parnamelength = length('Parameter:');
for j=1:length(pleGlobals.p_labels)
	if(length(pleGlobals.p_labels{j}) > parnamelength)
		parnamelength = length(pleGlobals.p_labels{j});
	end
end
partemplate = sprintf('%%%is ', parnamelength);
numlength = 10;
numtemplate = sprintf('%%%is ', numlength);
coverlength = 2;
covertemplate = sprintf('%%%is ', coverlength);

%% point-wise
fprintf(fid, '\nPLE %2i%% point-wise CI:\n', (1-pleGlobals.alpha_level)*100);
fprintf(fid, spacer);
fprintf(fid, '    # ');
fprintf(fid, partemplate, makestr('Parameter:', parnamelength));
fprintf(fid, numtemplate, makestr('Value:', numlength));
if(isfield(pleGlobals, 'p_true'))
    fprintf(fid, numtemplate, makestr('True:', numlength));
end
fprintf(fid, numtemplate, makestr('Min:', numlength));
fprintf(fid, numtemplate, makestr('Max:', numlength));
fprintf(fid, '|');
fprintf(fid, numtemplate, makestr('IDflag:', numlength));
fprintf(fid, numtemplate, makestr('LowerPL:', numlength));
fprintf(fid, numtemplate, makestr('UpperPL:', numlength));
fprintf(fid, numtemplate, makestr('Rel.PL:', numlength));
if(isfield(pleGlobals, 'p_true'))
    fprintf(fid, '|');
    fprintf(fid, covertemplate, makestr('?:', numlength));
end
fprintf(fid, '\n');

fprintf(fid, spacer);
for j=1:length(pleGlobals.p)
    idflag = '';
    if(pleGlobals.IDstatus_point(j)>1)
        idflag = pleGlobals.IDlabel{pleGlobals.IDstatus_point(j)};
    end
    
    fprintf(fid, ' %4i ', j);
    fprintf(fid, partemplate, makestr(pleGlobals.p_labels{j}, parnamelength));
    fprintf(fid, numtemplate, makestr(pleGlobals.p(j), numlength));
    if(isfield(pleGlobals, 'p_true'))
        fprintf(fid, numtemplate, makestr(pleGlobals.p_true(j), numlength));
    end
    if(pleGlobals.q_fit(j))
        fprintf(fid, numtemplate, makestr(pleGlobals.lb(j), numlength));
        fprintf(fid, numtemplate, makestr(pleGlobals.ub(j), numlength));
        fprintf(fid, '|');
        fprintf(fid, numtemplate, makestr(idflag, numlength));
        fprintf(fid, numtemplate, makestr(pleGlobals.conf_lb_point(j), numlength));
        fprintf(fid, numtemplate, makestr(pleGlobals.conf_ub_point(j), numlength));
        fprintf(fid, numtemplate, makestr(pleGlobals.conf_rel_point(j), numlength));
        if(isfield(pleGlobals, 'p_true'))
            fprintf(fid, '|');
            coverage_flag = ' ';
            if(pleGlobals.cover_point(j))
                coverage_flag = '*';
            end
            fprintf(fid, covertemplate, makestr(coverage_flag, numlength));
        end
    else
        fprintf(fid, '*** fixed');
    end
    fprintf(fid, '\n');
end
fprintf(fid, spacer);

%% simultaneous
if(pleGlobals.plot_simu)
    fprintf(fid, '\nPLE %2i%% simultaneous CI:\n', (1-pleGlobals.alpha_level)*100);
    fprintf(fid, spacer);
    fprintf(fid, '    # ');
    fprintf(fid, partemplate, makestr('Parameter:', parnamelength));
    fprintf(fid, numtemplate, makestr('Value:', numlength));
    if(isfield(pleGlobals, 'p_true'))
        fprintf(fid, numtemplate, makestr('True:', numlength));
    end
    fprintf(fid, numtemplate, makestr('Min:', numlength));
    fprintf(fid, numtemplate, makestr('Max:', numlength));
    fprintf(fid, '|');
    fprintf(fid, numtemplate, makestr('IDflag:', numlength));
    fprintf(fid, numtemplate, makestr('LowerPL:', numlength));
    fprintf(fid, numtemplate, makestr('UpperPL:', numlength));
    fprintf(fid, numtemplate, makestr('Rel.PL:', numlength));
    if(isfield(pleGlobals, 'p_true'))
        fprintf(fid, covertemplate, makestr('?:', numlength));
    end
    fprintf(fid, '\n');
    
    fprintf(fid, spacer);
    for j=1:length(pleGlobals.p)
        idflag = '';
        if(pleGlobals.IDstatus(j)>1)
            idflag = pleGlobals.IDlabel{pleGlobals.IDstatus(j)};
        end
        
        fprintf(fid, ' %4i ', j);
        fprintf(fid, partemplate, makestr(pleGlobals.p_labels{j}, parnamelength));
        fprintf(fid, numtemplate, makestr(pleGlobals.p(j), numlength));
        if(isfield(pleGlobals, 'p_true'))
            fprintf(fid, numtemplate, makestr(pleGlobals.p_true(j), numlength));
        end
        if(pleGlobals.q_fit(j))
            fprintf(fid, numtemplate, makestr(pleGlobals.lb(j), numlength));
            fprintf(fid, numtemplate, makestr(pleGlobals.ub(j), numlength));
            fprintf(fid, '|');
            fprintf(fid, numtemplate, makestr(idflag, numlength));
            fprintf(fid, numtemplate, makestr(pleGlobals.conf_lb(j), numlength));
            fprintf(fid, numtemplate, makestr(pleGlobals.conf_ub(j), numlength));
            fprintf(fid, numtemplate, makestr(pleGlobals.conf_rel(j), numlength));
            if(isfield(pleGlobals, 'p_true'))
                coverage_flag = ' ';
                if(pleGlobals.p_true(j) <= pleGlobals.conf_ub(j) && ...
                        pleGlobals.p_true(j) >= pleGlobals.conf_lb(j))
                    coverage_flag = '*';
                end
                fprintf(fid, covertemplate, makestr(coverage_flag, numlength));
            end
        else
            fprintf(fid, '*** fixed');
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, spacer);
end



function strout = makestr(sthin, width)
if(isnumeric(sthin))
    strout = sprintf('%+.4g', sthin);
elseif(ischar(sthin))
    strout = sprintf('%s', sthin);
    if(length(strout)>width)
        strout = strout(1:width);
    end
end



