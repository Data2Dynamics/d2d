% Print results of PLE
%
% plePrint(fid)
%
% fid:     file ID     [1=Standard Output to Console]

function plePrint(fid)

global ar

if(~isfield(ar,'ple') || isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
end

if(~exist('fid', 'var'))
    fid = 1;
end

if(isfield(ar.ple, 'p_true'))
    spacer = '-----------------------------------------------------------------------------------------------------------------------------------\n';
else
    spacer = '---------------------------------------------------------------------------------------------------------------\n';
end

parnamelength = length('Parameter:');
for j=1:length(ar.ple.p_labels)
	if(length(ar.ple.p_labels{j}) > parnamelength)
		parnamelength = length(ar.ple.p_labels{j});
	end
end
partemplate = sprintf('%%%is ', parnamelength);
numlength = 10;
numtemplate = sprintf('%%%is ', numlength);
coverlength = 2;
covertemplate = sprintf('%%%is ', coverlength);

%% point-wise
fprintf(fid, '\nPLE %2i%% point-wise CI:\n', (1-ar.ple.alpha_level)*100);
fprintf(fid, spacer);
fprintf(fid, '    # ');
fprintf(fid, partemplate, makestr('Parameter:', parnamelength));
fprintf(fid, numtemplate, makestr('Value:', numlength));
if(isfield(ar.ple, 'p_true'))
    fprintf(fid, numtemplate, makestr('True:', numlength));
end
fprintf(fid, numtemplate, makestr('Min:', numlength));
fprintf(fid, numtemplate, makestr('Max:', numlength));
fprintf(fid, '|');
fprintf(fid, numtemplate, makestr('IDflag:', numlength));
fprintf(fid, numtemplate, makestr('LowerPL:', numlength));
fprintf(fid, numtemplate, makestr('UpperPL:', numlength));
if(isfield(ar.ple, 'p_true'))
    fprintf(fid, '|');
    fprintf(fid, covertemplate, makestr('?:', numlength));
end
fprintf(fid, '\n');

fprintf(fid, spacer);
for j=1:length(ar.ple.p)
    idflag = '';
    if(ar.ple.IDstatus_point(j)>1)
        idflag = ar.ple.IDlabel{ar.ple.IDstatus_point(j)};
    end
    
    fprintf(fid, ' %4i ', j);
    fprintf(fid, partemplate, makestr(ar.ple.p_labels{j}, parnamelength));
    fprintf(fid, numtemplate, makestr(ar.ple.p(j), numlength));
    if(isfield(ar.ple, 'p_true'))
        fprintf(fid, numtemplate, makestr(ar.ple.p_true(j), numlength));
    end
    if(ar.qFit(j))
        fprintf(fid, numtemplate, makestr(ar.lb(j), numlength));
        fprintf(fid, numtemplate, makestr(ar.ub(j), numlength));
        fprintf(fid, '|');
        fprintf(fid, numtemplate, makestr(idflag, numlength));
        fprintf(fid, numtemplate, makestr(ar.ple.conf_lb_point(j), numlength));
        fprintf(fid, numtemplate, makestr(ar.ple.conf_ub_point(j), numlength));
        if(isfield(ar.ple, 'p_true'))
            fprintf(fid, '|');
            coverage_flag = ' ';
            if(ar.ple.cover_point(j))
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
if(ar.ple.plot_simu)
    fprintf(fid, '\nPLE %2i%% simultaneous CI:\n', (1-ar.ple.alpha_level)*100);
    fprintf(fid, spacer);
    fprintf(fid, '    # ');
    fprintf(fid, partemplate, makestr('Parameter:', parnamelength));
    fprintf(fid, numtemplate, makestr('Value:', numlength));
    if(isfield(ar.ple, 'p_true'))
        fprintf(fid, numtemplate, makestr('True:', numlength));
    end
    fprintf(fid, numtemplate, makestr('Min:', numlength));
    fprintf(fid, numtemplate, makestr('Max:', numlength));
    fprintf(fid, '|');
    fprintf(fid, numtemplate, makestr('IDflag:', numlength));
    fprintf(fid, numtemplate, makestr('LowerPL:', numlength));
    fprintf(fid, numtemplate, makestr('UpperPL:', numlength));
    if(isfield(ar.ple, 'p_true'))
        fprintf(fid, covertemplate, makestr('?:', numlength));
    end
    fprintf(fid, '\n');
    
    fprintf(fid, spacer);
    for j=1:length(ar.ple.p)
        idflag = '';
        if(ar.ple.IDstatus(j)>1)
            idflag = ar.ple.IDlabel{ar.ple.IDstatus(j)};
        end
        
        fprintf(fid, ' %4i ', j);
        fprintf(fid, partemplate, makestr(ar.ple.p_labels{j}, parnamelength));
        fprintf(fid, numtemplate, makestr(ar.ple.p(j), numlength));
        if(isfield(ar.ple, 'p_true'))
            fprintf(fid, numtemplate, makestr(ar.ple.p_true(j), numlength));
        end
        if(ar.qFit(j))
            fprintf(fid, numtemplate, makestr(ar.lb(j), numlength));
            fprintf(fid, numtemplate, makestr(ar.ub(j), numlength));
            fprintf(fid, '|');
            fprintf(fid, numtemplate, makestr(idflag, numlength));
            fprintf(fid, numtemplate, makestr(ar.ple.conf_lb(j), numlength));
            fprintf(fid, numtemplate, makestr(ar.ple.conf_ub(j), numlength));
            if(isfield(ar.ple, 'p_true'))
                coverage_flag = ' ';
                if(ar.ple.p_true(j) <= ar.ple.conf_ub(j) && ...
                        ar.ple.p_true(j) >= ar.ple.conf_lb(j))
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



