% Print results of PLE
%
% plePrint(fid)
%
% fid:     file ID     [1=Standard Output to Console]
% whichParams   indicates which parameter values should be printed
%          [] default: all
%          indices   a number or vector of indices
%          pattern   a string which is searched in p_labels as regexp
% 
% Examples:
% plePrint([],1:3)
% plePrint([],'fold')

function plePrint(fid, whichParams, doUnlog)
global ar

if ~exist('fid','var') || isempty(fid)
    fid = 1;
end
if ~exist('whichParams','var') || isempty(whichParams)
    whichParams = 1:length(ar.ple.p);
elseif ischar(whichParams)
    whichParams = find(~cellfun(@isempty,regexp(ar.ple.p_labels,whichParams)));
elseif(size(whichParams,1)>1)
        whichParams = js'; %should not be a row
end
if ~exist('doUnlog','var') || isempty(doUnlog)
    doUnlog = 0;
end

if(~isfield(ar,'ple') || isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
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
fprintf(fid, [' ',numtemplate], makestr('LowerCI:', numlength));
fprintf(fid, [numtemplate,' '], makestr('UpperCI:', numlength));
if(isfield(ar.ple, 'p_true'))
    fprintf(fid, numtemplate, makestr('True:', numlength));
end
fprintf(fid, numtemplate, makestr('Min:', numlength));
fprintf(fid, numtemplate, makestr('Max:', numlength));
fprintf(fid, '|');
fprintf(fid, numtemplate, makestr('IDflag:', numlength));
if(isfield(ar.ple, 'p_true'))
    fprintf(fid, '|');
    fprintf(fid, covertemplate, makestr('?:', numlength));
end
fprintf(fid, '\n');

fprintf(fid, spacer);
for j=whichParams
    idflag = '';
    if(ar.ple.IDstatus_point(j)>1)
        idflag = ar.ple.IDlabel{ar.ple.IDstatus_point(j)};
    end
    
    fprintf(fid, ' %4i ', j);
    fprintf(fid, partemplate, makestr(ar.ple.p_labels{j}, parnamelength));
    fprintf(fid, numtemplate, makestr(unlogTrsf(ar.ple.p(j),ar.ple.qLog10(j),doUnlog), numlength));
    if(isfield(ar.ple, 'p_true'))
        fprintf(fid, numtemplate, makestr(unlogTrsf(ar.ple.p_true(j),ar.ple.qLog10(j),doUnlog), numlength));
    end
    if(ar.qFit(j))
        fprintf(fid, ['[',numtemplate], makestr(unlogTrsf(ar.ple.conf_lb_point(j),ar.ple.qLog10(j),doUnlog), numlength));
        fprintf(fid, [numtemplate,']'], makestr(unlogTrsf(ar.ple.conf_ub_point(j),ar.ple.qLog10(j),doUnlog), numlength));
        fprintf(fid, numtemplate, makestr(unlogTrsf(ar.lb(j),ar.ple.qLog10(j),doUnlog), numlength));
        fprintf(fid, numtemplate, makestr(unlogTrsf(ar.ub(j),ar.ple.qLog10(j),doUnlog), numlength));
        fprintf(fid, '|');
        fprintf(fid, numtemplate, makestr(idflag, numlength));
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
    fprintf(fid, [' ',numtemplate], makestr('LowerCI:', numlength));
    fprintf(fid, [numtemplate,' '], makestr('UpperCI:', numlength));
    fprintf(fid, numtemplate, makestr('Min:', numlength));
    fprintf(fid, numtemplate, makestr('Max:', numlength));
    fprintf(fid, '|');
    fprintf(fid, numtemplate, makestr('IDflag:', numlength));
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
        fprintf(fid, numtemplate, makestr(unlogTrsf(ar.ple.p(j),ar.ple.qLog10(j),doUnlog), numlength));
        if(isfield(ar.ple, 'p_true'))
            fprintf(fid, numtemplate, makestr(unlogTrsf(ar.ple.p_true(j),ar.ple.qLog10(j),doUnlog), numlength));
        end
        if(ar.qFit(j))
            fprintf(fid, ['[',numtemplate], makestr(unlogTrsf(ar.ple.conf_lb(j),ar.ple.qLog10(j),doUnlog), numlength));
            fprintf(fid, [numtemplate,' '], makestr(unlogTrsf(ar.ple.conf_ub(j),ar.ple.qLog10(j),doUnlog), numlength));
            fprintf(fid, numtemplate, makestr(unlogTrsf(ar.lb(j),ar.ple.qLog10(j),doUnlog), numlength));
            fprintf(fid, numtemplate, makestr(unlogTrsf(ar.ub(j),ar.ple.qLog10(j),doUnlog), numlength));
            fprintf(fid, '|');
            fprintf(fid, numtemplate, makestr(idflag, numlength));
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

function val = unlogTrsf(val, qLog, doUnlog)
if doUnlog && qLog
    val = 10.^val;
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



