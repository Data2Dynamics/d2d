function arMergeFits
% useful if several workspaces e.g. from arFitClusterLHS with randomseed were
% saved
% select workspaces to be compared and check overall performance
% e.g. by arPlotChi2


global ar

if(isempty(ar))
    error('please initialize by arInit')
end

filenames = fileChooserMulti('./Results', true);

if(~iscell(filenames))
    filelist = fileList('./Results');
    filenames = filelist(filenames);
end

jcount = 0;
kcount = 1;
excount = 1;
arWaitbar(0);


for j=1:length(filenames)
    fname = ['./Results/' filenames{j} '/workspace.mat'];
    load(fname);
    % tmpple.ar could be used instead
    jcount = jcount + 1;
    chi2send = kcount - 1 + length(ar.chi2s);
    exitflagend = excount - 1 + length(ar.exitflag);
    chi2s(kcount:chi2send) = ar.chi2s;
    chi2sconstr(kcount:chi2send) = ar.chi2sconstr;
    chi2s_start(kcount:chi2send) = ar.chi2s_start;
    chi2sconstr_start(kcount:chi2send) = ar.chi2sconstr_start;
    optim_crit(kcount:chi2send) = ar.optim_crit;
    ps(kcount:chi2send,:) = ar.ps;
    ps_start(kcount:chi2send,:) = ar.ps_start;
    exitflag(excount:exitflagend) = ar.exitflag;
    kcount = kcount + length(ar.chi2s);
    excount = excount + length(ar.exitflag);
    ar.chi2s = chi2s;
    ar.chi2s_start = chi2s_start;
    ar.chi2sconstr_start = chi2sconstr_start;
    ar.ps = ps;
    ar.ps_start = ps_start;
    ar.optim_crit = optim_crit;
    ar.chi2sconstr = chi2sconstr;
    ar.exitflag = exitflag;
end

arWaitbar(-1);

end
