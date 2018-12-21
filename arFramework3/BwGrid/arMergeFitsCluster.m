function arMergeFitsCluster(scriptname)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end


myfilelist= fileList('./Results', [scriptname]);
myfilelist2= fileList(['./Results/' scriptname], [scriptname]);

if(isempty(myfilelist))
    error( [scriptname ' not found in Results'])
else
    disp([ num2str(length(myfilelist)) ' results found. Merging fits...' ])
    disp(myfilelist)
end
if ~isempty(myfilelist2)
	disp(myfilelist2)
	myfilelist = [myfilelist,myfilelist2]
	disp(myfilelist)
end

jcount = 0;
kcount = 1;
excount = 1;
arWaitbar(0);


for j=1:length(myfilelist)
    fname = ['./Results/' myfilelist{j} '/workspace.mat'];
    arWaitbar(j,length(myfilelist),['loading ' fname] );
    load(fname);
    % tmpple.ar could be used instead
    jcount = jcount + 1;
    singlefit = 0;
    if isfield(ar,'chi2s') && ~singlefit
	    chi2send = kcount - 1 + length(ar.chi2s);
	    exitflagend = excount - 1 + length(ar.exitflag);
	    chi2s(kcount:chi2send) = ar.chi2s;
	    chi2sconstr(kcount:chi2send) = ar.chi2sconstr;
	    chi2s_start(kcount:chi2send) = ar.chi2s_start;
	    chi2sconstr_start(kcount:chi2send) = ar.chi2sconstr_start;
            optim_crit(kcount:chi2send) = ar.optim_crit;
	    timing(kcount:chi2send) = ar.timing;
	    iter(kcount:chi2send) = ar.iter;
	    ps(kcount:chi2send,:) = ar.ps;
	    ps_errors(kcount:chi2send,:) = ar.ps_errors;
	    ps_start(kcount:chi2send,:) = ar.ps_start;
	    exitflag(excount:exitflagend) = ar.exitflag;
	    kcount = kcount + length(ar.chi2s);
	    excount = excount + length(ar.exitflag);
	    ar.chi2s = chi2s;
	    ar.chi2s_start = chi2s_start;
	    ar.chi2sconstr_start = chi2sconstr_start;
	    ar.ps = ps;
	    ar.ps_errors = ps_errors;
	    ar.ps_start = ps_start;
	    ar.optim_crit = optim_crit;
	    ar.timing = timing;
	    ar.iter = iter;
	    ar.chi2sconstr = chi2sconstr;
	    ar.exitflag = exitflag;
     else
	    singlefit = 1;
            chi2s(kcount) = ar.chi2s;
            chi2sconstr(kcount) = ar.chi2sconstr;
            chi2s_start(kcount) = ar.chi2s_start;
            chi2sconstr_start(kcount) = ar.chi2sconstr_start;
            ps(kcount,:) = ar.p;
            optim_crit(kcount) = ar.optim_crit;
            timing(kcount) = ar.timing;
            iter(kcount) = ar.iter;
	    exitflag(kcount) = ar.exitflag;
            kcount = kcount + 1;
            ar.chi2s = chi2s;
            ar.chi2sconstr = chi2sconstr;
            ar.chi2s_start = chi2s_start;
            ar.chi2sconstr_start = chi2sconstr_start;
            ar.ps = ps;
            ar.ps_errors = ps_errors;
            ar.ps_start = ps_start;
            ar.optim_crit = optim_crit;
            ar.timing = timing;
            ar.iter = iter;
	    ar.exitflag = exitflag;
     end
end

arWaitbar(-1);
fprintf('merged\n')

name = ['ClusterMerge_' scriptname '_' num2str(kcount-1) 'fits'];
arSave(name)
[~,ws]=fileparts(ar.config.savepath);
movefile(['Results/' ws],['Results/' name])
fprintf(['Saved to ' name '\n']);

[~,imin] = min(ar.chi2s + ar.chi2sconstr);
ar.p = ar.ps(imin,:);
arCalcMerit(true,[]);
fprintf('Found Minimum.\n')
arSimu
fprintf('Calculated Integrals.\n')

name = [name '_Min'];
arSave(name)
[~,ws]=fileparts(ar.config.savepath);
movefile(['Results/' ws],['Results/' name])
fprintf(['Minimum saved to ' name '\n']);
