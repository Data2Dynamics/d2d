% arMergeFits
%
% arMergeFits merges fits from different workspaces which are selected via
% fileChooser. 
% 
% This is useful if several workspaces e.g. from arFitClusterLHS 
% with different randomseeds were saved select workspaces to be compared and 
% check overall performance e.g. by arPlotChi2 or arPlotFits
% 
% See also arFits, arPlotChi2, arPlotFits, arMergeFitsCore

function arMergeFits

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

filenames = fileChooserMulti('./Results', true);

if(~iscell(filenames))
    filelist = fileList('./Results');
    filenames = filelist(filenames);
end

arWaitbar(0);


for j=1:length(filenames)
    fname = ['./Results/' filenames{j} '/workspace.mat'];
    tmp = load(fname,'ar');
    if isfield(tmp,'ar')
        ar = arMergeFitsCore(ar,tmp.ar);
    else
        warning('%s does not contain a D2D-struct ar',fname);
    end
    % tmpple.ar could be used instead
    
    
end

arWaitbar(-1);

end
