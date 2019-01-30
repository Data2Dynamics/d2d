% argout = arDeleteGraphicsHandles(ar2)
% This function deletes the graphcis handles in ar which are
% since R2014b objects that are displayed when loaded in a workspace.
%
% This function prevents this.
%
% The function can also be used to open new figures when calling arPlot a
% second time, e.g. usefull for a fast comparison of differnt parameter
% sets. 
%
% arPlot; arPlot 
% -> Second arPlot plots results in the same figure.
% arPlot; arDeleteGraphicsHandles; arPlot 
% -> All figures are opened twice.

function argout = arDeleteGraphicsHandles(ar2)

global ar
if nargin == 0 % in arSave a local copy of the ar-struct is saved and the global struct remains unchanged.  
    ar2 = ar;
end


if isstruct(ar2)
    for m=1:length(ar2.model)
        for p=1:length(ar2.model(m).plot)
            if isfield(ar2.model(m).plot(p),'fighandel_x')
                if ~isempty(ar2.model(m).plot(p).fighandel_x)
                    %delete(ar.model(m).plot(p).fighandel_x);% this line would close the plot
                    ar2.model(m).plot(p).fighandel_x = [];
                end
            end
            if isfield(ar2.model(m).plot(p),'fighandel_v')
                if ~isempty(ar2.model(m).plot(p).fighandel_v)
                    %delete(ar.model(m).plot(p).fighandel_v);% this line would close the plot
                    ar2.model(m).plot(p).fighandel_v = [];
                end
            end
            if isfield(ar2.model(m).plot(p),'fighandel_y')
                if ~isempty(ar2.model(m).plot(p).fighandel_y)
                    %delete(ar.model(m).plot(p).fighandel_y);  % this line would close the plot
                    ar2.model(m).plot(p).fighandel_y = [];
                end
            end
        end
    end
    if isfield(ar2.ple,'fighandel_multi')
        if ~isempty(ar2.ple.fighandel_multi)
            ar2.ple.fighandel_multi = [];
        end
    end
end
if nargout == 1
    argout = ar2;
end