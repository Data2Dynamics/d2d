function arSubplotStyle(g, labelfontsize, labelfonttype)

if(~exist('labelfontsize','var'))
    labelfontsize = 12;
end
if(~exist('labelfonttype','var'))
    labelfonttype = 'Arial';
end

set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);