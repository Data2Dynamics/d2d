function arSubplotStyle(g, labelfontsize, labelfonttype)

if(~exist('labelfontsize','var'))
    labelfontsize = 10;
end
if(~exist('labelfonttype','var'))
    labelfonttype = 'Arial';
end

set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);