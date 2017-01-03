function xytickset(labels, doY, rotation)

if(~exist('doY','var'))
    doY = false;
end
if(~exist('rotation','var'))
    rotation = 90;
end

if(doY)
    set(gca, 'YTick', 1:length(labels));
    set(gca, 'YTickLabel', labels);
    set(gca, 'YLim', [0 length(labels)+1]);
else
    set(gca, 'XTick', 1:length(labels));
    set(gca, 'XTickLabel', labels);
    set(gca, 'XTickLabelRotation', rotation);
    set(gca, 'XLim', [0 length(labels)+1]);
end