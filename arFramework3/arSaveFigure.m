function savePath = arSaveFigure(h, name, subdir)

savePath = [arSave subdir];

if(~exist(savePath, 'dir'))
	mkdir(savePath)
end

if(length(name)>60)
    name = name([1:29 (end-29):end]);
end

savePath = arPathConvert([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc', savePath);
% print('-dpng', savePath);
if(ispc)
    print('-dpdf', savePath);
elseif(ismac)
    system(['/usr/local/bin/ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
else
    system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
end


% if(exist('plot2svg','file')==2)
%     myAxisPermutation(h)
%     plot2svg([savePath '.svg'], h);
% end
% if(exist('matlab2tikz','file')==2)
%     set(h,'Units','in')
%     myAxisPermutation(h)
%     matlab2tikz([savePath '.tex'],'figurehandle',h,'showInfo', false, 'width','0.9\textwidth')
% end



% function myAxisPermutation(h)
% % Permute lines and patches
% haxes = findobj(h, 'type', 'axes');
% for iaxis = 1:length(haxes)
%     hlines = findobj(haxes(iaxis), 'type', 'line');
%     axischild = get(haxes(iaxis),'Children');
%     lineind = ismember(axischild,hlines);
%     hpatches = findobj(haxes(iaxis), 'type', 'patch');
%     patchind = find(ismember(axischild,hpatches));
%     patchind_sorted = [];
%     for ipatch = patchind'
%         if(isnumeric(get(axischild(ipatch),'FaceColor')))
%             % Patch with FaceColor to the background
%             patchind_sorted = [patchind_sorted ipatch];
%         else
%             % Patch without FaceColor to the front
%             patchind_sorted = [ipatch patchind_sorted];
%         end
%     end
%     % Find other handles (e.g. text)
%     otherhandles = ~ismember(1:length(axischild),[find(lineind)' patchind_sorted]);
%     % Lines to the front, then patches, then other handles
%     set(haxes(iaxis),'Children',[axischild(lineind); axischild(patchind_sorted); axischild(otherhandles)])
% end
