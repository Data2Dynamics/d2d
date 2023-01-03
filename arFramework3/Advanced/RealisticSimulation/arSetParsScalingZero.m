
function arSetParsScalingZero(range)

global ar

if ~exist('range','var') || isempty(range)
    range = 2;
end
arLink
try
arSimu(false,true)
end
for i=1:length(ar.p)
    if strncmp(ar.pLabel{i},'offset_',6)
        ar.p(i) = min(min(log10(ar.model.condition.xFineSimu(~isinf(log10(ar.model.condition.xFineSimu))))));
        ar.qLog10(i) = 1;
        ar.lb(i) = ar.p(i)-range;
        ar.ub(i) = ar.p(i)+range;
    end
    if strncmp(ar.pLabel{i},'scale_',6)
        ar.lb(i) = -range;
        ar.p(i) = 0;
        ar.ub(i) = range;
    end
end

% arLink
% arSimu(true,true)   % Integrate to get corret initial value
% for i=1:length(ar.p)
%     if (ar.qInitial(i) && ar.p(i)==0 && ~ar.qLog10(i))
%         for j=1:length(ar.model.data.y)
%             if strcmp(ar.pLabel{i}(6:end),ar.model.data.y{j}(1:end-4))
%                 ar.p(i) = log10(ar.model(1).data.yFineSimu(2,j));
%                 ar.lb(i) = ar.p(i)-range;
%                 ar.ub(i) = ar.p(i)+range;
%                 ar.qLog10(i) = 1;
%             end
%         end
%     end
% end
