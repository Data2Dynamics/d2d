function [sols, fits] = getMultipleResults(int)
% After optimization or post processing: return a set of 
% optimization solutions if the run has been finished. Returns
% a single intermediate solution if the run has not finished yet or an
% empty array if there is no intermediate solution yet.

if (isFinished(int)) 
    sArr = int.resultArr;
else
    sArr(1) = int.mp.getIntermediateResult();
end

fits=zeros(size(sArr,1),1);

for i=1:size(sArr,1)
    if (isempty(int.range)) % binary case
        sols(i,:)=convertUnsignedJE(int, sArr(i));
    else
        sols(i,:)=sArr(i,:);
    end;
    %disp(sols(i,:));
    if (isempty(int.args))
        fits(i) = feval(int.f, sols(i,:));
    else
        fits(i) = feval(int.f, sols(i,:), int.args);
    end
end