function [T,yExp] = artExpToVector(tT,y)

global ar

% Replace too small values
y(y<1e-8 & y>0)=1e-8;
y(y>-1e-6 & y<0)=1e-8;
if length(y(y<1e-7))<size(y,2)
    y(y<1e-7)=1e-7;
end
if length(y(y<1e-6))<size(y,2)
    y(y<1e-6)=1e-6;
end
if length(y(y<1e-5))<size(y,2)
    y(y<1e-5)=1e-5;
end
if length(y(y<1e-4))<size(y,2)
    y(y<1e-4)=1e-4;
end
if length(y(y<1e-3))<size(y,2)
    y(y<1e-3)=1e-3;
end

% Order data all in one matrix
T = unique(sort(tT(~isnan(tT(:)))));
yExp = nan(size(T,1),size(tT,2));
for i=1:size(y,1)
    for j=1:size(y,2)
        for k=1:size(T,1)
            if tT(i,j)==T(k)
                yExp(k,j) = y(i,j);
            end
        end
    end
end
while all(isnan(yExp(end,:)))
    yExp(end,:) = [];
    T(end) = [];
end

end