% plotApproxError(fits)
% 
%   Histogram of the approximation error defined as
% 
% Example
% arLoad
% fits = arApproximateTimeCoursesByTransientFunction2;
% plotApproxError(fits)

function approxErr = plotApproxError(fits,const_thresh)
if ~exist('const_thresh','var') || isempty(const_thresh)
    const_thresh = 1e-5;
end

approxErr = NaN(size(fits));
for i=1:length(fits)
    indsd = strmatch('sd_TF',fits{i}.pLabel,'exact');
    
    if range(fits{i}.data.yExp)>const_thresh        
        fits{i}.approxErr = 10.^fits{i}.p(indsd)./range(fits{i}.data.yExp);
    else
        fits{i}.approxErr = NaN;
    end
    approxErr(i) = fits{i}.approxErr;
end

[~,rf2] = sort(approxErr);

hist(approxErr(rf2),100)
set(gca,'FontSize',18)
title('Performance')
ylabel('No. of fits')
xlabel('approximation error [1/range]')
set(gcf,'Position',[500.0000  342.0000  900.6000  420.0000])
