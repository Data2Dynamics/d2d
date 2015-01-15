% Workaround for bug in Matlab's copyfile in R2014b running on Mac Yosemite

function copyfile(source,target)

system(['cp -r -p ' source ' ' target]);
