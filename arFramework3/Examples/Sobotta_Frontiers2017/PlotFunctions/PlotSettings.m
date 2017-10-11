% Which type of error visualization are we using?
% 1 = bars on the data
% 2 = lines on the model
% 3 = patches
global errorBarMode;
errorBarMode = 2;
 
global allTicks;
allTicks = 1;

% Make sure the simulations are up to date
arSimu(false, true, true);
arCalcMerit(false);
arGetMerit;