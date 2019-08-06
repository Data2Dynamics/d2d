% Setup
% arSave

%% adding path to folder 'TransientFunction_library'
% addpath('E:\clemens\Repositories\d2d\arFramework3\Examples\ToyModels\TransientFunction')

%% Analysis
arLoad

if exist('Analysis.log','file')
    delete('Analysis.log');
end
diary('Analysis.log');

fits = arApproximateTimeCoursesByTransientFunction2;
WriteMatlabFunction(fits)

diary off
