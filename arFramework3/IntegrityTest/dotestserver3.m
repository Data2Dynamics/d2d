% ausführen in Documents/MATLAB/d2d_testing/d2d
addpath('~/Documents/MATLAB/d2d_testing/d2d/arFramework3');
addpath('~/Documents/MATLAB/d2d_testing/d2d/arFramework3/Examples');

p = gcp('nocreate');
if isempty(p)
    parpool('local', 4);
end

doTests()

diaryfile = fileread('matlab_output.txt');

if isempty(strfind(diaryfile, 'Testing complete!'))
    test_status = 'tests not finished';
elseif isempty(strfind(diaryfile, '0 tests failed'))
   % msg = extractAfter(diaryfile, 'Testing complete!');
    whichTests = extractBetween(diaryfile, 'The following tests failed: ', '*');
    whichTests = strtrim(regexprep(whichTests{1},'\n+',''));
    test_status = ['errors in ' whichTests];
    system('cp ~/monitor_allknecht/failed.png ~/monitor_allknecht/current_status.png')
elseif isempty(strfind(diaryfile, '0 tests skipped'))
    test_status = 'some tests were skipped';
else
    test_status = 'working';
    system('cp ~/monitor_allknecht/working.png ~/monitor_allknecht/current_status.png')
end

logfileID = fopen('log.txt', 'a');
fprintf(logfileID, [datestr(datetime) '\t' test_status '\n']);
fclose(logfileID);
quit;
%sendmail('adrianhauber@gmail.com', 'd2d status alert', msg, '~/Documents/MATLAB/arFramework3/Examples/diary');