% set silent to false for debugging
silent = true;

fprintf( 2, 'TEST FOR IMPORT OF PETAB TOY MODELS\n' );

cases = {'0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008',...
    '0009', '0010', '0011', '0012', '0013', '0014','0015'};
Ncases = numel(cases);
try
    parpool(4)
end

for i = 1:Ncases
    cd(cases{i})
    fprintf( 2, ['\n\nCase ' cases{i} '...\n'] );
    try
        arInit
        arImportPEtab({'_model','_observables','_measurements','_conditions','_parameters'})
        
        arSimu(true,true,true)
        arCalcMerit
        [~,scores] = arGetMerit;
        
        chi2(i) = scores.chi2_res;
        llh(i) = scores.loglik/(-2);
        
        % check simulations
        simus = tdfread('_simulations.tsv');
        simus = struct2table(simus);
        simus2 = simus;
        
        simus2.simulation = NaN(size(simus2.simulation));
        
        q = 1;
        while q <= size(simus2,1)
            myobs = simus2.observableId(q,:);
            mycond = simus2.simulationConditionId(q,:);
            mytime = simus2.time(q);
            
            dataid = find(ismember({ar.model.data.name}, mycond));
            obsid = find(ismember(ar.model.data(dataid).y, myobs));
            
            timeid = ar.model.data(dataid).tExp == mytime;
            mysimus = ar.model.data(dataid).yExpSimu(timeid, obsid);
            simus2.simulation(q:q+size(mysimus,1)-1) = mysimus;
            
            qq = 1;
            if size(mysimus,1) > 0
                qq = size(mysimus,1);
            end
            q = q + qq;
        end
        
        abs(simus.simulation - simus2.simulation);
        SimuDiff(i) = sum(abs(simus.simulation - simus2.simulation));
        SolTable = ReadOutSolutionFile(cases{i});
        
        chi2Solution(i) = SolTable.Value(SolTable.Variable=='chi2:');
        llhSolution(i) = SolTable.Value(SolTable.Variable=='llh:');
        
        TolChi2(i) = SolTable.Value(SolTable.Variable=='tol_chi2:');
        TolLLH(i) = SolTable.Value(SolTable.Variable=='tol_llh:');
        TolSimu(i) = SolTable.Value(SolTable.Variable=='tol_simulations:');
        Error{i} = 'none';
        ErrorFile{i} = 'none';
        ErrorLine{i} = 'none';
    catch Err
        Error{i} = Err.message;
        ErrorFile{i} = Err.stack.file;
        ErrorLine{i} = Err.stack.line;
        chi2(i) = nan;
        llh(i) = nan;
        SimuDiff(i) = nan;
        
        chi2Solution(i) = nan;
        llhSolution(i) = nan;
        
        TolChi2(i) = nan;
        TolLLH(i) = nan;
        TolSimu(i) = nan;
    end
    
    cd ..
end

Chi2Diff = abs(chi2-chi2Solution);
LLHDiff = abs(llh-llhSolution);

SimuCheck = SimuDiff<TolSimu;
Chi2Check = Chi2Diff<TolChi2;
LLHCheck = LLHDiff<TolLLH;

Working = double(SimuCheck.*Chi2Check.*LLHCheck);

if ~silent
    Table = table(cases',Working',SimuCheck',Chi2Check',LLHCheck',SimuDiff',Chi2Diff',LLHDiff',chi2',llh',Error',ErrorFile',ErrorLine','VariableNames',{'Case','Working','Simu','Chi2','LLH','SimuDiff','Chi2Diff','LLHDiff','chi2','llh','Error','ErrorFile','ErrorLine'})
end

if sum(Working) == numel(Working)
    fprintf( 2, 'PASSED\n' );
else
    fprintf( 2, 'Errors in test case(s) %s\n', strjoin(cases(logical(Working)),', '));
    error( 'FAILED');
end

function solution = ReadOutSolutionFile(caseName)

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, 2; 5, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["Variable", "Value"];
opts.VariableTypes = ["string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "Variable", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Variable", "EmptyFieldRule", "auto");
% Import the data
solution = readtable(['_' caseName '_solution.yaml'], opts);
end


function opts = delimitedTextImportOptions(varargin)
%Create options for importing delimited text data
% opts = delimitedTextImportOptions('Prop1',val1,'Prop2',val2,...) creates
%        options for importing a delimited text file.
%
%   DelimitedTextImportOptions Properties:
%
%                     DataLines - The lines in the text file where the data is located.
%             VariableNamesLine - Where the variable names are located
%                RowNamesColumn - Where the row names are located
%             VariableUnitsLine - Where the variable units are located
%      VariableDescriptionsLine - Where the variable descriptions are located
%                 VariableNames - Names of the variables in the file
%         SelectedVariableNames - Names of the variables to be imported
%                 VariableTypes - The import types of the variables
%               VariableOptions - Advanced options for variable import
%         PreserveVariableNames - Whether or not to convert variable names
%                                 to valid MATLAB identifiers.
%               ImportErrorRule - Rules for interpreting nonconvertible or bad data
%                   MissingRule - Rules for interpreting missing or unavailable data
%              ExtraColumnsRule - What to do with extra columns of data that appear
%                                 after the expected variables
%     ConsecutiveDelimitersRule - What to do with consecutive
%                                 delimiters that appear in the file
%         LeadingDelimitersRule - What to do with delimiters at the beginning of a
%                                 line
%                 EmptyLineRule - What to do with empty lines in the file
%                    Delimiter  - Symbol(s) indicating the end of data fields in the
%                                 file
%                    Whitespace - Characters to be treated as whitespace.
%                    LineEnding - Symbol(s) indicating the end of a line in the file
%                  CommentStyle - Symbol(s) designating text to ignore
%                      Encoding - Text encoding of the file to be imported
%                  NumVariables - The number of variables to import
%
% See Also
%   detectImportOptions, readtable,
%   matlab.io.text.DelimitedTextImportOptions

% Copyright 2018-2019 MathWorks, Inc.

    opts = matlab.io.text.DelimitedTextImportOptions(varargin{:});
end
