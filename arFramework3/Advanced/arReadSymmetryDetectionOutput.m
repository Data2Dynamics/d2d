% [group,variable,trsf] = arReadSymmetryDetectionOutput(fileResult) 
% 
%   Reads the output of Benjamin Mertk's Symmetry-Detection-Tool.
%   This tool has to be called externally (e.g. via command line).
% 
% Example:
% 
% [group,variable,trsf] = arReadSymmetryDetectionOutput('jak2_stat5_feedbacks__model.csv_result.txt')

function  [group,variable,trsf] = arReadSymmetryDetectionOutput(fileResult)

fid = fopen(fileResult,'r');

variable = cell(0);
trsf = cell(0);
group = cell(0);
try
    str = textscan(fid, '%s\t%s\t%s\t%s\n', 1,'Delimiter','\t');
catch
    error('Result file %s not found or wrong format. Does python script symetry detection work? ',fileResult);
end
while ~isempty(str{1})
    type = str{2}{1};
    if strcmp(type,'Type: scaling')==1
        group{end+1} = str{1}{1};
        variable{end+1} = str{3}{1};
        trsf{end+1} = str{4}{1};
    end
    str = textscan(fid, '%s\t%s\t%s\t%s\n', 1,'Delimiter','\t');
end
fclose(fid);
