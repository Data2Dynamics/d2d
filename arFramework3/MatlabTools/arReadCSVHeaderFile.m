% [header, data, dataCell] = arReadCSVHeaderFile(filename, delimiter, [quoted])
%
% Reads in CSV
%
% quoted [false] if true also quoted text (delimiters, white-space) can be
%                read in
% 
% [header, data] = arReadCSVHeaderFile('Test.csv',';',true)

function [header, data, dataCell] = arReadCSVHeaderFile(filename, delimiter, quoted)

fid = fopen(filename, 'r');

if(quoted)
    strquoted = '%q';
else
    strquoted = '%s';
end

arFprintf( 6, 'Reading CSV file ' );
C = textscan(fid,'%s\n',1,'Delimiter','');
C = textscan(C{1}{1},strquoted,'Delimiter',delimiter);
C = C{1};
header = C';

data = nan(0, length(header));
dataCell = cell(0, length(header));

rcount = 1;
C = textscan(fid,'%s\n',1,'Delimiter','');
while(~isempty(C{1}))
    C = textscan(C{1}{1},strquoted,'Delimiter',delimiter);
    C = strrep(C{1}',',','.');
    for j=1:length(header)
        if(j>length(C))
            data(rcount,j) = str2double('');
            dataCell{rcount,j} = '';
        else
            data(rcount,j) = str2double(C{j});
            dataCell{rcount,j} = C{j};
        end
    end
    C = textscan(fid,'%s\n',1,'Delimiter','');
    rcount = rcount + 1;
end
arFprintf( 6, '[ OK ]\n' );

fclose(fid);