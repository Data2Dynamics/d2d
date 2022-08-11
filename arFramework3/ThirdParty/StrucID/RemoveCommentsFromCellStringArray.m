function trimmedTxt=RemoveCommentsFromCellStringArray(strArray)
% removes comments from a cell array of strings
indices=regexp(strArray,'%','once');
% remove lines that start with a comment
isOne=cellfun(@(x) isequal(x,1),indices);
strArray(isOne)=[]; indices(isOne)=[];

trimmedTxt=cellfun(@(x,y)strtrim(x(1:y-1)),strArray,indices,'UniformOutput',false);
indexEmpty=cellfun(@isempty,trimmedTxt);
trimmedTxt(indexEmpty) = strArray(indexEmpty);