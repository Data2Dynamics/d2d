
function str = strRepsCellLines(str)
str = upper(str);
str = strrep(str, ' ', '');
str = strrep(str, '-', '');
str = strrep(str, '_', '');
str = strrep(str, '.', '');
str = strrep(str, ',', '');
str = strrep(str, ';', '');
str = strrep(str, ':', '');
str = strrep(str, '/', '');
str = strrep(str, '\', '');
str = strrep(str, 'NCI', '');
str = strrep(str, 'NIH', '');
