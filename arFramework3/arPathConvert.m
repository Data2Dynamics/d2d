function str = arPathConvert(str)
str = strrep(str, ' ', '\ ');
str = strrep(str, '(', '\(');
str = strrep(str, ')', '\)');
