function str = arExtendStr(str, n)
nd = n - length(str);
S = ' ';
str = [str S(ones(1,nd))];
