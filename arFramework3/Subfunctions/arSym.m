% s = arSym(str)
%
%   This function implements the functionality of sym() in a
%   version-independent manner.
% 
%   This function converts a string into a symbolic variable.
%   It is required to support matlab's deprecated functionality
%   of using sym.m with arbitrary strings.
%
%   Example: sym('a+b') worked before R2018a but was replaced by str2sym('a+b')
%   starting from R2018a.


% The performance of verLessThan compared to version.m was compared using
% the following commands: [Result = alost the same speed]
% tic
% for i=1:1000
%     tmp = verLessThan('matlab','9.4');
% end
% toc
%
% tic
% for i=1:1000
%     v = version;
%     tmp = str2num(v(1:3))<9.4;
% end
% toc

function s = arSym(str)

try
    if verLessThan('matlab','9.4')
        ws=warning('query','symbolic:sym:sym:DeprecateExpressions');
        warning('off','symbolic:sym:sym:DeprecateExpressions');
        s = sym(str);
        warning(ws.state,'symbolic:sym:sym:DeprecateExpressions');
    elseif isempty(str)
        s = sym([]);
        s = reshape(s,size(str)); % required e.g. for arSym(ones(0,10))
    else
        if ischar(str)
            s = arSym(str);
        elseif iscell(str)
%             try
%                 try
%                     s = cell2sym(str);  % this only works for a cell of numbers
%                 catch
%                     str = regexprep(str,'\[(\d+)\]','_$1_')
%                     s = cell2sym(str);  % this only works for a cell of numbers
%                 end
%                 
%             catch 
%                 s = cell(size(str));                
                for i=1:length(str(:))
                    s(i) = arSym(str{i});
                end
                s = reshape(s,size(str));
%                 if length(s)==1
%                     s = s{1};
%                 end
%             end
        else
            s = sym(str);
        end
    end
catch ERR
    str
    rethrow(ERR)
end
