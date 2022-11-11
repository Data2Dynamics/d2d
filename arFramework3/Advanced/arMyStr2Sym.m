% out = arMyStr2Sym(s) 
%
% Cast to symbolic expression based on input type. Required since
% MATLAB2018a due to changes in sym(#).
%
%   s        String which is converted to symbolic expression
%   out      Symbolic expression 

function out = arMyStr2Sym(s)

persistent matver % keeping the value from the last call
if isempty(matver)
    matver = ver('MATLAB');  % calling this function every time is too time-consuming
end
    % The explicit cast to string is necessary for MATLAB R2017a at least,
    % otherwise double will convert the string on a char by char basis.
    if(str2double(matver.Version) >= 9.1)
        if(isa(s,'double'))
            out = sym(s);
        elseif(isa(s,'char'))
            % The use of symengine for evaluating strings resulted in
            % problematic symbolic expressions which were not properly
            % usable for later substitions
            % Therefore the code was reverted to the standard sym function
            % from Matlab
            % This works fine for Matlab version 2018b
            try
                out = sym(s);
            catch
                try
                    out = evalin(symengine,s);
                catch ERR
                    error(('Failed to cast to symbolic expression'))
                end
            end        
        elseif(isa(s,'cell'))
            size_s = size(s);
            out = sym(zeros(size_s));
            for i = 1:size_s(1)
                for j = 1:size_s(2)
                    if ~strcmp(s{i,j}, 'matrix([])')
                        out(i,j) = evalin(symengine,s{i,j});
                    end
                end
            end
        else
            try
                out = sym(s);
            catch
                try
                    out = evalin(symengine,s);
                catch ERR
                    error(('Failed to cast to symbolic expression'))
                end
            end
        end
    else
        out = sym(s);
    end
end