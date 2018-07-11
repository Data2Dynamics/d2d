function out = arMyStr2Sym(s)

%   Cast to symbolic expression based on input type. Required since
%   MATLAB2018a due to changes in sym(#).

persistent matver % keeping the value from the last call
if isempty(matver)
    matver = ver('MATLAB');  % calling this function every time is too time-consuming
end

    if(double(matver.Version) >= 9.4)
        if(isa(s,'double'))
            out = sym(s);
        elseif(isa(s,'char'))
            out = evalin(symengine,s);
        elseif(isa(s,'cell'))
            size_s = size(s);
            out = sym(zeros(size_s));
            for i = 1:size_s(1)
                for j = 1:size_s(2)
                    out(i,j) = evalin(symengine,s{i,j});
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