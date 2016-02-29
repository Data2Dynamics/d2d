function val = get(m,prop)

if ~exist('prop','var') 
    prop = [];
end

if length(m) > 1
    for i=1:length(m)
        a = get(m(i),prop);
        val(i) = {a};
    end
    return
end

if isempty(prop)
   val = fields(m); 
   return
end

try
    val = m.(prop);
catch
    error([mfilename,': Can not get property ', prop, '.'])
end

