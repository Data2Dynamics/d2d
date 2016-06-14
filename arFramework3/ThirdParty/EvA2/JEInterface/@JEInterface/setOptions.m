function int = setOptions(int, usrOpts)
% Set the optimization options for the interface. The existing options are
% overwritten by the given setings.
%   parameters:
%       int:    an interface instance
%       usrOpts:    a JE options structure

fn=fieldnames(usrOpts);

options = int.opts;

try
    for i=1:length(fn)
        % make sure all option fields and values are valid
        % fn(i)
        % ischar(fn(i))
        makeOptions(int, char(fn(i)), usrOpts.(char(fn(i))));
        options.(char(fn(i))) = usrOpts.(char(fn(i)));
    end
catch ME
    error('invalid option "%s"... check makeOptions to learn about available options', char(fn(i)));

end

int.opts = options;