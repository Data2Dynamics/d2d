function retInt = optimizeWith(int, optType, varargin)
% Start a EvA2 optimization run with specific optimizer parameter settings.
%       optimize(interface, optType, [, outputFilePrefix ] [, memName, memVal]* )
% where
%       interface: instance of JEInterface
%       optType: integer indicating the type of the optimization strategy
%       to use.
%       resultFilePrefix: (optional) char prefix for an optional verbose
%           output file
%       memName: character name of a member of the optimizer
%       memVal: value to set for the member.
% Note that the parameter settings will not be stored and remain valid only for 
% the current optimization run.


% This just reads the parameters, sets them as interface members and then
% calls the standard optimize method.

if (int.finished == 0) 
    error('please wait for the current run to finish');
end
if (nargin < 2) 
    error('invalid number of arguments!');
end

int.optParams = [];
int.optParamValues = [];

if (nargin == 2) 
    % standard case, just call optimize
    retInt = optimize(int, optType);
elseif (nargin == 3)
    % standard case, just call optimize
    retInt = optimize(int, optType, varargin{1});
else
    if (mod(nargin, 2) == 1) % odd arguments! There is an outputfilePrefix!
        nextMem = 2;
        output = varargin{1};
    else
        nextMem = 1;
        output = '';
    end
    pairCnt = (nargin-nextMem-1)/2;
    parValues(1) = {-1};
%    parArr = ones(2*pairCnt,1); % parameter array
    for i=1:pairCnt
        % load parameter/value pairs into an array
        parNames(i) = cellstr(varargin{nextMem});
        parValues(i) = {varargin{nextMem+1}};

%        if (~isstr(parNames(i))) 
%            error('invalid argument, member names must be char');
%        end
        nextMem = nextMem+2;
    end
    int.optParams = parNames;
    int.optParamValues = parValues;
    retInt = optimize(int, optType, output);
    retInt.optParams = [];
    retInt.optParamValues = [];
end
