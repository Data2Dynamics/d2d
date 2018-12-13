% varargout = arEvaluate(fun,varargin)
% 
% Wrapper function which can make your code safer by preventing unintended 
% changes in the global ar variable if an unhandled error occurs.
% 
%   fun         function name or function handle evaluated by feval
%   varargin    these arguments are passed to the function fun
% 
% This function ensures that whenever an error occurs that the original
% global variable is restored. As a first step, a copy is made (but not
% using arDeepCopy).
%
% It prevents changes on the global ar which are not reset to the original
% state due to an error.
% 
% 
% Example:
% 1) Define a function altering ar and then generating an error
% 
% 2a) Call the function directly => ar is changed
% 2b) Call the fucntion via arEvaluate(fun,args...) => ar not changed:
% 
% function arTestFun 
% global ar; 
% ar.p = zeros(size(ar.p);
% error('Now any error should occur');
% 
% % function call
% ar.p(:)=1;
% arTestFun
% ar.p
% 
% % function call
% ar.p(:)=1;
% arEvaluate(@arTestFun)
% ar.p
% 
% 
% See also arPop, arPush

function varargout = arEvaluate(fun,varargin)

global ar
arIn = ar;
arIn.p = arIn.p + 0.0; % enforce (weakly) making a copy in the storage

try
    args = cell(1,nargout);
    switch nargout
        case 0
            feval(fun,varargin{:});
        case 1
            args{1} = feval(fun,varargin{:});
        case 2
            [args{1},args{2}] = feval(fun,varargin{:});
        case 3
            [args{1},args{2},args{3}] = feval(fun,varargin{:});
        case 4
            [args{1},args{2},args{3},args{4}] = feval(fun,varargin{:});
        case 5
            [args{1},args{2},args{3},args{4},args{5}] = feval(fun,varargin{:});
        case 6
            [args{1},args{2},args{3},args{4},args{5},args{6}] = feval(fun,varargin{:});
        otherwise
            error('Up to now arEvaluate is only implemented for up to six output variables.');
    end
            
catch err
    ar = arIn;
%     err.message = [err.message,' (fun=',char(fun),')'];
    rethrow(err);
end

varargout = cell(size(args));
for i=1:length(args)
    varargout{i} = args{i};
end

