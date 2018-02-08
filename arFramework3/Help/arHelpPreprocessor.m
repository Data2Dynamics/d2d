% Preprocessor directives can be used in both types of def files. They make it 
% easier to quickly activate or deactivate parts of models and allow
% splitting complicated expressions into multiple blocks.
%
% The preprocessor can be used to incorporate basic compile time logic by 
% defining or un-defining labels using #define and #undefine. Code blocks 
% can then be bracketed by #ifdef, #else, #endif statements in order to 
% selectively activate/deactivate specific parts of the model. An example use:
% 
%   Examples:
%     #define HAS_HILL
%     #ifdef HAS_HILL
%         // This code will be evaluated
%     #else
%         // This code will not be evaluated
%     #end
%   
%     #undefine HAS_HILL
%     #ifdef HAS_HILL
%         // This code will not be evaluated
%     #else
%         // This code will be evaluated
%     #end
% 
%     #ifndef HAS_HILL
%         // This code will be evaluated
%     #else
%         // This code will not be evaluated
%     #end
% 
%     //#define TEST
%     #ifdef TEST
%         // This code will not be evaluated
%     #end
%
% One can also define subexpressions using the preprocessor, but one must
% be careful that the preprocessor acts as though they were find/replace
% operations and has no knowledge of mathematical rules. It is therefore
% prudent to always pad any formula define with round brackets.
% 
%   Examples:
%     #define FUNC1 ( is*a*test )
%     #define MYINPUT ( this * FUNC1 )
%
%   // Here the input would be defined as ( this * ( is*a*test ) )
%   TestInput   C   ng/ml   conc.   "MYINPUT"

help arHelpPreprocessor
