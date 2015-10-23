%   arPrint
%   arPrint(3:4)
%   arPrint('sd')
%   arPrint({'para1','para2'})
% 
% print parameter values
% 
%   js      Indices of the parameters to be displayed. 'All' refers to all
%           parameters
%           (see Example below)
%
% optional arguments:
%           'initial'                   - only show initial conditions
%           'fitted'                    - only show fitted parameters
%           'constant'                  - only show constants
%           'fixed'                     - only show fixed parameters
%           'dynamic'                   - only show dynamic parameters
%           'observation'               - only show non-dynamic parameters
%           'error'                     - only show error model parameters
%           'closetobound'              - show the parameters near bounds
%           'lb' followed by value      - only show values above lb
%           'ub' followed by value      - only show values below lb
%           Combinations of these flags are possible
% 
% Examples:
% arPrint('turn')
%            name                      lb       value       ub          10^value        fitted   prior
% #  22|DI | geneA_turn              |       -5      -0.41         +3 | 1      +0.39 |       1 | uniform(-5,3) 
% #  31|DI | geneB_turn              |       -5       -2.7         +3 | 1    +0.0022 |       1 | uniform(-5,3) 
% #  40|DI | geneC_turn              |       -5      -0.91         +3 | 1      +0.12 |       1 | uniform(-5,3) 
%      |   |                         |                                |              |         |      
% #  49|DI | geneD_turn              |       -5       -2.1         +3 | 1     +0.008 |       1 | uniform(-5,3) 
% #  58|DI | geneE_turn              |       -5       -2.5         +3 | 1     +0.003 |       1 | uniform(-5,3) 
% 
% arPrint(15:17)
% Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model
% 
% Example:
%            name                      lb       value       ub          10^value        fitted   prior
% #  15|D  | geneA_deg1              |       -2         -2         +3 | 1      +0.01 |       1 | uniform(-2,3) 
% #  16|D  | geneA_deg2              |       -2      +0.36         +3 | 1       +2.3 |       1 | uniform(-2,3) 
% #  17|D  | geneA_deg3              |       -2       -1.2         +3 | 1     +0.063 |       1 | uniform(-2,3) 
%
% Other examples:
% arPrint(1:4, 'constant')                  - Show the constants among the
%                                             first four parameters
% arPrint('all', 'fitted', 'dynamic')       - Show all fitted dynamic parameters
% arPrint('gene', 'observation')            - Show all observation parameters
%                                             containing the string "gene" in the name

function varargout = arPrint(js, varargin)
global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('js','var') || isempty(js))
    js = 1:length(ar.p);
elseif(islogical(js))
    js = find(js);
elseif(isnumeric(js))
    if(size(js,1)>1)
        js = js'; %should not be a row
    end
    if(sum(js==1 | js==0) ==length(js) && length(js)>1)
        js = find(js);
    end
elseif(ischar(js))
    if ( strcmp( js, 'all' ) )
        js = 1:length(ar.p);
    else
        js = find(~cellfun(@isempty,regexp(ar.pLabel,js)));
        if isempty(js)
            disp('Pattern not found in ar.pLabel');
            return;
        end
    end
elseif(iscell(js)) % cell of pNames
    [~,js] = intersect(ar.pLabel,js);
else
    error('Argument has to be a string or an array of indices.')
end

if(sum(isnan(js))>0 || sum(isinf(js))>0 || min(js)<1 || max(js-round(js))>eps)
    js %#ok
    warning('arPrint.m: argument js is not plausible (should be an array of indices).')
else
    if(size(js,1)~=1)
        js = js';
    end
end

pTrans = ar.p;
pTrans(ar.qLog10==1) = 10.^pTrans(ar.qLog10==1);
 
% Determine which parameters are close to their bounds
qFit = ar.qFit == 1;
ar.qCloseToBound(size(ar.qFit)) = 0;
ar.qCloseToBound(qFit) = (ar.p(qFit) - ar.lb(qFit)) < ar.config.par_close_to_bound*(ar.p(qFit) - ar.lb(qFit)) | ...
    (ar.ub(qFit) - ar.p(qFit)) < ar.config.par_close_to_bound*(ar.p(qFit) - ar.lb(qFit));

% Additional options
if ( nargin > 1 )
    opts = argSwitch( {'closetobound', 'initial', 'fixed', 'fitted', 'dynamic', 'constant', 'observation', 'error', 'lb', 'ub'}, varargin{:} );

    if ( opts.constant && opts.fitted )
        error( 'Incompatible flag constant and fitted' );
    end
    if ( opts.fixed && opts.fitted )
        error( 'Incompatible flag fixed and fitted' );
    end    
    if ( opts.constant && opts.fixed )
        error( 'Incompatible flag constant and fixed' );
    end
    if ( opts.observation && opts.dynamic )
        error( 'Incompatible flag observation and dynamic' );
    end    

    if ( opts.fixed )
        js = js( ar.qFit( js ) == 0 );
    end
    if ( opts.fitted )
        js = js( ar.qFit( js ) == 1 );
    end
    if ( opts.constant )
        js = js( ar.qFit( js ) == 2 );
    end
    if ( opts.dynamic )
        js = js( ar.qDynamic( js ) == 1 );
    end
    if ( opts.observation )
        js = js( ar.qDynamic( js ) == 0 );
    end
    if ( opts.error )
        js = js( ar.qError( js ) == 1 );
    end    
    if ( opts.initial )
        js = js( ar.qInitial( js ) == 1 );
    end
    if ( opts.ub )
        js = js( ar.p(js) < opts.ub );
    end
    if ( opts.lb )
        js = js( ar.p(js) > opts.lb );
    end
    if ( opts.closetobound )
        js = js( ar.qCloseToBound(js) );
    end
end

if nargout>0
    varargout{1} = js;
end

maxlabellength = max(cellfun(@length, ar.pLabel(js)));

fprintf('Parameters: # = free, C = constant, D = dynamic, I = initial value, E = error model\n\n');
printHead;
for j=1:length(js)
    printPar(js(j), ar.qCloseToBound(js(j)));
	if(mod(j,10)==0 && j<length(js))
		fprintf(['     |   | ' arExtendStr('', maxlabellength) ' |                                |              |         |      \n']);
	end
end

    function printHead
        strhead = ['     |   | ' arExtendStr('name', maxlabellength) ' | lb       value       ub        | 10^value      | fitted | prior\n'];
        strhead = strrep(strhead, '|', ' ');
        fprintf(strhead);
    end

    function printPar(j, qclosetobound)
        strdyn = ' ';
        if(ar.qDynamic(j))
            strdyn = 'D';
        end
        strerr = ' ';
        if(ar.qError(j))
            strerr = 'E';
        end
        strinit = ' ';
        if(ar.qInitial(j))
            strinit = 'I';
        end
        
        if(qclosetobound)
            outstream = 2;
        else
            outstream = 1;
        end
        if(ar.qFit(j)==2)
            fit_flag = 'C';
        else
            fit_flag = '#';
        end
        fprintf(outstream, '%s%4i|%s%s%s| %s | %+8.2g   %+8.2g   %+8.2g | %i   %+8.2g | %7i | %s \n', ...
            fit_flag, j, strdyn, strinit, strerr, arExtendStr(ar.pLabel{j}, maxlabellength), ar.lb(j), ar.p(j), ar.ub(j), ar.qLog10(j), pTrans(j), ar.qFit(j), priorStr(j));
        
    end

    function str = priorStr(j)
        if(ar.type(j) == 0)
            str = sprintf('uniform(%g,%g)', ar.lb(j), ar.ub(j));
        elseif(ar.type(j) == 1)
            str = sprintf('normal(%g,%g^2)', ar.mean(j), ar.std(j));
        elseif(ar.type(j) == 2)
            str = sprintf('uniform(%g,%g) with soft bounds', ar.lb(j), ar.ub(j));
        elseif(ar.type(j) == 3)
            str = sprintf('L1(%g,%g)', ar.mean(j), ar.std(j));
        end
    end
end

function [opts] = argSwitch( switches, varargin )

    for a = 1 : length(switches)
        opts.(switches{a}) = 0;
    end

    a = 1;
    while (a <= length(varargin))
        if ( max( strcmp( varargin{a}, switches ) ) == 0 )
            str = sprintf( 'Legal switch arguments are:\n' );
            str = [str sprintf( '%s\n', switches{:} ) ];%#ok<AGROW>
            error( 'Invalid switch argument was provided. Provided %s, %s', varargin{a}, str );
        else
            fieldname = varargin{a};
        end
        
        val = 1;
        if ( length(varargin) > a )
            if isnumeric( varargin{a+1} )
                val = varargin{a+1};
                a = a + 1;
            end
        end
        
        opts.(fieldname) = val;
        a = a + 1;
    end
end



