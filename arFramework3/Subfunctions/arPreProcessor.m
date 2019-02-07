% function fid = preprocess(fid)
%
% The preprocessor can be used to incorporate basic compile time logic by 
% defining or un-defining labels using #define and #undefine. Code blocks 
% can then be bracketed by #ifdef, #else, #endif statements in order to 
% selectively activate/deactivate specific parts of the model. 
% 
% This function will apply the preprocessor on a file. This will look for
% #define/#undefine/#ifdef/#ifndef/#else blocks. The result is a file where
% defined names have been replaced by their definitions and blocks guarded
% by ifdef's or ifndef's are conditionally kept or discarded based on the
% state of that flag at the location of parsing.
%
%   fid     Handle to be preprocessed.
%
% Examples:
% #define HAS_HILL
% 
% #ifdef HAS_HILL
%     // This code will be evaluated
% #else
%     // This code will not be evaluated
% #end
% 
% #undefine HAS_HILL
% 
% #ifdef HAS_HILL
%     // This code will not be evaluated
% #else
%     // This code will be evaluated
% #end
% 
% #ifndef HAS_HILL
%     // This code will be evaluated
% #else
%     // This code will not be evaluated
% #end
% 
% //#define TEST
% #ifdef TEST
%     // This code will not be evaluated
% #end
% 
% #define VOLUME 5
% // Any instances of VOLUME after this define will be replaced by 5
%

function fid = preprocess(fid)

    global ar;
    
    % List of definitions that have been defined
    activeDefines = {};
    namedDefines = struct;
       
    % Defines that are relevant for the current codepiece
    defineStack = {};
    defineRequirement = [];

    C{1} = 'dummy';
    out = '';
    writing = 1;
    done = 0;
    directives = { '#ifdef', '#ifndef', '#else', '#endif', '#define', '#def', '#undefine', '#undef', '#end' };
    outStr = [];
    
    % Find newlines
    fid.str = regexprep(fid.str, '\/\/[^\n]+\n', '\n');
    newlines = regexp(fid.str, '\n|\r\n|\r');
    newlines = [newlines numel(fid.str)];
    statusChange = 1;
    while( ~done )
        %[line, fid] = arTextScan(fid, '%[^\n]', 'CommentStyle', ar.config.comment_string );
        nl = newlines( find( newlines > fid.pos, 1 ) );
        
        line = fid.str( fid.pos : nl );
        if ( fid.pos < numel(fid.str) )
            outLine = [];

            % Find whitespace in line
            a = 1;
            cur = 1;
            ws = regexp(line, '\s');
            ws = [ws numel(line)+1];
            while ( a <= numel( ws ) )
                [cur, outLine, parseBlock, a] = getParseBlock( line, cur, ws, outLine, a, writing, directives, namedDefines );

                % Preprocessor directives!
                if strcmpi( parseBlock, '#ifdef' )
                    [cur, outLine, parseBlock, a] = getParseBlock( line, cur, ws, outLine, a, 0, directives, namedDefines ); % Parse, but never write stuff that belongs to the precompiler directive
                    if isempty( parseBlock )
                        error( '#ifdef without condition on line %d: %s', fid.nlines, line );
                    end
                    defineStack{end+1} = parseBlock;
                    defineRequirement(end+1) = 1;           %  1 indicates that the define must be active for this block to be transcribed
                    statusChange = 1;
                end
                % Preprocessor directives!
                if strcmpi( parseBlock, '#ifndef' )
                    [cur, outLine, parseBlock, a] = getParseBlock( line, cur, ws, outLine, a, 0, directives, namedDefines ); % Parse, but never write stuff that belongs to the precompiler directive
                    if isempty( parseBlock )
                        error( '#ifndef without condition on line %d: %s', fid.nlines, line );
                    end
                    defineStack{end+1} = parseBlock;
                    defineRequirement(end+1) = -1;          %  -1 indicates that the define may not be active for this block to be transcribed
                    statusChange = 1;
                end            
                if strcmpi( parseBlock, '#else' )
                    defineRequirement(end) = -1 * defineRequirement(end);   % Else inverts the criterion
                    statusChange = 1;
                end
                if ( strcmpi( parseBlock, '#endif' ) || ( strcmpi( parseBlock, '#end' ) ) )
                    defineStack(end) = [];
                    defineRequirement(end) = [];
                    statusChange = 1;
                end

                if strcmpi( parseBlock, '#define' ) || strcmpi( parseBlock, '#def' )
                    [cur, outLine, parseBlock, a] = getParseBlock( line, cur, ws, outLine, a, 0, directives, namedDefines ); % Parse, but never write stuff that belongs to the precompiler directive
                    [cur, outLine, parseBlock2, a] = getParseBlock( line, cur, ws, outLine, a, 0, directives, namedDefines ); % Parse, but never write stuff that belongs to the precompiler directive
                    
                    if ~isempty( parseBlock2 )
                        parseBlock2 = arStrip(line(ws(a-2):end));
                        a = a + numel( parseBlock2 );
                    end
                    
                    % Is this being parsed or not? Defines can be conditional!
                    if ( writing )
                        activeDefines = union( activeDefines, parseBlock );
                        if ( parseBlock2 )
                            namedDefines.(parseBlock) = parseBlock2;
                        end
                        statusChange = 1;
                    end
                end
                if strcmpi( parseBlock, '#undefine' ) || strcmpi( parseBlock, '#undef' )
                    [cur, outLine, parseBlock, a] = getParseBlock( line, cur, ws, outLine, a, 0, directives, namedDefines ); % Parse, but never write stuff that belongs to the precompiler directive
                    if isempty( parseBlock )
                        error( '#ifdef without condition on line %d: %s', fid.nlines, line );
                    end
                    
                    % Is this being parsed or not? Defines can be conditional!
                    if ( writing )
                        activeDefines = setdiff( activeDefines, parseBlock );
                        if isfield( namedDefines, parseBlock )
                            namedDefines = rmfield( namedDefines, parseBlock );
                        end
                        statusChange = 1;
                    end
                end            

                if ( statusChange )
                    statusChange = 0;
                    
                    % Which defines have to be active for the current code block?
                    mustBeActive = defineStack((defineRequirement == 1));
                    mustBeInactive = defineStack((defineRequirement == -1));

                    % Any incompatible defines? (Unreachable code)
                    inco = intersect( mustBeActive, mustBeInactive );
                    if ( numel( inco ) > 0 )
                        incomp = sprintf( '%s ', inco{:} );
                        error( 'Incompatible #define / #undefine %s at line %d: %s', incomp, fid.nlines, line );
                    end

                    % Are we going to transcribe where we are?
                    writing = ( numel( intersect( activeDefines, mustBeActive ) ) == numel( mustBeActive ) ) && ... 
                              ( numel( intersect( activeDefines, mustBeInactive ) ) == 0 );
                end
            end
            outStr = sprintf( '%s%s', outStr, outLine );
        else
            done = 1;
        end
        
        % Move position
        if ~isempty( nl )
            fid.pos = nl + 1;
        end
    end
    
    fid.str = outStr;

    % Reset position
    fid.pos = 1;
    fid.nlines = 0;    
    
% This function is used by the preprocessor to parse a block of code
% It takes a line, a character position pointer (cur) w.r.t. line, a white space
% position counter (a), a list of whitespace locations (ws) and an outLine (output) 
% which is being assembled for output. The writing flag gives whether the
% output has to be written to the outline or not (depends on defines
% currently active). Directives contains a list of valid preprocessor
% directives, which will not be streamed to the output.
function [cur, outLine, parseBlock, a] = getParseBlock( line, cur, ws, outLine, a, writing, directives, namedDefines )

    % This will remove the leading and trailing whitespace. Note
    % that for the output we wish to preserve this whitespace. The
    % preprocessor should not touch anything (including whitespace)
    % that isn't preprocessed out
    parseBlock = '';
    while ( isempty( parseBlock ) && (a <= numel(ws) ) )
        block = line( cur : ws(a)-1 );
        parseBlock = strtrim(block);
        cur = ws(a);
        
        % Does the parse block start with a #?
        if ( numel( parseBlock ) > 0 )
            if ( parseBlock(1) == '#' )
                if ~ismember( parseBlock, directives )
                    valid = sprintf( '%s ', directives{:} );
                    error( 'Invalid preprocessor directive %s\nValid directives are: %s', parseBlock, valid );
                else
                    % Found a valid directive. Don't write this to the
                    % output stream
                    writing = 0;
                end
            else
                % Replace defined identifiers
                % Make sure long names get replaced first
                names = fieldnames( namedDefines );
                [~,I] = sort(cellfun(@length, names), 'descend');
                names = names(I);
                escapedNames = strcat( '\<', names, '\>' );
                
                done = false;
                iter = 0;
                while ( ~done )
                    newblock = block;
                    if ~isempty( names )
                        for b = 1 : numel( names )
                            newblock = regexprep( newblock, escapedNames{b}, namedDefines.(names{b}) );
                        end
                    end
                    if ( strcmp( newblock, block ) )
                        done = true;
                    else
                        block = newblock;
                    end
                    iter = iter + 1;
                    if ( iter > 100 )
                        error( 'Failure to run preprocessor. Do you have a recursive #define?' );
                    end
                end
            end
        end
        
        if ( writing )
            outLine = [outLine block]; %#ok
        else
            outLine = [outLine regexprep(block, '\S', ' ')]; %#ok
        end
        a = a + 1;
    end
