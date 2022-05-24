% In this code, the read commands of arLoadModel are executed in the same
% way as in arLoadModel.m 
% 
% In this code example, the content ist written in an unmodified manner. 
% This Code can serve as template for custom model-def conversions. One has
% to add conversions at the right place.

function arConvertModelDef_Template(name)

global ar
ModelPath = 'Models/';

% load model from mat-file
if(~exist(ModelPath,'dir'))
    error('folder %s does not exist. Possible reason: Did you rename the result folder by hand?',ModelPath)
end
if strcmp(strrep(name,' ',''),name)~=1
    name
    error('File names should not contain empty spaces. Please remove it.');
end
if(~exist([ModelPath,  name '.def'],'file'))
    error('model definition file %s.def does not exist in folder %s', name, ModelPath)
end

fid = fopen([ModelPath,  name '.def'], 'r');
fidOut = fopen([ModelPath, name,'_logX.def'],'w');


matVer = arVer;

% DESCRIPTION
[str, fid] = arTextScan(fid, '%s', 1, 'CommentStyle', ar.config.comment_string);
if(isempty(strfind(str{1},'DESCRIPTION')))
    arParsingError( fid, 'parsing model %s for DESCRIPTION', ar.model(m).name);
end
myFprintf(fidOut,'%s\n',str{1}{:});

% read comments
[str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(str{1},'PREDICTOR'))
    myFprintf(fidOut,'"%q"\n',str{1}{:});
    [str, fid] = arTextScan(fid, '%q', 1, 'CommentStyle', ar.config.comment_string);
end
myFprintf(fidOut,'\n%q\n',str{1}{:});

% PREDICTOR
[C, fid] = arTextScan(fid, '%s %q %q %q %n %n\n',1, 'CommentStyle', ar.config.comment_string);
arValidateInput( C, 'predictor', 'identifier for independent variable', 'unit type', 'unit', 'label for plotting' );
myFprintf(fidOut,'%s %q %q %q %d %d\n',C{1}{1},C{2}{1},C{3}{1},C{4}{1},C{5},C{6});

myFprintf(fidOut,'\n');
% COMPARTMENTS
[C, fid] = arTextScan(fid, '%s %q %q %q %f\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'STATES'))
    if(~strcmp(C{1},'COMPARTMENTS'))
        myFprintf(fidOut,'%s %q %q %q %f \n',C{1}{1},C{2}{1},C{3}{1},C{4}{1},C{5});
        arValidateInput( C, 'compartments', 'compartment', 'unit type (i.e. V)', 'unit (i.e. pl)', 'label (i.e. "vol.")' );
    else
        myFprintf(fidOut,'%s\n',C{1}{1}); % 1st case
    end
    [C, fid] = arTextScan(fid, '%s %q %q %q %f\n',1, 'CommentStyle', ar.config.comment_string);
end

% STATES
myFprintf(fidOut,'\n%s\n',C{1}{1}); % STATES
[C, fid] = arTextScan(fid, '%s %s %s %s %s %n %q %n\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'INPUTS'))
    if ( strcmp( C{1}, 'REACTIONS' ) )
        arParsingError( fid,  'Missing field INPUTS. This section should be specified after STATES and before REACTIONS. See: "Setting up models"' );
    end
    
    arValidateInput( C, 'state', 'unique identifier', 'unit type (i.e. C)', 'unit (i.e. nM)', 'label for plots (i.e. "conc.")' );
    myFprintf(fidOut,'%s %q %q %q %s %n %q %n\n',C{1}{1},C{2}{1},C{3}{1},C{4}{1},C{5}{1},C{6},C{7}{1},C{8});
        
    [C, fid] = arTextScan(fid, '%s %s %s %s %s %n %q %n\n',1, 'CommentStyle', ar.config.comment_string);
end

% INPUTS
myFprintf(fidOut,'\n%s\n',C{1}{1}); % INPUTS
[C, fid] = arTextScan(fid, '%s %q %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'REACTIONS') && ~strcmp(C{1},'REACTIONS-AMOUNTBASED') && ~strcmp(C{1},'ODES'))
    if(~strcmp(C{1},''))
        arValidateInput( C, 'input', 'unique input name', 'unit type (i.e. C)', 'unit (i.e. "units/cell")', 'plain text label for plots ("conc.")', 'input function' );
        myFprintf(fidOut,'%s %q %q %q %q %q\n',C{1}{1},C{2}{1},C{3}{1},C{4}{1},C{5}{1},C{6}{1});
    end
    [C, fid] = arTextScan(fid, '%s %q %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
end

% REACTIONS/ODEs
myFprintf(fidOut,'\n%s\n',C{1}{1}); % REACTIONS
if(strcmp(C{1},'REACTIONS') || strcmp(C{1},'REACTIONS-AMOUNTBASED'))
    % Read single line (easier to trace errors back to their line)
    [ line, remainder, fid ] = readLine( fid, ar.config.comment_string );
    [str, remainder] = grabtoken( remainder, ' %s', 1 );
    myFprintf(fidOut,'%s \n',line);
    
    while(~strcmp(str{1},'INVARIANTS') && ~strcmp(str{1},'DERIVED'))
        [ line, remainder, fid ] = readLine( fid, ar.config.comment_string );        
        [str, remainder] = grabtoken( remainder, '%s', 1 );
        myFprintf(fidOut,'%s \n',line);
    end
elseif(strcmp(C{1},'ODES'))
    myFprintf(fidOut,'\n');
    [str, fid] = arTextScan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
    myFprintf(fidOut,'%q\n',str{1}{1});
    while(~strcmp(str{1},'INVARIANTS') && ~strcmp(str{1},'DERIVED'))
        [str, fid] = arTextScan(fid, '%q\n',1, 'CommentStyle', ar.config.comment_string);
        myFprintf(fidOut,'%q \n',str{1}{1});
    end
end


[C, fid] = arTextScan(fid, '%s %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
while(~strcmp(C{1},'CONDITIONS') && ~strcmp(C{1},'SUBSTITUTIONS') && ~strcmp(C{1},'OBSERVABLES'))
    myFprintf(fidOut,'%s %q %q %q %q\n',C{1}{1},C{2}{1},C{3}{1},C{4}{1},C{5}{1});
    [C, fid] = arTextScan(fid, '%s %q %q %q %q\n',1, 'CommentStyle', ar.config.comment_string);
end


% OBSERVABLES
if(strcmp(C{1},'OBSERVABLES'))
    myFprintf(fidOut,'\n%s\n',C{1}{1}); % OBSERVABLES
    [C, fid] = arTextScan(fid, '%s %q %s %s %n %n %s %s\n',1, 'CommentStyle', ar.config.comment_string);
    while(~strcmp(C{1},'ERRORS'))
        myFprintf(fidOut,'%s %s %q %q %n %n %q %q\n',C{1}{1},C{2}{1},C{3}{1},C{4}{1},C{5},C{6},C{7}{1},C{8}{1});
        [C, fid] = arTextScan(fid, '%s %s %s %s %n %n %s %s\n',1, 'CommentStyle', ar.config.comment_string);
    end
    
    % ERRORS
    myFprintf(fidOut,'\n%s\n',C{1}{1}); % ERRORS
    [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
    while(~(strcmp(C{1},'CONDITIONS') || strcmp(C{1},'SUBSTITUTIONS')))
        myFprintf(fidOut,'%s %s\n',C{1}{1},C{2}{1});
        [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
    end
end
myFprintf(fidOut,'\n%s\n',C{1}{1}); % CONDITIONS or SUBSTITUTIONS

% SUBSTITUTIONS (beta)
substitutions = 0;
if ( strcmp(C{1},'SUBSTITUTIONS') )
    if(str2double(matVer.Version)>=8.4)
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
    else
        [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16);
    end
    myFprintf(fidOut,'%s %q\n',C{1}{1},C{2}{1});
    substitutions = 1;
    
    % Fetch desired substitutions
    while(~isempty(C{1}) && ~strcmp(C{1},'CONDITIONS'))
        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %q\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
        myFprintf(fidOut,'%s %q\n',C{1}{1},C{2}{1});
    end
end

% CONDITIONS
if(str2double(matVer.Version)>=8.4)
    [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
else
    [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16);
end
myFprintf(fidOut,'%s %s\n',C{1}{1},C{2}{1});

if ( substitutions )
    % Fetch desired conditions
    while(~isempty(C{1}) && ~(strcmp(C{1},'PARAMETERS') || strcmp(C{1}, 'RANDOM')))
        myFprintf(fidOut,'%s %s\n',C{1}{1},C{2}{1});
        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
    end
    
else
    % Old code path
    while(~isempty(C{1}) && ~(strcmp(C{1},'PARAMETERS') || strcmp(C{1}, 'RANDOM')))
        myFprintf(fidOut,'%s %s\n',C{1}{1},C{2}{1});
        if(str2double(matVer.Version)>=8.4)
            [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
        else
            [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string, 'BufSize', 2^16-1);
        end
    end
end


if ( strcmp(C{1}, 'RANDOM' ) )
    myFprintf(fidOut,'%s\n',C{1});
    [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
    myFprintf(fidOut,'%s %s\n',C{1}{1},C{2}{1});
    while(~isempty(C{1}) && ~strcmp(C{1},'PARAMETERS'))
        [C, fid] = arTextScan(fid, '%s %s\n',1, 'CommentStyle', ar.config.comment_string);
        myFprintf(fidOut,'%s %s\n',C{1}{1},C{2}{1});
    end
end

% PARAMETERS
[C, fid] = arTextScan(fid, '%s %f %n %n %f %f\n',1, 'CommentStyle', ar.config.comment_string);
while(~isempty(C{1}))
    myFprintf(fidOut,'%s %f %n %n %f %f\n',C{1}{1},C{2},C{3},C{4},C{5},C{6});
    [C, fid] = arTextScan(fid, '%s %f %n %n %f %f\n',1, 'CommentStyle', ar.config.comment_string);
end

fclose(fidOut);
if ~isstruct( fid )
    fclose(fid);
end


function [ line, remainder, fid ] = readLine( fid, commentStyle )
line = '';
while ( isempty(line) )
    [line, fid] = arTextScan(fid, '%[^\n]' );
    line = strtrim(line{1}{1});
    Q = strfind( line, commentStyle );
    if ( ~isempty(Q) )
        line = line(1:Q-1);
    end
end
remainder = line;


function myFprintf(fidOut,format,varargin)
% fprintf(format,varargin{:})
format = strrep(strrep(strrep(format,'%q','%s'),'%n','%d'),' ','\t');
fprintf(fidOut,format,varargin{:});

function [str, remainder] = grabtoken( inputString, varargin )
    if ( isempty( inputString ) )
        str{1} = {};
        remainder = '';
        return;
    end

	[str, pos] = textscan(inputString, varargin{:});
	remainder = inputString(pos+1:end);    
