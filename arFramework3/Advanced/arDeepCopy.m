% copy = arDeepCopy(in)
%
% Recursively deepcopy a structure to avoid MATLAB's shallow copying behaviour
%
% in   - arStruct
% copy - arStruct
%
% Sometimes you really wish to compare subfields of two ar structures; one of which
% you are currently still simulating with. What can happen is that MATLAB makes a shallow 
% copy of some subfield that is filled by the mex simulation files. In this
% case, your shallow copy is also overwritten. For such cases, you can use
% this function to make a deep copy of the ar struct into another variable.
%
%   Deeply copies: substructures, cell arrays, matrices, numbers, strings and logicals
%   Does not deepcopy: matlab handles (such as figures)
%   This function is *very* slow, but will do the job
%
% Example
%   arModel = arDeepCopy(ar)

function copy = arDeepCopy(in)

	if ( isstruct( in ) )
        if ( length( in ) > 1 )
            for s = 1 : length( in )
                copy(s) = arDeepCopy(in(s));%#ok<AGROW>
            end
        else
            names = fieldnames(in);
            copy = struct;
            if ~isempty(names) && numel(in)>0
                for a = 1 : length( names )
                    copy.(names{a}) = arDeepCopy(in.(names{a}));
                end
            elseif ~isempty(names) && numel(in)==0
                copy = in;
            end
        end
    else
        if ( iscell( in ) )
            copy = cell(size(in));
            if ( ~isempty( in ) )
                for b = 1 : length( in )
                    copy{b} = arDeepCopy(in{b});
                end
            else
                copy = {};
            end
        elseif ( isnumeric( in ) )
            copy = in + 0;
        elseif ( islogical( in ) )
            copy = (in == 1);
        elseif ( ischar( in ) )
            copy = strcat( in, '' );
        elseif ( isstruct( in ) )
            copy = arDeepCopy(in);
        else
            % For all other datatypes, we do a shallow copy. This includes
            % things like figure handles and such.
            copy = in;
        end
    end
    
end