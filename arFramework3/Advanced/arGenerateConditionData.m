% name = arGenerateConditionData( name, [varargin] )
%
% Generate a dummy datafile for a specific condition
% 
%       name        Name of the dummy data
%       varargin    condition variables [none]

function name = arGenerateConditionData( name, varargin )
    fid = fopen( sprintf( 'Data/%s.csv', name ), 'w' );
    fprintf( fid, 't, dummy\n' );
    fprintf( fid, '0, 0\n');
    fclose(fid);
    
    copyfile(which('data_template.def'),['./Data/' name '.def']);
    fid = fopen( sprintf( 'Data/%s.def', name ), 'a' );
    for a = 1 : 2 : length( varargin )
        if ( isnumeric( varargin{a+1} ) )
            varargin{a+1} = num2str(varargin{a+1});
        end
        fprintf( fid, '\n%s "%s"', varargin{a}, varargin{a+1} );
    end
    fprintf( fid, '\n' );
    fclose(fid);
end

