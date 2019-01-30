% str = arAssembleReactionStr(source, target, [simbio, sourceCoeffs, targetCoeffs])
%
% Transform input strings to reaction equation string ('source -> target')
%
% source ['']   cell array of reactants, eg {'A','B'}
% target ['']   cell array of products, eg {'C'}
% simbio [false] if true, replace '' to 'null'
% sourceCoeffs  vector of source coefficients
% targetCoeffs  vector of target coefficients / how many molecules
% 
% 
% Example:
%  arAssembleReactionStr({'A','B'},{'C'}) results in 'A + B -> C'
%  arAssembleReactionStr({'A','B'},{'C'},false,[1 2]) results in 'A + 2 B -> C'

function str = arAssembleReactionStr(source, target, simbio, sourceCoeffs, targetCoeffs)

    if(~exist('simbio', 'var'))
        simbio = false;
    end

    if(isempty(source))
        s = '';
    else
        if ( nargin > 3 )
            s = [ coeff(sourceCoeffs(1)) source{1} ];
            for jj=2:length(source)
                s = [s ' + ' coeff(sourceCoeffs(jj)) source{jj}]; %#ok<AGROW>
            end
        else
            s = source{1};
            for jj=2:length(source)
                s = [s ' + ' source{jj}]; %#ok<AGROW>
            end    
        end
    end

    if(isempty(target))
        t = '';
    else
        if ( nargin > 4 )
            t = [ coeff(targetCoeffs(1)) target{1} ];
            for jj=2:length(target)
                t = [t ' + ' coeff(targetCoeffs(jj)) target{jj}]; %#ok<AGROW>
            end
        else
            t = target{1};
            for jj=2:length(target)
                t = [t ' + ' target{jj}]; %#ok<AGROW>
            end    
        end
    end    
    if(simbio)
        if(isempty(s))
            s = 'null';
        end
        if(isempty(t))
            t = 'null';
        end
    end

    str = sprintf('%s -> %s', s, t);

function str = coeff( in )

    if ( in ~= 1 )
        str = [num2str(in) ' '];
    else
        str = '';
    end
