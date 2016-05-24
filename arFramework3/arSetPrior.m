% Function to set a prior on a parameter
%
% arSetPrior( par, type of prior, mean, std )
%   par     - parameter or name/mask
%   type    - bounds/normal/soft bounds/L1
%   mean    - mean
%   std     - std

function arSetPrior( par, type, mean, std )
    global ar;

    if ischar(par)
        par = arFindPar(par);
    end
    
    if ischar(type)
        switch(lower(type))
            case 'bounds'
                type = 0;
            case 'normal'
                type = 1;
            case 'gaussian'
                type = 1;                
            case 'quadratic'
                type = 1; 
            case 'normal bounds'
                type = 2;
            case 'soft bounds'
                type = 2;
            case 'L1'
                type = 3;
            case 'lasso'
                type = 3;                
        end
    end
    
    if ~isnumeric( par )
        error( 'Did not provide valid parameter ID or name' );
    end
    
    for a = 1 : length( par )
        ar.type(par(a)) = type;
        ar.mean(par(a)) = mean;
        ar.std(par(a)) = std;
    end
end