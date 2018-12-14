% arSetPrior( par, type, mean, std )
%
% Function to set a prior on a parameter
%
%   par     parameter name or number (position in ar.p)
%   type    type of prior as number or string:
%           0 or 'bounds'
%           1 or 'normal'/'gaussian'/'quadratic'
%           2 or 'normal bounds'/'soft bounds'
%           3 or 'L1'/'lasso'
%   mean    mean value (center of prior)
%   std     std value of prior
%     
% Examples:
% 
% arSetPrior('kon', 2, 1.1, 0.2)
% 
% arSetPrior(7, 'soft bounds', -3, 0.2)
% 
% See also arPrint, arSetPars, arSetParsPattern

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