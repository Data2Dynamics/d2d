% a = arSym2str(b)
% 
% convert sym array to string array
% 
%   This function is used in arCompileAll (later in arCcode) and was called
%   sym2str previously.
% 
% This function is intended to catch all version dependencies when symbolic
% expressions are converted to strings.

function a = arSym2str(b)
    nonzeros = logical(b~=0);
    a = cell(size(b));
    for j=1:size(b,1)
        for jj = 1:size(b,2)
            if ( nonzeros(j,jj) )
                a{j,jj} = char(b(j,jj));
            else
                a{j,jj} = '0';
            end
        end
    end

    
%   The alternative function commented out below seems to be a version for
%   older matlab releases:
% 
%     % convert sym array to string array
% function a = sym2str_old(b)
%     a = cell(size(b));
%     for j=1:length(b)
%         a{j} = char(b(j));
%     end
%     
