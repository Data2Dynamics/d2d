% Determines index position of each element of <A>
% in <B>. Similar to ismember but with index numbers.
%
%       pos_ia_in_b = position_of_ia_in_b(A, B)
%
% The length of the index vector will be 
% length(B). If non-unique matches exist, error.
% If elements are missing index will be zero used.

function pos_ia_in_b = position_of_ia_in_b(A, B)

pos_ia_in_b = zeros(size(B));

for j=1:length(B)
    jk = find(ismember(A, B{j}));
    if(length(jk)==1)
        pos_ia_in_b(j) = jk;
    elseif(length(jk)>1)
%         error('Multiple entries of B(%i) = %s in A\n', ...
%             j, B{j});
        warning('Multiple entries of B(%i) = %s in A, using first entry!\n', ...
            j, B{j});
        pos_ia_in_b(j) = jk(1);
    end
end