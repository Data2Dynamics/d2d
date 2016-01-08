% determines index position of each element of <A>
% in <B>. The length of the index vector will be 
% length(B). If non-unique matches exist, error.
% If elements are missing index will be zero used.

function pos_ia_in_b = position_of_ia_in_b(A, B)

pos_ia_in_b = zeros(size(B));

for j=1:length(B)
    jk = find(ismember(A, B{j}));
    if(length(jk)==1)
        pos_ia_in_b(j) = jk;
    elseif(length(jk)>1)
        error('Multiple entries of B(%i) = %s in A\n', ...
            j, B{j});
    end
end