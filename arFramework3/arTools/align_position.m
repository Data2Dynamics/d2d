% aligns <data> according to labels <A> and <B>
%
%    data_aligned = align_position(data, A, B)
%
% such that <data_aligned> will have length(B).
% For 2D data, looks for matching dimensions.
% If the length of the data dimensions are equal, 
% uses first matching dimension.

function data_aligned = align_position(data, A, B)

if(nargin==0)
    A = {'B' 'A' 'D' 'G'};
    B = {'A' 'F' 'C' 'D' 'G' 'M' 'B'};
    data = rand(length(B),length(A));
end

if(~isvector(A))
    error('A must be a vector');
end
if(~isvector(B))
    error('B must be a vector');
end

pos_ia_in_b = position_of_ia_in_b(A, B);

data_size = size(data);
if(length(data_size)>2)
    error('function not supported in more than 2D data');
end

matching_dimensions = find(data_size==length(A));
if(isempty(matching_dimensions))
    error('one dimension of <data> must match length(A)');
end
align_dimension = matching_dimensions(1);

data_size_new = data_size;
data_size_new(align_dimension) = length(B);

if(iscell(data))
    data_aligned = cell(data_size_new);
elseif(isnumeric(data))
    data_aligned = nan(data_size_new);
else
    error('unsupported data type %s', class(data));
end

% do align
for j=1:length(pos_ia_in_b)
    if(pos_ia_in_b(j)~=0)
        if(matching_dimensions==1)
            data_aligned(j,:) = data(pos_ia_in_b(j),:);
        elseif(matching_dimensions==2)
            data_aligned(:,j) = data(:,pos_ia_in_b(j));
        end
    end
end

if(nargin==0)
    disp(A)
    disp(B)
    disp(data)
    disp(data_aligned)
end





