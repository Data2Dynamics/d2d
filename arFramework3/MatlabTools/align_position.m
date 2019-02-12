% data_aligned = align_position(data, dataLabels, refLabels, matching_dimensions)
%
% aligns <data> according to labels <dataLabels> and <refLabels>, if
% omitted, index positions will be returned
% such that <data_aligned> will have length(refLabels).
% For 2D data, looks for matching dimensions,
% unless <matching_dimensions> is specified.
% 
% Example:
%   align_position(1:4,{'a','A','a','b'},{'a','b','c','d','A','B'})
%   align_position(1:6,{'a','b','c','d','A','B'},{'a','A','a','b'})
% 
% Example: Calling the function without arguments profides an illustrative example:
%   align_position

function data_aligned = align_position(data, dataLabels, refLabels, matching_dimensions)

if(nargin==0)
    dataLabels = {'B' 'A' 'D' 'G'};
    refLabels = {'A' 'F' 'C' 'D' 'G' 'M' 'B'};
    data = rand(length(refLabels),length(dataLabels));
end

if(isempty(data))
    data = 1:length(dataLabels);
end

if(~isvector(dataLabels))
    error('dataLabels must be a vector');
end
if(~isvector(refLabels))
    error('refLabels must be a vector');
end

pos_ia_in_b = position_of_ia_in_b(dataLabels, refLabels);

data_size = size(data);
if(length(data_size)>2)
    error('function not supported in more than 2D data');
end

if(~exist('matching_dimensions','var'))
    matching_dimensions = find(data_size==length(dataLabels));
    if(isempty(matching_dimensions))
        error('one dimension of <data> must match length(dataLabels)');
    end
else
    if(size(data, matching_dimensions) ~= length(dataLabels))
        error('dimension %i of <data> must match length(dataLabels)', matching_dimensions);
    end
end
if(length(matching_dimensions)>1)
    error('multiple matching dimensions, please specify matching_dimensions');
end
align_dimension = matching_dimensions(1);

data_size_new = data_size;
data_size_new(align_dimension) = length(refLabels);

if(iscell(data))
    data_aligned = repmat({'n/a'},data_size_new);
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
    disp('Example illustrating the function:')
    disp('data=')
    disp(data)
    disp('dataLabels=')
    disp(dataLabels)
    disp('refLabels=')
    disp(refLabels)
    disp('data_aligned=')
    disp(data_aligned)
end





